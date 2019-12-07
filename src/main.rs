/*
Input:
 - n_filt: filter length
 - j_type: type of filter
    - 1 = multiple passband / stopband filter
    - 2 = differentiator
    - 3 = hilbert transform filter
 - n_bands: number of bands
 - l_grid: grid density (will be set to 16 unless specified otherwise by a positive constant)
 - edge (2 * n_bands): band-edge array: lower and upper edges for each band (max 10 bands)
 - fx (n_bands): desired function array (or desired slope if differentiator) for each band
 - wtx (n_bands): weight function array in each band. For a differentiator, the weight function is
                  inversely proportional to f.
*/

mod tests;

const PI: f64 = std::f64::consts::PI;
const PI2: f64 = PI * 2.0;

// The program is set up for a maximum length of 128, but this upper limit can be changed by
// redimensioning the arrays IEXT, AD, ALPHA, X, Y, and H to be NF_MAX/2 + 2.
// The DES, GRID, and WT arrays must be dimensioned 16(NF_MAX / 2 + 2).

pub struct Band {
    pub lower_edge: f32,
    pub upper_edge: f32,
    pub desired_value: f32,
    pub weight: f32,
}

const NF_MAX: usize = 128;

#[derive(Copy, Clone)]
#[allow(dead_code)]
enum JType {
    MultipleBand,   // 1
    Differentiator, // 2
    Hilbert,        // 3
}

// This output captures all of the output printed by the Fortran code from the whitepaper.
// The format isn't great, so TODO: refactor.
#[derive(Clone)]
pub struct ParksMcclellanOutput {
    pub filter_length: usize,
    pub impulse_response: Vec<f32>,
    pub lower_band_edges: Vec<f32>,
    pub upper_band_edges: Vec<f32>,
    pub desired_values: Vec<f32>,
    pub weightings: Vec<f32>,
    pub deviations: Vec<f32>,
    pub deviation_dbs: Vec<f32>,
    pub extremal_frequencies: Vec<f32>,
}

// Note: l_grid is defaulted to 16 in the Parks-McClellan paper.
fn design(n_filt: usize, j_type: JType, bands: &Vec<Band>, l_grid: usize) -> ParksMcclellanOutput {
    // The filter length must be greater than zero, and less/equal to NF_MAX.
    assert!(n_filt >= 3);
    assert!(n_filt <= NF_MAX);

    // n_bands must be greater than 0.
    let n_bands = bands.len();
    assert!(n_bands > 0);

    // l_grid must be greater than 0.
    assert!(l_grid > 0);

    // The algorithm uses these internally ...
    let mut edges = vec![0.0f32; n_bands * 2];
    let mut fx = vec![0.0f32; n_bands];
    let mut wtx = vec![0.0f32; n_bands];
    // ... so let's fill it up with our `bands` input. TODO: refactor.
    for (i, band) in bands.iter().enumerate() {
        edges[i*2] = band.lower_edge;
        edges[i*2+1] = band.upper_edge;
        fx[i] = band.desired_value;
        wtx[i] = band.weight;
    }

    let neg = match j_type {
        JType::MultipleBand => 0,
        JType::Differentiator => 1,
        JType::Hilbert => 1,
    };

    let mut n_odd = n_filt / 2;
    n_odd = n_filt - 2 * n_odd;

    let mut nfcns = n_filt / 2;

    if n_odd == 1 && neg == 0 {
        nfcns += 1;
    }

    // Set up the dense grid. The number of points in the grid is
    // (filter_length + 1) * (grid_density) / 2.
    let mut grid = [0.0f32; 1045];
    grid[0] = edges[0];
    let mut del_f = (l_grid as f32) * (nfcns as f32);
    del_f = 0.5 / del_f;

    if neg != 0 {
        if edges[0] < del_f {
            grid[0] = del_f;
        }
    }

    let mut des = [0.0f32; 1045];
    let mut wt = [0.0f32; 1045];

    let mut j = 1;
    let mut L = 1;
    let mut l_band = 1;

    'b: loop {
        let f_up = edges[L]; // edges[L+1-1]
        let mut temp: f32;
        'a: loop {
            temp = grid[j-1];
            // Calculate the desired magnitude response,
            // and the weight function on the grid.
            des[j-1] = eff(temp, &fx, &wtx, l_band, j_type);
            wt[j-1] = wate(temp, &fx, &wtx, l_band, j_type);
            j += 1;
            grid[j-1] = temp + del_f;
            if grid[j-1] > f_up {
                break 'a;
            }
        }

        grid[j-2] = f_up; // grid[j-1-1]
        des[j-2] = eff(f_up, &fx, &wtx, l_band, j_type);
        wt[j-2] = wate(f_up, &fx, &wtx, l_band, j_type);
        l_band = l_band + 1;
        L += 2;
        if l_band > n_bands {
            break 'b;
        }
        grid[j-1] = edges[L-1];
    }

    let mut n_grid = j-1;

    if neg == n_odd {
        if grid[n_grid-1] > (0.5 - del_f) {
            n_grid = n_grid - 1;
        }
    }

    // Set up a new approximation problem which is equivalent to the original problem.
    if neg == 0 {
        if n_odd == 1 {
            // NOP (?)
        } else {
            for j in 1..(n_grid+1) {
                let change: f32 = (PI * grid[j-1] as f64).cos() as f32;
                des[j-1] = des[j-1] / change;
                wt[j-1] = wt[j-1] * change;
            }
        }
    }
    else {
        if n_odd == 1 {
            for j in 1..(n_grid+1) {
                let change = (PI2 * grid[j-1] as f64).sin() as f32;
                des[j-1] = des[j-1] / change;
                wt[j-1] = wt[j-1] * change;
            }
        } else {
            for j in 1..(n_grid+1) {
                let change = (PI * grid[j-1] as f64).sin() as f32;
                des[j-1] = des[j-1] / change;
                wt[j-1] = wt[j-1] * change;
            }
        }
    }

    // Initial guess for the extremal frequencies: equally spaced along the grid.
    let mut iext = [0; 66];
    let temp = ((n_grid-1) as f32) / (nfcns as f32);
    for j in 1..(nfcns+1) {
        let xt = j-1;
        iext[j-1] = ((xt as f32) * temp + 1.0) as i64;
    }
    iext[nfcns] = n_grid as i64; // iext[nfcns+1-1]
    let nm1 = nfcns - 1;
    let nz = nfcns + 1;

    // Call the remez exchange algorithm to do the approximation problem.
    let mut alpha = [0.0f32; 66];
    let mut dev: f64 = 0.0;
    remez(nfcns, &mut iext, n_grid, &mut grid, &wt, &mut des, &mut alpha, &mut dev);

    // Calculate the impulse response.
    let mut h = [0.0f32; 66];
    if neg == 0 {
        if n_odd == 0 {
            h[0] = 0.25 * alpha[nfcns-1];
            for j in 2..(nm1+1) {
                let nzmj = nz-j;
                let nf2j = nfcns+2-j;
                h[j-1] = 0.25 * (alpha[nzmj-1] + alpha[nf2j-1]);
            }
            h[nfcns-1] = 0.5 * alpha[0] + 0.25 * alpha[1];
        } else {
            for j in 1..(nm1+1) {
                let nzmj = nz-j;
                h[j-1] = 0.5 * alpha[nzmj-1];
            }
            h[nfcns-1] = alpha[0];
        }
    } else {
        if n_odd == 0 {
            h[0] = 0.25 * alpha[nfcns-1];
            for j in 2..(nm1+1) {
                let nzmj = nz-j;
                let nf2j = nfcns+2-j;
                h[j-1] = 0.25 * (alpha[nzmj-1] - alpha[nf2j-1]);
            }
            h[nfcns-1] = 0.5 * alpha[0] - 0.25 * alpha[1];
        } else {
            h[0] = 0.25 * alpha[nfcns-1];
            if nm1 > 0 {
                // Fortran treats indexing to the "zeroth" element
                // as a no-op, I guess? (remember it's 1-indexed)
                h[1] = 0.25 * alpha[nm1-1];
            }
            for j in 3..(nm1+1) {
                let nzmj = nz-j;
                let nf3j = nfcns+3-j;
                h[j-1] = 0.25 * (alpha[nzmj-1] - alpha[nf3j-1]);
            }
            h[nfcns-1] = 0.5 * alpha[0] - 0.25 * alpha[2];
            h[nz-1] = 0.0;
        }
    }

    // Program output section

    let mut parks_mcclellan_output = ParksMcclellanOutput {
        filter_length: 0,
        impulse_response: vec![],
        lower_band_edges: vec![],
        upper_band_edges: vec![],
        desired_values: vec![],
        weightings: vec![],
        deviations: vec![],
        deviation_dbs: vec![],
        extremal_frequencies: vec![],
    };

    match j_type {
        JType::MultipleBand => println!("BANDPASS FILTER"),
        JType::Differentiator => println!("DIFFERENTIATOR"),
        JType::Hilbert => println!("HILBERT TRANSFORMER"),
    }
    println!("Filter length: {}", n_filt);
    parks_mcclellan_output.filter_length = n_filt;
    println!();
    println!("***** IMPULSE RESPONSE *****");
    for j in 1..(nfcns+1) {
        println!("{:e}", h[j-1]);
        parks_mcclellan_output.impulse_response.push(h[j-1]);
    }
    println!();
    if neg == 1 && n_odd == 1 {
        parks_mcclellan_output.impulse_response.push(0.0);
    }

    for k in 0..n_bands {
        println!("Band {}:", k);
        println!("    lower edge: {}", edges[2*k]);
        parks_mcclellan_output.lower_band_edges.push(edges[2*k]);
        println!("    upper edge: {}", edges[2*k+1]);
        parks_mcclellan_output.upper_band_edges.push(edges[2*k+1]);
        match j_type {
            JType::MultipleBand => {
                println!("    desired value: {}", fx[k]);
                parks_mcclellan_output.desired_values.push(fx[k]);
            },
            JType::Differentiator => {
                println!("    desired slope: {}", fx[k]);
                parks_mcclellan_output.desired_values.push(fx[k]);
            },
            JType::Hilbert => {
                println!("    desired value: {}", fx[k]);
                parks_mcclellan_output.desired_values.push(fx[k]);
            },
        }
        println!("    weighting: {}", wtx[k]);
        parks_mcclellan_output.weightings.push(wtx[k]);

        let deviation = (dev/(wtx[k] as f64)) as f32;
        println!("    deviation: {}", deviation);
        parks_mcclellan_output.deviations.push(deviation);
        match j_type {
            JType::MultipleBand => {
                let deviation_db: f32 = 20.0 * (deviation + fx[k]).log10();
                println!("    deviation in dB: {}", deviation_db);
                parks_mcclellan_output.deviation_dbs.push(deviation_db);
            },
            _ => {}
        }
    }

    println!("Extremal frequencies (maxima of the error curve)");
    for j in 1..(nz+1) {
        let ix = iext[j-1];
        grid[j-1] = grid[ix as usize - 1];
        println!("{}", grid[j-1]);
        parks_mcclellan_output.extremal_frequencies.push(grid[j-1]);
    }

    parks_mcclellan_output
}

// Function to calculate the desired magnitude response as a function of frequency.
// An arbitrary function of frequency can be approximated if the user replaces this function
// with the appropriate code to evaluate the ideal magnitude.
// Note that the parameter `freq` is the value of **normalized** frequency needed for evaluation.
fn eff(freq: f32, fx: &Vec<f32>, _wtx: &Vec<f32>, l_band: usize, j_type: JType) -> f32 {
    match j_type {
        JType::Differentiator => {
            fx[l_band-1] * freq
        },
        _ => {
            fx[l_band-1]
        },
    }
}

// Function to calculate the weight function as a function of frequency.
// Similar to the function `eff`, this function can be replaced by a user-written function
// to calculate any desired weighting function.
fn wate(freq: f32, fx: &Vec<f32>, wtx: &Vec<f32>, l_band: usize, j_type: JType) -> f32 {
    match j_type {
        JType::Differentiator => {
            if fx[l_band-1] < 0.0001 {
                wtx[l_band-1]
            } else {
                wtx[l_band-1] / freq
            }
        },
        _ => {
            wtx[l_band-1]
        }
    }
}

// This function implements the Remez Exchange Algorithm for the weighted Chebyshev approximation
// of a continuous function with a sum of cosines.
// Inputs to the subroutine are:
//  - a dense grid which replaces the frequency axis
//  - the desired function on this grid
//  - the number of cosines
//  - an initial guess of the extremal frequencies
// The program minimizes the Chebyshev error by determining the best location of the extremal
// frequencies (points of maximum error) and then calculates the coefficients of the best
// approximation.
// This function is a giant mess of GOTO spaghetti in the Fortran code,
// so I'm just going to model it as a state machine. I'll clean it up later.
fn remez(
    nfcns: usize,
    iext: &mut [i64; 66],
    n_grid: usize,
    grid: &mut [f32; 1045],
    wt: &[f32; 1045],
    des: &[f32; 1045],
    alpha: &mut [f32; 66],
    dev: &mut f64,
) {
    // Function-scoped data
    let mut x = [0.0f64; 66];
    let mut y = [0.0f64; 66];
    let mut ad = [0.0; 66];
    let mut nu = 0;
    let mut jchnge = 0;
    let mut k1 = 0;
    let mut knz = 0;
    let mut klow = 0;
    let mut nut = 0;
    let mut non_loop_j = 0;
    let mut kup = 0;
    let mut L = 0;
    let mut err = 0.0f32;
    let mut nut1 = 0;
    let mut luck = 0;
    let mut aa = 0.0f32;
    let mut bb = 0.0f32;

    let mut ynz = None;
    let mut comp = None;
    let mut y1 = None;

    let itrmax = 25; // Max number of iterations
    let mut devl = -1.0f32;
    let nz = nfcns + 1;
    let nzz = nfcns + 2;
    let mut niter = 0;

    let mut state = 100;
    loop {
        match state {
            100 => {
                iext[nzz-1] = n_grid as i64 + 1;
                niter += 1;
                println!("niter: {}", niter);
                if niter > itrmax {
                    state = 400; continue; // GOTO 400
                }
                for j in 1..(nz+1) {
                    let jxt = iext[j-1];
                    let mut dtemp: f64 = grid[(jxt-1) as usize] as f64;
                    dtemp = (dtemp * PI2).cos();
                    x[j-1] = dtemp;
                }
                let jet = (nfcns-1) / 15 + 1;
                for j in 1..(nz+1) {
                    ad[j-1] = d_func(&x, j, nz, jet);
                }

                let mut dnum = 0.0;
                let mut dden = 0.0;
                let mut k = 1;
                for j in 1..(nz+1) {

                    L = iext[j-1];
                    let dtemp = ad[j-1] * des[(L-1) as usize] as f64;
//                    println!("L: {}", L);
//                    println!(" ADJ = {}", ad[j-1]);
//                    println!("DESL = {}", des[(L-1) as usize]);
//                    println!(" WTL = {}", wt[(L-1) as usize]);
                    dnum += dtemp;
                    let dtemp = (k as f64) * ad[j-1] / wt[(L-1) as usize] as f64;
                    dden += dtemp;
//                    println!("DNUM = {}", dnum);
//                    println!("DDEN = {}", dden);
                    k = -k;
                }
                *dev = dnum / dden;
                println!("DEVIATION: {}", *dev);

                nu = 1;
                if *dev > 0.0 {
                    nu = -1;
                }
                *dev = -(nu as f64) * *dev;
                k = nu;

                for j in 1..(nz+1) {
                    L = iext[j-1];
                    let dtemp = (k as f64) * *dev / wt[(L-1) as usize] as f64;
                    y[j-1] = des[(L-1) as usize] as f64 + dtemp;
                    k = -k;
                }

                if *dev <= devl as f64 {
                    println!("***** FAILURE TO CONVERGE *****");
                    println!("Number of iterations: {}", niter);
                    println!("If the number of iterations is greater than 3,");
                    println!("the design might be correct, but should be verified by FFT.");
                    state = 400; continue; // GOTO 400
                }

                devl = *dev as f32;
                jchnge = 0;
                k1 = iext[0];
                knz = iext[nz-1];
                klow = 0;
                nut = -nu;
                non_loop_j = 1;
                state = 200; continue; // fall through to 200

            },

            // Search for the extremal frequencies of the best approximation
            200 => {
                if non_loop_j == nzz {
                    ynz = comp;
                }
                if non_loop_j >= nzz {
                    state = 300; continue; // GOTO 300
                }
                kup = iext[non_loop_j]; // iext[j+1-1]
                L = iext[non_loop_j-1] + 1;
                nut = -nut;
                if non_loop_j == 2 {
                    y1 = comp;
                }
                comp = Some(*dev);
                if L >= kup {
                    state = 220; continue; // GOTO 220
                }
                err = gee(grid, &x, &y, &ad, L, nz) as f32;
                err = (err - des[(L-1) as usize]) * wt[(L-1) as usize];
                let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                if dtemp <= 0.0 {
                    state = 220; continue; // GOTO 220
                }
                comp = Some((nut as f64) * (err as f64));
                state = 210; continue; // fall through to 210
            },
            210 => {
                L = L + 1;
                if L >= kup {
                    state = 215; continue; // GOTO 215
                }
                err = gee(grid, &x, &y, &ad, L, nz) as f32;
                err = (err - des[(L-1) as usize]) * wt[(L-1) as usize];
                let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                if dtemp <= 0.0 {
                    state = 215; continue; // GOTO 215
                }
                comp = Some((nut as f64) * (err as f64));
                state = 210; continue; // GOTO 210 (self loop)
            },
            215 => {
                iext[non_loop_j-1] = L - 1;
                non_loop_j = non_loop_j + 1;
                klow = L - 1;
                jchnge = jchnge + 1;
                state = 200; continue; // GOTO 200
            },
            220 => {
                L = L - 1;
                state = 225; continue; // fall through to 225
            },
            225 => {
                L = L - 1;
                if L <= klow {
                    state = 250; continue; // GOTO 250
                }
                err = gee(grid, &x, &y, &ad, L, nz) as f32;
                err = (err - des[(L-1) as usize]) * wt[(L-1) as usize];
                let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                if dtemp > 0.0 {
                    state = 230; continue; // GOTO 230
                }
                if jchnge <= 0 {
                    state = 225; continue; // GOTO 225 (self loop)
                }
                state = 260; continue; // GOTO 260
            },
            230 => {
                comp = Some((nut as f64) * (err as f64));
                state = 235; continue; // fall through to 235
            },
            235 => {
                L = L - 1;
                if L <= klow {
                    state = 240; continue; // GOTO 240
                }
                err = gee(grid, &x, &y, &ad, L, nz) as f32;
                err = (err - des[(L-1) as usize]) * wt[(L-1) as usize];
                let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                if dtemp <= 0.0 {
                    state = 240; continue; // GOTO 240
                }
                comp = Some((nut as f64) * (err as f64));
                state = 235; continue; // GOTO 235 (self loop)
            },
            240 => {
                klow = iext[non_loop_j-1];
                iext[non_loop_j-1] = L+1;
                non_loop_j += 1;
                jchnge += 1;
                state = 200; continue; // GOTO 200
            },
            250 => {
                L = iext[non_loop_j-1] + 1;
                if jchnge > 0 {
                    state = 215; continue; // GOTO 215
                }
                state = 255; continue; // fall through to 255
            },
            255 => {
                L += 1;
                if L >= kup {
                    state = 260; continue; // GOTO 260
                }
                err = gee(grid, &x, &y, &ad, L, nz) as f32;
                err = (err - des[(L-1) as usize]) * wt[(L-1) as usize];
                let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                if dtemp <= 0.0 {
                    state = 255; continue; // GOTO 255 (self loop)
                }
                comp = Some((nut as f64) * (err as f64));
                state = 210; continue; // GOTO 210
            },
            260 => {
                klow = iext[non_loop_j-1];
                non_loop_j += 1;
                state = 200; continue; // GOTO 200
            },
            300 => {
                if non_loop_j > nzz {
                    state = 320; continue; // GOTO 320
                }
                if k1 > iext[0] {
                    k1 = iext[0];
                }
                if knz < iext[nz-1] {
                    knz = iext[nz-1];
                }
                nut1 = nut;
                nut = -nu;
                L = 0;
                kup = k1;
                comp = Some(ynz.unwrap() * 1.00001);
                luck = 1;
                state = 310; continue; // fall through to 310
            },
            310 => {
                L = L+1;
                if L >= kup {
                    state = 315; continue; // GOTO 315
                }
                err = gee(grid, &x, &y, &ad, L, nz) as f32;
                err = (err - des[(L-1) as usize]) * wt[(L-1) as usize];
                let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                if dtemp <= 0.0 {
                    state = 310; continue; // GOTO 310 (self loop)
                }
                comp = Some((nut as f64) * (err as f64));
                non_loop_j = nzz;
                state = 210; continue; // GOTO 210
            },
            315 => {
                luck = 6;
                state = 325; continue; // GOTO 325
            },
            320 => {
                if luck > 9 {
                    state = 350; continue; // GOTO 350
                }
                if comp.unwrap() > y1.unwrap() {
                    y1 = comp;
                }
                k1 = iext[nzz-1];
                state = 325; continue; // fall through to 325
            },
            325 => {
                L = n_grid as i64 + 1;
                klow = knz;
                nut = -nut1;
                comp = Some(y1.unwrap() * 1.00001);
                state = 330; continue; // fall through to 330
            },
            330 => {
                L = L-1;
                if L <= klow {
                    state = 340; continue; // GOTO 340
                }
                err = gee(grid, &x, &y, &ad, L, nz) as f32;
                err = (err - des[(L-1) as usize]) * wt[(L-1) as usize];
                let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                if dtemp <= 0.0 {
                    state = 330; continue; // GOTO 330 (self loop)
                }
                non_loop_j = nzz;
                comp = Some((nut as f64) * (err as f64));
                luck += 10;
                state = 235; continue; // GOTO 235
            },
            340 => {
                if luck == 6 {
                    state = 370; continue; // GOTO 370
                }
                for j in 1..(nfcns+1) {
                    let nzzmj = nzz - j;
                    let nzmj = nz - j;
                    iext[nzzmj-1] = iext[nzmj-1];
                }
                iext[0] = k1;
                state = 100; continue; // GOTO 100 (next iteration!)
            },
            350 => {
                let kn = iext[nzz-1];
                for j in 1..(nfcns+1) {
                    iext[j-1] = iext[j] // j+1-1
                }
                iext[nz-1] = kn;
                state = 100; continue; // GOTO 100 (next iteration!)
            },
            370 => {
                if jchnge > 0 {
                    state = 100; continue; // GOTO 100 (next iteration!)
                }
                state = 400; continue; // fall through to 400
            },

            // Calculation of the coefficients of the best approximation using the inverse DFT
            // By here, we're done iterating the Remez exchange algorithm, and we're mostly just
            // calculating "outputs".
            400 => {
                let mut a = [0.0f64; 66];

                let nm1 = nfcns - 1;
                let fsh: f32 = 1.0e-06;
                let gtemp: f32 = grid[0]; // grid[1-1]
                x[nzz-1] = -2.0;
                let cn = 2 * nfcns - 1;
                let delf = 1.0f32 / (cn as f32);
                L = 1;
                let mut kkk = 0;
                if grid[0] < 0.01 && grid[n_grid-1] > 0.49 {
                    kkk = 1;
                }
                if nfcns <= 3 {
                    kkk = 1;
                }

                if kkk != 1 {
                    let dtemp = (PI2 * grid[0] as f64).cos();
                    let dnum = (PI2 * grid[n_grid - 1] as f64).cos();
                    aa = (2.0 / (dtemp - dnum)) as f32;
                    bb = (-(dtemp + dnum) / (dtemp - dnum)) as f32;
                }

                // jump label 405 + 1
                for j in 1..(nfcns+1) {
                    let mut ft = (j-1) as f32;
                    ft = ft * delf;
                    let mut xt: f32 = (PI2 * ft as f64).cos() as f32;
                    if kkk != 1 {
                        xt = (xt - bb) / aa;
                        let xt1 = (1.0 - xt.powi(2)).sqrt();
                        ft = (xt1.atan2(xt) as f64 / PI2) as f32;
                    }

                    loop {
                        let xe: f32 = x[(L-1) as usize] as f32;
                        if xt > xe {
                            if (xt-xe) < fsh {
                                a[j-1] = y[(L-1) as usize];
                                break; // Stop looping
                            } else {
                                grid[0] = ft;
                                a[j - 1] = gee(grid, &x, &y, &ad, 1, nz);
                                break; // Stop looping
                            }
                        } else if (xe-xt) < fsh {
                            a[j-1] = y[(L-1) as usize];
                            break; // Stop looping
                        } else {
                            L = L+1;
                            continue; // Loop again
                        }
                    }

                    if L > 1 {
                        L -= 1;
                    }
                }

                // jump label 430
                grid[0] = gtemp;
                let dden = PI2 / (cn as f64);

                for j in 1..(nfcns+1) {
                    let mut dtemp = 0.0;
                    let mut dnum = (j-1) as f64;
                    dnum = dnum * dden;
                    if nm1 >= 1 {
                        for k in 1..(nm1+1) {
                            let dak = a[k]; // a[k+1-1]
                            let dk = k as f64;
                            dtemp = dtemp + dak * (dnum*dk).cos();
                        }
                    }
                    dtemp = 2.0 * dtemp + a[0];
                    alpha[j-1] = dtemp as f32;
                }

                for j in 2..(nfcns+1) {
                    alpha[j-1] = 2.0 * alpha[j-1] / (cn as f32);
                }
                alpha[0] = alpha[0] / (cn as f32);

                if kkk == 1 {
                    if nfcns > 3 {
                        return;
                    }
                    alpha[nfcns] = 0.0; // alpha[nfcns+1-1]
                    alpha[nfcns+1] = 0.0; // alpha[nfcns+2-1]
                    return;
                }

                let mut p = [0.0f64; 65];
                let mut q = [0.0f64; 65];

                p[0] = (2.0 * alpha[nfcns-1] * bb + alpha[nm1-1]) as f64;
                p[1] = (2.0 * aa * alpha[nfcns-1]) as f64;
                q[0] = (alpha[nfcns-3] - alpha[nfcns-1]) as f64;

                for j in 2..(nm1+1) {
                    if j >= nm1 {
                        aa = 0.5 * aa;
                        bb = 0.5 * bb;
                    }
                    p[j] = 0.0; // p[j+1-1]
                    for k in 1..(j + 1) {
                        a[k-1] = p[k-1];
                        p[k-1] = 2.0 * (bb as f64) * a[k - 1];
                    }
                    p[1] = p[1] + a[0] * 2.0 * (aa as f64);
                    let jm1 = j - 1;
                    for k in 1..(jm1 + 1) {
                        p[k - 1] = p[k - 1] + q[k - 1] + (aa as f64) * a[k]; // a[k+1-1]
                    }
                    let jp1 = j + 1;
                    for k in 3..(jp1 + 1) {
                        p[k - 1] = p[k - 1] + (aa as f64) * a[k - 2] // a[k-1-1]
                    }

                    if j != nm1 {
                        for k in 1..(j + 1) {
                            q[k - 1] = -a[k - 1];
                        }
                        let nf1j = nfcns - 1 - j;
                        q[0] = q[0] + (alpha[nf1j - 1] as f64);
                    }
                }

                for j in 1..(nfcns+1) {
                    alpha[j-1] = (p[j-1] as f32);
                }

                if nfcns > 3 {
                    return;
                }
                alpha[nfcns] = 0.0; // alpha[nfcns+1-1]
                alpha[nfcns+1] = 0.0; // alpha[nfcns+2-1]
                return;
            },
            _ => panic!("Unknown state {}", state),
        }
    }
}

// Function to calculate the lagrange interpolation coefficients
// for use in the function `gee`
fn d_func(x: &[f64; 66], k: usize, n: usize, m: usize) -> f64 {
    let mut d = 1.0;
    let q = x[k-1];

    let mut l = 1;
    while l <= m {
        let mut j = l;
        while j <= n {
            if j != k {
                d = 2.0 * d * (q - x[j-1]);
            }
            j += m;
        }
        l += 1;
    }
    d = 1.0 / d;
    d
}

// Function to evaluate the frequency response using the lagrange interpolation formula
// in the barycentric form
fn gee(grid: &[f32; 1045], x: &[f64; 66], y: &[f64; 66], ad: &[f64; 66], k: i64, n: usize) -> f64 {
    let mut p = 0.0;
    let mut d = 0.0;

    let mut xf = (grid[(k-1) as usize] as f64);
    xf = (PI2 * xf).cos();
    for j in 1..(n+1) {
        let mut c = xf - x[j-1];
        c = ad[j-1] / c;
        d += c;
        p += c*y[j-1];
    }

    p/d
}

fn main() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design(10, JType::MultipleBand, &bands, 16);
}