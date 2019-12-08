#[cfg(test)] mod tests;

pub mod dense_grid;
use dense_grid::DenseGrid;

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
pub enum JType {
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
    // The filter length must be [3, NF_MAX].
    assert!(n_filt >= 3);
    assert!(n_filt <= NF_MAX);

    // Must have at least one band
    assert!(bands.len() > 0);

    // l_grid must be greater than 0.
    assert!(l_grid > 0);

    let neg = match j_type {
        JType::MultipleBand => false,
        JType::Differentiator => true,
        JType::Hilbert => true,
    };

    let odd = n_filt % 2 != 0;

    // The number of unique coefficients in the resulting filter
    // (half will be reflected around the midpoint, but we need to include the center for odd N)
    let mut num_coefficients = n_filt / 2;
    if odd && !neg {
        num_coefficients += 1;
    }
    let num_coefficients = num_coefficients;

    // Set up the dense grid.
    let mut grid = DenseGrid::new(&bands, j_type, l_grid, num_coefficients, neg, odd);
    grid.reformat_grid_for_remez(neg, odd);
    let grid = grid;


    // Initial guess for the extremal frequencies: equally spaced along the grid.
    let mut iext = [0; 66];
    let temp = ((grid.n_grid()-1) as f32) / (num_coefficients as f32);
    for j in 0..num_coefficients {
        iext[j] = ((j as f32) * temp + 1.0) as i64;
    }
    iext[num_coefficients] = grid.n_grid() as i64; // iext[nfcns+1-1]
    let nm1 = num_coefficients - 1;
    let nz = num_coefficients + 1;

    // Call the remez exchange algorithm to do the approximation problem.
    let mut alpha = [0.0f32; 66];
    let mut dev: f64 = 0.0;
    remez(num_coefficients, &grid, &mut iext, &mut alpha, &mut dev);

    // Calculate the impulse response.
    let mut h = [0.0f32; 66];
    if !neg {
        if !odd {
            h[0] = 0.25 * alpha[num_coefficients -1];
            for j in 2..(nm1+1) {
                let nzmj = nz-j;
                let nf2j = num_coefficients +2-j;
                h[j-1] = 0.25 * (alpha[nzmj-1] + alpha[nf2j-1]);
            }
            h[num_coefficients -1] = 0.5 * alpha[0] + 0.25 * alpha[1];
        } else {
            for j in 1..(nm1+1) {
                let nzmj = nz-j;
                h[j-1] = 0.5 * alpha[nzmj-1];
            }
            h[num_coefficients -1] = alpha[0];
        }
    } else {
        if !odd {
            h[0] = 0.25 * alpha[num_coefficients -1];
            for j in 2..(nm1+1) {
                let nzmj = nz-j;
                let nf2j = num_coefficients +2-j;
                h[j-1] = 0.25 * (alpha[nzmj-1] - alpha[nf2j-1]);
            }
            h[num_coefficients -1] = 0.5 * alpha[0] - 0.25 * alpha[1];
        } else {
            h[0] = 0.25 * alpha[num_coefficients -1];
            if nm1 > 0 {
                // Fortran treats indexing to the "zeroth" element
                // as a no-op, I guess? (remember it's 1-indexed)
                h[1] = 0.25 * alpha[nm1-1];
            }
            for j in 3..(nm1+1) {
                let nzmj = nz-j;
                let nf3j = num_coefficients +3-j;
                h[j-1] = 0.25 * (alpha[nzmj-1] - alpha[nf3j-1]);
            }
            h[num_coefficients -1] = 0.5 * alpha[0] - 0.25 * alpha[2];
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
    for j in 1..(num_coefficients +1) {
        println!("{:e}", h[j-1]);
        parks_mcclellan_output.impulse_response.push(h[j-1]);
    }
    println!();
    if neg && odd {
        parks_mcclellan_output.impulse_response.push(0.0);
    }

    for k in 0..bands.len() {
        println!("Band {}:", k);
        println!("    lower edge: {}", bands[k].lower_edge);
        parks_mcclellan_output.lower_band_edges.push(bands[k].lower_edge);
        println!("    upper edge: {}", bands[k].upper_edge);
        parks_mcclellan_output.upper_band_edges.push(bands[k].upper_edge);
        match j_type {
            JType::MultipleBand => {
                println!("    desired value: {}", bands[k].desired_value);
                parks_mcclellan_output.desired_values.push(bands[k].desired_value);
            },
            JType::Differentiator => {
                println!("    desired slope: {}", bands[k].desired_value);
                parks_mcclellan_output.desired_values.push(bands[k].desired_value);
            },
            JType::Hilbert => {
                println!("    desired value: {}", bands[k].desired_value);
                parks_mcclellan_output.desired_values.push(bands[k].desired_value);
            },
        }
        println!("    weighting: {}", bands[k].weight);
        parks_mcclellan_output.weightings.push(bands[k].weight);

        let deviation = (dev/(bands[k].weight as f64)) as f32;
        println!("    deviation: {}", deviation);
        parks_mcclellan_output.deviations.push(deviation);
        match j_type {
            JType::MultipleBand => {
                let deviation_db: f32 = 20.0 * (deviation + bands[k].desired_value).log10();
                println!("    deviation in dB: {}", deviation_db);
                parks_mcclellan_output.deviation_dbs.push(deviation_db);
            },
            _ => {}
        }
    }

    println!("Extremal frequencies (maxima of the error curve)");
    for j in 1..(nz+1) {
        let ix = iext[j-1];
        let temp = grid.get_grid(ix as usize - 1);
        println!("{}", temp);
        parks_mcclellan_output.extremal_frequencies.push(temp);
    }

    parks_mcclellan_output
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
    grid: &DenseGrid,
    iext: &mut [i64; 66],
    alpha: &mut [f32; 66],
    dev: &mut f64,
) {
    // Function-scoped data
    let mut x = [0.0f64; 66];
    let mut y = [0.0f64; 66];
    let mut ad = [0.0; 66];
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


    'state_100: loop {
        iext[nzz-1] = grid.n_grid() as i64 + 1;
        niter += 1;
        println!("niter: {}", niter);
        if niter > itrmax {
            break 'state_100;
        }
        for j in 1..(nz+1) {
            let jxt = iext[j-1];
            let mut dtemp: f64 = grid.get_grid((jxt-1) as usize) as f64;
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
            let dtemp = ad[j-1] * grid.get_des((L-1) as usize) as f64;
            dnum += dtemp;
            let dtemp = (k as f64) * ad[j-1] / grid.get_wt((L-1) as usize) as f64;
            dden += dtemp;
            k = -k;
        }
        *dev = dnum / dden;
        println!("DEVIATION: {}", *dev);

        let nu = if *dev > 0.0 { -1 } else { 1 };
        *dev = -(nu as f64) * *dev;
        k = nu;

        for j in 1..(nz+1) {
            L = iext[j-1];
            let dtemp = (k as f64) * *dev / grid.get_wt((L-1) as usize) as f64;
            y[j-1] = grid.get_des((L-1) as usize) as f64 + dtemp;
            k = -k;
        }

        if *dev <= devl as f64 {
            println!("***** FAILURE TO CONVERGE *****");
            println!("Number of iterations: {}", niter);
            println!("If the number of iterations is greater than 3,");
            println!("the design might be correct, but should be verified by FFT.");
            break 'state_100;
        }

        devl = *dev as f32;
        jchnge = 0;
        k1 = iext[0];
        knz = iext[nz-1];
        klow = 0;
        nut = -nu;
        non_loop_j = 1;

        // Search for the extremal frequencies of the best approximation
        'state_200: loop {
            if non_loop_j == nzz {
                ynz = comp;
            }
            if non_loop_j >= nzz {
                if non_loop_j > nzz {
                    if luck > 9 {
                        let kn = iext[nzz-1];
                        for j in 1..(nfcns+1) {
                            iext[j-1] = iext[j] // j+1-1
                        }
                        iext[nz-1] = kn;
                        continue 'state_100;
                    }
                    if comp.unwrap() > y1.unwrap() {
                        y1 = comp;
                    }
                    k1 = iext[nzz-1];
                    L = grid.n_grid() as i64 + 1;
                    klow = knz;
                    nut = -nut1;
                    comp = Some(y1.unwrap() * 1.00001);
                    'loop_11: loop {
                        L = L-1;
                        if L <= klow {
                            if luck == 6 {
                                if jchnge > 0 {
                                    continue 'state_100;
                                }
                                break 'state_100;
                            }
                            for j in 1..(nfcns+1) {
                                let nzzmj = nzz - j;
                                let nzmj = nz - j;
                                iext[nzzmj-1] = iext[nzmj-1];
                            }
                            iext[0] = k1;
                            continue 'state_100;
                        }
                        err = grid.gee(None, &x, &y, &ad, L, nz) as f32;
                        err = (err - grid.get_des((L-1) as usize)) * grid.get_wt((L-1) as usize);
                        let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                        if dtemp <= 0.0 {
                            continue 'loop_11;
                        }
                        non_loop_j = nzz;
                        luck += 10;
                        common_loop_function_01(grid, iext, &x, &y, &ad, &mut non_loop_j, &mut klow, &mut jchnge, &mut L, &mut err, &mut comp, nut, nz);
                        continue 'state_200;
                    }
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
                'loop_06: loop {
                    L = L+1;
                    if L >= kup {
                        luck = 6;
                        L = grid.n_grid() as i64 + 1;
                        klow = knz;
                        nut = -nut1;
                        comp = Some(y1.unwrap() * 1.00001);
                        'loop_07: loop {
                            L = L-1;
                            if L <= klow {
                                if luck == 6 {
                                    if jchnge > 0 {
                                        continue 'state_100;
                                    }
                                    break 'state_100;
                                }
                                for j in 1..(nfcns+1) {
                                    let nzzmj = nzz - j;
                                    let nzmj = nz - j;
                                    iext[nzzmj-1] = iext[nzmj-1];
                                }
                                iext[0] = k1;
                                continue 'state_100;
                            }
                            err = grid.gee(None, &x, &y, &ad, L, nz) as f32;
                            err = (err - grid.get_des((L-1) as usize)) * grid.get_wt((L-1) as usize);
                            let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                            if dtemp <= 0.0 {
                                continue 'loop_07;
                            }
                            non_loop_j = nzz;
                            luck += 10;
                            common_loop_function_01(grid, iext, &x, &y, &ad, &mut non_loop_j, &mut klow, &mut jchnge, &mut L, &mut err, &mut comp, nut, nz);
                            continue 'state_200;
                        }
                    }
                    err = grid.gee(None, &x, &y, &ad, L, nz) as f32;
                    err = (err - grid.get_des((L-1) as usize)) * grid.get_wt((L-1) as usize);
                    let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                    if dtemp <= 0.0 {
                        continue 'loop_06;
                    }
                    comp = Some((nut as f64) * (err as f64));
                    non_loop_j = nzz;
                    common_loop_function_02(grid, iext, &x, &y, &ad, &mut non_loop_j, &mut klow, &mut jchnge, &mut L, &mut err, &mut comp, kup, nut, nz);
                    continue 'state_200;
                }
            }
            kup = iext[non_loop_j];
            L = iext[non_loop_j-1] + 1;
            nut = -nut;
            if non_loop_j == 2 {
                y1 = comp;
            }
            comp = Some(*dev);
            if L >= kup {
                L = L - 1;
                'loop_03: loop {
                    L = L - 1;
                    if L <= klow {
                        L = iext[non_loop_j-1] + 1;
                        if jchnge > 0 {
                            iext[non_loop_j-1] = L - 1;
                            non_loop_j = non_loop_j + 1;
                            klow = L - 1;
                            jchnge = jchnge + 1;
                            continue 'state_200;
                        }
                        'loop_05: loop {
                            L += 1;
                            if L >= kup {
                                klow = iext[non_loop_j-1];
                                non_loop_j += 1;
                                continue 'state_200;
                            }
                            err = grid.gee(None, &x, &y, &ad, L, nz) as f32;
                            err = (err - grid.get_des((L-1) as usize)) * grid.get_wt((L-1) as usize);
                            let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                            if dtemp <= 0.0 {
                                continue 'loop_05;
                            }
                            comp = Some((nut as f64) * (err as f64));
                            common_loop_function_02(grid, iext, &x, &y, &ad, &mut non_loop_j, &mut klow, &mut jchnge, &mut L, &mut err, &mut comp, kup, nut, nz);
                            continue 'state_200;
                        }
                    }
                    err = grid.gee(None, &x, &y, &ad, L, nz) as f32;
                    err = (err - grid.get_des((L-1) as usize)) * grid.get_wt((L-1) as usize);
                    let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                    if dtemp > 0.0 {
                        common_loop_function_01(grid, iext, &x, &y, &ad, &mut non_loop_j, &mut klow, &mut jchnge, &mut L, &mut err, &mut comp, nut, nz);
                        continue 'state_200;
                    }
                    if jchnge <= 0 {
                        continue 'loop_03;
                    }
                    klow = iext[non_loop_j-1];
                    non_loop_j += 1;
                    continue 'state_200;
                }
            }
            err = grid.gee(None, &x, &y, &ad, L, nz) as f32;
            err = (err - grid.get_des((L-1) as usize)) * grid.get_wt((L-1) as usize);
            let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
            if dtemp <= 0.0 {
                L = L - 1;
                'loop_13: loop {
                    L = L - 1;
                    if L <= klow {
                        L = iext[non_loop_j-1] + 1;
                        if jchnge > 0 {
                            iext[non_loop_j-1] = L - 1;
                            non_loop_j = non_loop_j + 1;
                            klow = L - 1;
                            jchnge = jchnge + 1;
                            continue 'state_200;
                        }
                        'loop_15: loop {
                            L += 1;
                            if L >= kup {
                                klow = iext[non_loop_j-1];
                                non_loop_j += 1;
                                continue 'state_200;
                            }
                            err = grid.gee(None, &x, &y, &ad, L, nz) as f32;
                            err = (err - grid.get_des((L-1) as usize)) * grid.get_wt((L-1) as usize);
                            let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                            if dtemp <= 0.0 {
                                continue 'loop_15;
                            }
                            comp = Some((nut as f64) * (err as f64));
                            common_loop_function_02(grid, iext, &x, &y, &ad, &mut non_loop_j, &mut klow, &mut jchnge, &mut L, &mut err, &mut comp, kup, nut, nz);
                            continue 'state_200;
                        }
                    }
                    err = grid.gee(None, &x, &y, &ad, L, nz) as f32;
                    err = (err - grid.get_des((L-1) as usize)) * grid.get_wt((L-1) as usize);
                    let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                    if dtemp > 0.0 {
                        common_loop_function_01(grid, iext, &x, &y, &ad, &mut non_loop_j, &mut klow, &mut jchnge, &mut L, &mut err, &mut comp, nut, nz);
                        continue 'state_200;
                    }
                    if jchnge <= 0 {
                        continue 'loop_13;
                    }
                    klow = iext[non_loop_j-1];
                    non_loop_j += 1;
                    continue 'state_200;
                }
            }
            comp = Some((nut as f64) * (err as f64));
            common_loop_function_02(grid, iext, &x, &y, &ad, &mut non_loop_j, &mut klow, &mut jchnge, &mut L, &mut err, &mut comp, kup, nut, nz);
            continue 'state_200;
        }
    }

    // Calculation of the coefficients of the best approximation using the inverse DFT
    // By here, we're done iterating the Remez exchange algorithm, and we're mostly just
    // calculating "outputs".
    let mut a = [0.0f64; 66];

    let nm1 = nfcns - 1;
    let fsh: f32 = 1.0e-06;
    x[nzz-1] = -2.0;
    let cn = 2 * nfcns - 1;
    let delf = 1.0f32 / (cn as f32);
    L = 1;
    let mut kkk = 0;
    if grid.get_grid(0) < 0.01 && grid.get_grid(grid.n_grid()-1) > 0.49 {
        kkk = 1;
    }
    if nfcns <= 3 {
        kkk = 1;
    }

    if kkk != 1 {
        let dtemp = (PI2 * grid.get_grid(0) as f64).cos();
        let dnum = (PI2 * grid.get_grid(grid.n_grid()-1) as f64).cos();
        aa = (2.0 / (dtemp - dnum)) as f32;
        bb = (-(dtemp + dnum) / (dtemp - dnum)) as f32;
    }

    for j in 1..(nfcns+1) {
        let mut ft = (j-1) as f32;
        ft = ft * delf;
        let mut xt: f32 = (PI2 * ft as f64).cos() as f32;
        if kkk != 1 {
            xt = (xt - bb) / aa;
            let xt1 = (1.0 - xt.powi(2)).sqrt();
            ft = (xt1.atan2(xt) as f64 / PI2) as f32;
        }

        'loop_01: loop {
            let xe: f32 = x[(L-1) as usize] as f32;
            if xt > xe {
                if (xt-xe) < fsh {
                    a[j-1] = y[(L-1) as usize];
                    break 'loop_01;
                } else {
                    a[j - 1] = grid.gee(Some(ft), &x, &y, &ad, 1, nz);
                    break 'loop_01;
                }
            } else if (xe-xt) < fsh {
                a[j-1] = y[(L-1) as usize];
                break 'loop_01;
            } else {
                L = L+1;
                continue 'loop_01;
            }
        }

        if L > 1 {
            L -= 1;
        }
    }

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
        alpha[j-1] = p[j-1] as f32;
    }

    if nfcns > 3 {
        return;
    }
    alpha[nfcns] = 0.0; // alpha[nfcns+1-1]
    alpha[nfcns+1] = 0.0; // alpha[nfcns+2-1]
    return;
}

// Function to calculate the lagrange interpolation coefficients
// for use in the function `grid.gee`
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

fn common_loop_function_01(
    grid: &DenseGrid,
    iext: &mut [i64; 66],
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    non_loop_j: &mut usize,
    klow: &mut i64,
    jchnge: &mut i64,
    L: &mut i64,
    err: &mut f32,
    comp: &mut Option<f64>,
    nut: i32,
    nz: usize,
) {
    *comp = Some((nut as f64) * (*err as f64));
    loop {
        *L = *L - 1;
        if L <= klow {
            *klow = iext[*non_loop_j-1];
            iext[*non_loop_j-1] = *L+1;
            *non_loop_j += 1;
            *jchnge += 1;
            break;
        }
        *err = grid.gee(None, x, y, ad, *L, nz) as f32;
        *err = (*err - grid.get_des((*L-1) as usize)) * grid.get_wt((*L-1) as usize);
        let dtemp = (nut as f64) * (*err as f64) - comp.unwrap();
        if dtemp <= 0.0 {
            *klow = iext[*non_loop_j-1];
            iext[*non_loop_j-1] = *L+1;
            *non_loop_j += 1;
            *jchnge += 1;
            break;
        }
        *comp = Some((nut as f64) * (*err as f64));
    }
}

fn common_loop_function_02(
    grid: &DenseGrid,
    iext: &mut [i64; 66],
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    non_loop_j: &mut usize,
    klow: &mut i64,
    jchnge: &mut i64,
    L: &mut i64,
    err: &mut f32,
    comp: &mut Option<f64>,
    kup: i64,
    nut: i32,
    nz: usize,
) {
    loop {
        *L = *L + 1;
        if *L >= kup {
            iext[*non_loop_j-1] = *L - 1;
            *non_loop_j = *non_loop_j + 1;
            *klow = *L - 1;
            *jchnge = *jchnge + 1;
            break;
        }
        *err = grid.gee(None, x, y, ad, *L, nz) as f32;
        *err = (*err - grid.get_des((*L-1) as usize)) * grid.get_wt((*L-1) as usize);
        let dtemp = (nut as f64) * (*err as f64) - comp.unwrap();
        if dtemp <= 0.0 {
            iext[*non_loop_j-1] = *L - 1;
            *non_loop_j = *non_loop_j + 1;
            *klow = *L - 1;
            *jchnge = *jchnge + 1;
            break;
        }
        *comp = Some((nut as f64) * (*err as f64));
    }
}

fn main() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0.0,
        upper_edge: 0.1,
        desired_value: 0.0,
        weight: 10.0,
    });
    bands.push(Band {
        lower_edge: 0.2,
        upper_edge: 0.35,
        desired_value: 1.0,
        weight: 1.0,
    });
    bands.push(Band {
        lower_edge: 0.425,
        upper_edge: 0.5,
        desired_value: 0.0,
        weight: 10.0,
    });

    let pm_output = design(32, JType::MultipleBand, &bands, 16);
}