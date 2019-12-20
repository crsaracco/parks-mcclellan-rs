#[cfg(test)] mod tests;

pub mod dense_grid;
use dense_grid::DenseGrid;

pub mod extremal_frequencies;
use extremal_frequencies::ExtremalFrequencies;

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
    // (not sure why we don't have to for "positive"-symmetry)
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
    let mut extremal_frequencies = ExtremalFrequencies::new(num_coefficients);
    extremal_frequencies.initialize_guess(&grid);

    let nm1 = num_coefficients - 1;
    let last_coefficient_index = num_coefficients + 1;

    // Call the remez exchange algorithm to do the approximation problem.
    let (alpha, deviation) = remez(num_coefficients, &grid, &mut extremal_frequencies);

    // Calculate the impulse response.
    let mut h = [0.0f32; 66];
    if !neg {
        if !odd {
            h[0] = 0.25 * alpha[num_coefficients -1];
            for j in 2..(nm1+1) {
                let nzmj = last_coefficient_index -j;
                let nf2j = num_coefficients +2-j;
                h[j-1] = 0.25 * (alpha[nzmj-1] + alpha[nf2j-1]);
            }
            h[num_coefficients -1] = 0.5 * alpha[0] + 0.25 * alpha[1];
        } else {
            for j in 1..(nm1+1) {
                let nzmj = last_coefficient_index -j;
                h[j-1] = 0.5 * alpha[nzmj-1];
            }
            h[num_coefficients -1] = alpha[0];
        }
    } else {
        if !odd {
            h[0] = 0.25 * alpha[num_coefficients -1];
            for j in 2..(nm1+1) {
                let nzmj = last_coefficient_index -j;
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
                let nzmj = last_coefficient_index -j;
                let nf3j = num_coefficients +3-j;
                h[j-1] = 0.25 * (alpha[nzmj-1] - alpha[nf3j-1]);
            }
            h[num_coefficients -1] = 0.5 * alpha[0] - 0.25 * alpha[2];
            h[last_coefficient_index -1] = 0.0;
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

        let deviation = (deviation /(bands[k].weight as f64)) as f32;
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
    for j in 1..(last_coefficient_index +1) {
        let grid_index = extremal_frequencies.get_grid_index(j-1);
        let temp = grid.get_grid(grid_index as usize - 1);
        println!("{}", temp);
        parks_mcclellan_output.extremal_frequencies.push(temp);
    }

    parks_mcclellan_output
}

// This function implements the Remez Exchange Algorithm for the weighted Chebyshev approximation
// of a continuous function with a sum of cosines.
// Inputs to the subroutine are:
//  - the number of cosines (= number of coefficients of the filter)
//  - a dense grid which replaces the frequency axis, along with the desired function on that grid
//  - an initial guess of the extremal frequencies
// The program minimizes the Chebyshev error by determining the best location of the extremal
// frequencies (points of maximum error) and then calculates the coefficients of the best
// approximation.
fn remez(
    num_coefficients: usize,
    grid: &DenseGrid,
    extremal_frequencies: &mut ExtremalFrequencies,
) -> ([f32; 66], f64) {
    // TODO: maybe move these into "state_100" and return them
    let mut x = [0.0f64; 66];
    let mut y = [0.0f64; 66];
    let mut ad = [0.0f64; 66];

    // Iterate to find the best approximation
    let deviation = state_100(&grid, num_coefficients, &mut x, &mut y, &mut ad, extremal_frequencies);

    // Calculate the coefficients of the best approximation using the inverse DFT.
    let alpha = calculate_alpha(num_coefficients, &grid, &mut x, &mut y, &mut ad);

    (alpha, deviation)
}

fn state_100(
    grid: &DenseGrid,
    num_coefficients: usize,
    x: &mut [f64; 66],
    y: &mut [f64; 66],
    ad: &mut [f64; 66],
    extremal_frequencies: &mut ExtremalFrequencies,
) -> f64 {
    let mut deviation = 0.0;

    let mut niter = 0;
    let itrmax = 25; // Max number of iterations

    let last_coefficient_index = num_coefficients + 1;
    let coefficient_off_end_index = num_coefficients + 2;

    let mut devl = -1.0f32;
    let mut luck = 0;
    let mut nut1 = 0;
    let mut ynz = None;
    let mut comp = None;
    let mut y1 = None;

    'state_100: loop {
        extremal_frequencies.set_grid_index(last_coefficient_index, grid.n_grid() as i64 + 1);
        niter += 1;
        println!("niter: {}", niter);
        if niter > itrmax {
            break 'state_100;
        }
        for j in 1..(last_coefficient_index +1) {
            let jxt = extremal_frequencies.get_grid_index(j-1);
            let mut dtemp: f64 = grid.get_grid((jxt-1) as usize) as f64;
            dtemp = (dtemp * PI2).cos();
            x[j-1] = dtemp;
        }
        let jet = (num_coefficients -1) / 15 + 1;
        for j in 1..(last_coefficient_index +1) {
            ad[j-1] = d_func(&x, j, last_coefficient_index, jet);
        }

        let mut dnum = 0.0;
        let mut dden = 0.0;
        let mut k = 1;
        for j in 1..(last_coefficient_index+1) {
            let local_ell = extremal_frequencies.get_grid_index(j-1);
            let dtemp = ad[j-1] * grid.get_des((local_ell-1) as usize) as f64;
            dnum += dtemp;
            let dtemp = (k as f64) * ad[j-1] / grid.get_wt((local_ell-1) as usize) as f64;
            dden += dtemp;
            k = -k;
        }
        deviation = dnum / dden;
        println!("DEVIATION: {}", deviation);

        let nu = if deviation > 0.0 { -1 } else { 1 };
        deviation = -(nu as f64) * deviation;
        k = nu;

        for j in 1..(last_coefficient_index +1) {
            let local_ell = extremal_frequencies.get_grid_index(j-1);
            let dtemp = (k as f64) * deviation / grid.get_wt((local_ell-1) as usize) as f64;
            y[j-1] = grid.get_des((local_ell-1) as usize) as f64 + dtemp;
            k = -k;
        }

        if deviation <= devl as f64 {
            println!("***** FAILURE TO CONVERGE *****");
            println!("Number of iterations: {}", niter);
            println!("If the number of iterations is greater than 3,");
            println!("the design might be correct, but should be verified by FFT.");
            break 'state_100;
        }

        devl = deviation as f32;
        let mut jchnge = 0;
        let mut k1 = extremal_frequencies.get_grid_index(0);
        let mut knz = extremal_frequencies.get_grid_index(last_coefficient_index-1);
        let mut klow = 0;
        let mut nut = -nu;
        let mut non_loop_j = 1; // Only ever increases inside 'state_200

        // Search for the extremal frequencies of the best approximation
        'state_200: loop {
            if non_loop_j == coefficient_off_end_index {
                ynz = comp;
            }
            if non_loop_j >= coefficient_off_end_index {
                if non_loop_j > coefficient_off_end_index {
                    if luck > 9 {
                        extremal_frequencies.shift_grid_indexes_left();
                        continue 'state_100;
                    }
                    if comp.unwrap() > y1.unwrap() {
                        y1 = comp;
                    }
                    k1 = extremal_frequencies.get_grid_index(last_coefficient_index);
                    let mut local_ell = grid.n_grid() as i64 + 1;
                    klow = knz;
                    nut = -nut1;
                    comp = Some(y1.unwrap() * 1.00001);
                    'loop_11: loop {
                        local_ell -= 1;
                        if local_ell <= klow {
                            if luck == 6 {
                                if jchnge > 0 {
                                    continue 'state_100;
                                }
                                break 'state_100;
                            }
                            extremal_frequencies.shift_grid_indexes_right(k1);
                            continue 'state_100;
                        }
                        let err = calculate_err(grid, None, &x, &y, &ad, local_ell, last_coefficient_index);
                        let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                        if dtemp <= 0.0 {
                            continue 'loop_11;
                        }
                        break 'loop_11;
                    }

                    loop {
                        let err = calculate_err(grid, None, &x, &y, &ad, local_ell, last_coefficient_index);
                        comp = Some((nut as f64) * (err as f64));
                        local_ell -= 1;
                        if local_ell <= klow {
                            break;
                        }
                        if (nut as f64) * (err as f64) - comp.unwrap() <= 0.0 {
                            break;
                        }
                    }
                    non_loop_j = coefficient_off_end_index;
                    klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                    extremal_frequencies.set_grid_index(non_loop_j-1, local_ell+1);
                    non_loop_j += 1;
                    jchnge += 1;
                    luck += 10;
                    continue 'state_200;
                }
                let zeroth_grid_index = extremal_frequencies.get_grid_index(0);
                if k1 > zeroth_grid_index {
                    k1 = zeroth_grid_index;
                }
                let almost_last_grid_index = extremal_frequencies.get_grid_index(last_coefficient_index-1);
                if knz < almost_last_grid_index {
                    knz = almost_last_grid_index;
                }
                nut1 = nut;
                nut = -nu;
                let mut ell = 0;
                let kup = k1;
                comp = Some(ynz.unwrap() * 1.00001);
                luck = 1;
                'loop_06: loop {
                    ell += 1;
                    if ell >= kup {
                        luck = 6;
                        ell = grid.n_grid() as i64 + 1;
                        klow = knz;
                        nut = -nut1;
                        comp = Some(y1.unwrap() * 1.00001);

                        ell = ell-1;
                        if ell <= klow {
                            if jchnge > 0 {
                                continue 'state_100;
                            }
                            break 'state_100;
                        }
                        let mut err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);

                        'loop_07: while (nut as f64) * (err as f64) - comp.unwrap() <= 0.0 {
                            ell = ell-1;
                            if ell <= klow {
                                if jchnge > 0 {
                                    continue 'state_100;
                                }
                                break 'state_100;
                            }
                            err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                        }

                        non_loop_j = coefficient_off_end_index;
                        luck += 10;
                        comp = Some((nut as f64) * (err as f64));

                        loop {
                            ell -= 1;
                            if ell <= klow {
                                break;
                            }
                            err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                            let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                            if dtemp <= 0.0 {
                                break;
                            }
                            comp = Some((nut as f64) * (err as f64));
                        }

                        klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                        extremal_frequencies.set_grid_index(non_loop_j-1, ell+1);
                        non_loop_j += 1;
                        jchnge += 1;
                        continue 'state_200;
                    }
                    let mut err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                    let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                    if dtemp <= 0.0 {
                        continue 'loop_06;
                    }
                    non_loop_j = coefficient_off_end_index;

                    comp = Some((nut as f64) * (err as f64));

                    loop {
                        ell = ell + 1;
                        if ell >= kup {
                            break;
                        }
                        err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                        let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                        if dtemp <= 0.0 {
                            break;
                        }
                        comp = Some((nut as f64) * (err as f64));
                    }

                    extremal_frequencies.set_grid_index(non_loop_j-1, ell-1);
                    non_loop_j = non_loop_j + 1;
                    klow = ell - 1;
                    jchnge = jchnge + 1;

                    continue 'state_200;
                }
            }
            let kup = extremal_frequencies.get_grid_index(non_loop_j);
            let mut ell = extremal_frequencies.get_grid_index(non_loop_j-1) + 1;
            nut = -nut;
            if non_loop_j == 2 {
                y1 = comp;
            }
            comp = Some(deviation);
            if ell >= kup {
                ell = ell - 1;
                'loop_03: loop {
                    ell = ell - 1;
                    if ell <= klow {
                        ell = extremal_frequencies.get_grid_index(non_loop_j-1) + 1;
                        if jchnge > 0 {
                            extremal_frequencies.set_grid_index(non_loop_j-1, ell-1);
                            non_loop_j = non_loop_j + 1;
                            klow = ell - 1;
                            jchnge = jchnge + 1;
                            continue 'state_200;
                        }
                        'loop_05: loop {
                            ell += 1;
                            if ell >= kup {
                                klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                                non_loop_j += 1;
                                continue 'state_200;
                            }
                            let mut err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                            let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                            if dtemp <= 0.0 {
                                continue 'loop_05;
                            }
                            comp = Some((nut as f64) * (err as f64));

                            loop {
                                ell = ell + 1;
                                if ell >= kup {
                                    break;
                                }
                                err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                                let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                                if dtemp <= 0.0 {
                                    break;
                                }
                                comp = Some((nut as f64) * (err as f64));
                            }

                            extremal_frequencies.set_grid_index(non_loop_j-1, ell-1);
                            non_loop_j = non_loop_j + 1;
                            klow = ell - 1;
                            jchnge = jchnge + 1;
                            continue 'state_200;
                        }
                    }
                    let mut err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                    let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                    if dtemp > 0.0 {
                        comp = Some((nut as f64) * (err as f64));

                        loop {
                            ell -= 1;
                            if ell <= klow {
                                break;
                            }
                            err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                            let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                            if dtemp <= 0.0 {
                                break;
                            }
                            comp = Some((nut as f64) * (err as f64));
                        }

                        klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                        extremal_frequencies.set_grid_index(non_loop_j-1, ell+1);
                        non_loop_j += 1;
                        jchnge += 1;
                        continue 'state_200;
                    }
                    if jchnge <= 0 {
                        continue 'loop_03;
                    }
                    klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                    non_loop_j += 1;
                    continue 'state_200;
                }
            }
            let mut err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
            let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
            if dtemp <= 0.0 {
                ell = ell - 1;
                'loop_13: loop {
                    ell = ell - 1;
                    if ell <= klow {
                        ell = extremal_frequencies.get_grid_index(non_loop_j-1) + 1;
                        if jchnge > 0 {
                            extremal_frequencies.set_grid_index(non_loop_j-1, ell-1);
                            non_loop_j = non_loop_j + 1;
                            klow = ell - 1;
                            jchnge = jchnge + 1;
                            continue 'state_200;
                        }
                        'loop_15: loop {
                            ell += 1;
                            if ell >= kup {
                                klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                                non_loop_j += 1;
                                continue 'state_200;
                            }
                            err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                            let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                            if dtemp <= 0.0 {
                                continue 'loop_15;
                            }
                            comp = Some((nut as f64) * (err as f64));

                            loop {
                                ell = ell + 1;
                                if ell >= kup {
                                    break;
                                }
                                err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                                let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                                if dtemp <= 0.0 {
                                    break;
                                }
                                comp = Some((nut as f64) * (err as f64));
                            }

                            extremal_frequencies.set_grid_index(non_loop_j-1, ell-1);
                            non_loop_j = non_loop_j + 1;
                            klow = ell - 1;
                            jchnge = jchnge + 1;
                            continue 'state_200;
                        }
                    }
                    err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                    let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                    if dtemp > 0.0 {
                        comp = Some((nut as f64) * (err as f64));

                        loop {
                            ell -= 1;
                            if ell <= klow {
                                break;
                            }
                            err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                            let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                            if dtemp <= 0.0 {
                                break;
                            }
                            comp = Some((nut as f64) * (err as f64));
                        }

                        klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                        extremal_frequencies.set_grid_index(non_loop_j-1, ell+1);
                        non_loop_j += 1;
                        jchnge += 1;
                        continue 'state_200;
                    }
                    if jchnge <= 0 {
                        continue 'loop_13;
                    }
                    klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                    non_loop_j += 1;
                    continue 'state_200;
                }
            }
            comp = Some((nut as f64) * (err as f64));

            loop {
                ell = ell + 1;
                if ell >= kup {
                    break;
                }
                err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                let dtemp = (nut as f64) * (err as f64) - comp.unwrap();
                if dtemp <= 0.0 {
                    break;
                }
                comp = Some((nut as f64) * (err as f64));
            }

            extremal_frequencies.set_grid_index(non_loop_j-1, ell-1);
            non_loop_j = non_loop_j + 1;
            klow = ell - 1;
            jchnge = jchnge + 1;
            continue 'state_200;
        }
    }

    deviation
}

fn calculate_alpha(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &mut [f64; 66],
    y: &mut [f64; 66],
    ad: &mut [f64; 66],
) -> [f32; 66] {
    let last_coefficient_index = num_coefficients + 1;

    let mut a = [0.0f64; 66];
    let mut aa = 0.0f32;
    let mut bb = 0.0f32;
    let mut kkk = 0;

    let fsh: f32 = 1.0e-06;
    x[last_coefficient_index] = -2.0;
    let cn = 2 * num_coefficients - 1;
    let delf = 1.0f32 / (cn as f32);
    let mut ell = 1;

    if grid.get_grid(0) < 0.01 && grid.get_grid(grid.n_grid()-1) > 0.49 {
        kkk = 1;
    }
    if num_coefficients <= 3 {
        kkk = 1;
    }

    if kkk != 1 {
        let dtemp = (PI2 * grid.get_grid(0) as f64).cos();
        let dnum = (PI2 * grid.get_grid(grid.n_grid()-1) as f64).cos();
        aa = (2.0 / (dtemp - dnum)) as f32;
        bb = (-(dtemp + dnum) / (dtemp - dnum)) as f32;
    }

    for j in 1..(num_coefficients +1) {
        let mut ft = (j-1) as f32;
        ft = ft * delf;
        let mut xt: f32 = (PI2 * ft as f64).cos() as f32;
        if kkk != 1 {
            xt = (xt - bb) / aa;
            let xt1 = (1.0 - xt.powi(2)).sqrt();
            ft = (xt1.atan2(xt) as f64 / PI2) as f32;
        }

        'loop_01: loop {
            let xe: f32 = x[(ell -1) as usize] as f32;
            if xt > xe {
                if (xt-xe) < fsh {
                    a[j-1] = y[(ell -1) as usize];
                    break 'loop_01;
                } else {
                    a[j - 1] = grid.gee(Some(ft), &x, &y, &ad, 1, last_coefficient_index);
                    break 'loop_01;
                }
            } else if (xe-xt) < fsh {
                a[j-1] = y[(ell -1) as usize];
                break 'loop_01;
            } else {
                ell = ell +1;
                continue 'loop_01;
            }
        }

        if ell > 1 {
            ell -= 1;
        }
    }

    let mut alpha = [0.0f32; 66];

    let cn = 2 * num_coefficients - 1;
    let dden = PI2 / (cn as f64);
    let nm1 = num_coefficients - 1;

    for j in 1..(num_coefficients + 1) {
        let mut dtemp = 0.0;
        let mut dnum = (j-1) as f64;
        dnum = dnum * dden;
        if nm1 >= 1 {
            for k in 1..(nm1+1) {
                let dak = a[k];
                let dk = k as f64;
                dtemp = dtemp + dak * (dnum*dk).cos();
            }
        }
        dtemp = 2.0 * dtemp + a[0];
        alpha[j-1] = dtemp as f32;
    }

    for j in 2..(num_coefficients +1) {
        alpha[j-1] = 2.0 * alpha[j-1] / (cn as f32);
    }
    alpha[0] = alpha[0] / (cn as f32);

    if kkk == 1 {
        if num_coefficients > 3 {
            return alpha;
        }
        alpha[num_coefficients] = 0.0; // alpha[nfcns+1-1]
        alpha[num_coefficients +1] = 0.0; // alpha[nfcns+2-1]
        return alpha;
    }

    let mut p = [0.0f64; 65];
    let mut q = [0.0f64; 65];

    p[0] = (2.0 * alpha[num_coefficients -1] * bb + alpha[nm1-1]) as f64;
    p[1] = (2.0 * aa * alpha[num_coefficients -1]) as f64;
    q[0] = (alpha[num_coefficients -3] - alpha[num_coefficients -1]) as f64;

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
            let nf1j = num_coefficients - 1 - j;
            q[0] = q[0] + (alpha[nf1j - 1] as f64);
        }
    }

    for j in 1..(num_coefficients +1) {
        alpha[j-1] = p[j-1] as f32;
    }

    if num_coefficients > 3 {
        return alpha;
    }
    alpha[num_coefficients] = 0.0;
    alpha[num_coefficients +1] = 0.0;

    alpha
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

fn calculate_err(
    grid: &DenseGrid,
    zeroth_value_override: Option<f32>,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    ell: i64,
    last_coefficient_index: usize,
) -> f32 {
    let err = grid.gee(zeroth_value_override, x, y, ad, ell, last_coefficient_index) as f32;
    (err - grid.get_des((ell - 1) as usize)) * grid.get_wt((ell - 1) as usize)
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

    let _pm_output = design(32, JType::MultipleBand, &bands, 16);
}