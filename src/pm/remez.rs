use super::endpoints_search;
use super::extremal_frequencies::ExtremalFrequencies;
use super::extremal_frequency_search;
use super::dense_grid::DenseGrid;
use super::lagrange_interpolation;
use super::PI2;

// This function implements the Remez Exchange Algorithm for the weighted Chebyshev approximation
// of a continuous function with a sum of cosines.
// Inputs to the subroutine are:
//  - the number of cosines (= number of coefficients of the filter)
//  - a dense grid which replaces the frequency axis, along with the desired function on that grid
//  - an initial guess of the extremal frequencies
// The program minimizes the Chebyshev error by determining the best location of the extremal
// frequencies (points of maximum error) and then calculates the coefficients of the best
// approximation.
pub fn remez(
    num_coefficients: usize,
    grid: &DenseGrid,
) -> ([f32; 66], f64, ExtremalFrequencies) {
    // TODO: maybe move these into "remez_iterate" and return them?
    let mut x = [0.0f64; 66];
    let mut y = [0.0f64; 66];
    let mut ad = [0.0f64; 66];

    // Initial guess for the extremal frequencies: equally spaced along the grid.
    let mut extremal_frequencies = ExtremalFrequencies::new(num_coefficients);
    extremal_frequencies.initialize_guess(&grid);

    // Iterate to find the best approximation
    let deviation = remez_iterate(num_coefficients, &grid, &mut extremal_frequencies, &mut x, &mut y, &mut ad);

    // Calculate the coefficients of the best approximation using the inverse DFT.
    let alpha = calculate_alpha(num_coefficients, &grid, &mut x, &mut y, &mut ad);

    (alpha, deviation, extremal_frequencies)
}

fn remez_iterate(
    num_coefficients: usize,
    grid: &DenseGrid,
    extremal_frequencies: &mut ExtremalFrequencies,
    x: &mut [f64; 66],
    y: &mut [f64; 66],
    ad: &mut [f64; 66],
) -> f64 {
    let last_coefficient_index = num_coefficients + 1;

    let mut deviation = 0.0;
    let mut devl = -1.0f32;

    // TODO: make num_iterations max a parameter
    for _ in 0..25 {
        extremal_frequencies.set_grid_index(last_coefficient_index, grid.n_grid() as i64 + 1);

        // Lagrange interpolation
        *x = lagrange_interpolation::lagrange_x_coordinates(num_coefficients, grid, extremal_frequencies);
        *ad = lagrange_interpolation::lagrange_interpolation_coefficients(num_coefficients, x);
        deviation = lagrange_interpolation::deviation(num_coefficients, grid, ad, extremal_frequencies);
        // println!("DEVIATION: {}", deviation);
        let nu = if deviation > 0.0 { -1 } else { 1 };
        deviation = deviation.abs();
        *y = lagrange_interpolation::lagrange_y_coordinates(num_coefficients, grid, nu, deviation, extremal_frequencies);

        // If deviation this time is less than (or equal to) the deviation last time (devl),
        // then we've failed to converge. Break out of iteration, and calculate filter coefficients.
        if deviation <= devl as f64 {
            // TODO: represent this as a Result<T> type, or maybe your own enum type.
            // println!("***** FAILURE TO CONVERGE *****");
            // println!("Number of iterations: {}", num_iterations);
            // println!("If the number of iterations is greater than 3,");
            // println!("the design might be correct, but should be verified by FFT.");
            break;
        }

        // Record deviation this time, to use for the above comparison next time.
        devl = deviation as f32;

        // Search for the extremal frequencies of the best approximation
        let recalculation_result = recalculate_extremal_frequencies(
            num_coefficients,
            grid,
            nu,
            x,
            y,
            ad,
            deviation,
            extremal_frequencies
        );

        if recalculation_result.is_err() {
            break;
        }
    }

    deviation
}

// If Ok(()): keep iterating.
// If Err(()): stop iterating.
fn recalculate_extremal_frequencies(
    num_coefficients: usize,
    grid: &DenseGrid,
    nu: i32,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    deviation: f64,
    extremal_frequencies: &mut ExtremalFrequencies,
) -> Result<(), ()> {
    use endpoints_search::EndpointSearchResult;

    let left_nu = nu;
    let mut right_nu = -nu;
    let mut comp = 0.0;

    // Capture k1 and knz before the extremal frequencies are changed by `find_extremal_frequencies`
    let old_k1 = extremal_frequencies.get_grid_index(0);
    let old_knz = extremal_frequencies.get_grid_index(num_coefficients);

    // Find the extremal frequencies
    let extremal_frequencies_changed = extremal_frequency_search::find_extremal_frequencies(
        num_coefficients,
        grid,
        x,
        y,
        ad,
        deviation,
        &mut right_nu,
        &mut comp,
        extremal_frequencies,
    );

    // Search for the endpoints
    let endpoint_search_result = endpoints_search::endpoints_search(
        num_coefficients,
        grid,
        x,
        y,
        ad,
        left_nu,
        right_nu,
        old_k1,
        old_knz,
        deviation,
        extremal_frequencies_changed,
        extremal_frequencies,
    );

    match endpoint_search_result {
        EndpointSearchResult::KeepIteratingRemez => return Ok(()),
        EndpointSearchResult::StopIteratingRemez => return Err(()),
    };
}

// TODO: Desperately needs to be refactored.
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
    let mut kk = false;

    let fsh: f32 = 1.0e-06;
    x[last_coefficient_index] = -2.0;
    let cn = 2 * num_coefficients - 1;
    let delf = 1.0f32 / (cn as f32);
    let mut ell = 1;

    if grid.get_grid(0) < 0.01 && grid.get_grid(grid.n_grid()-1) > 0.49 {
        kk = true;
    }
    if num_coefficients <= 3 {
        kk = true;
    }

    if !kk {
        let dtemp = (PI2 * grid.get_grid(0) as f64).cos();
        let dnum = (PI2 * grid.get_grid(grid.n_grid()-1) as f64).cos();
        aa = (2.0 / (dtemp - dnum)) as f32;
        bb = (-(dtemp + dnum) / (dtemp - dnum)) as f32;
    }

    for j in 1..(num_coefficients +1) {
        let mut ft = (j-1) as f32;
        ft = ft * delf;
        let mut xt: f32 = (PI2 * ft as f64).cos() as f32;
        if !kk {
            xt = (xt - bb) / aa;
            let xt1 = (1.0 - xt.powi(2)).sqrt();
            ft = (xt1.atan2(xt) as f64 / PI2) as f32;
        }

        loop {
            let xe: f32 = x[(ell -1) as usize] as f32;
            if xt > xe {
                if (xt-xe) < fsh {
                    a[j-1] = y[(ell -1) as usize];
                    break;
                } else {
                    a[j - 1] = grid.frequency_response(Some(ft), &x, &y, &ad, 1, last_coefficient_index);
                    break;
                }
            } else if (xe-xt) < fsh {
                a[j-1] = y[(ell -1) as usize];
                break;
            } else {
                ell = ell +1;
                continue;
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

    if kk {
        if num_coefficients > 3 {
            return alpha;
        }
        alpha[num_coefficients] = 0.0;
        alpha[num_coefficients +1] = 0.0;
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
        p[j] = 0.0;
        for k in 1..(j + 1) {
            a[k-1] = p[k-1];
            p[k-1] = 2.0 * (bb as f64) * a[k - 1];
        }
        p[1] = p[1] + a[0] * 2.0 * (aa as f64);
        let jm1 = j - 1;
        for k in 1..(jm1 + 1) {
            p[k - 1] = p[k - 1] + q[k - 1] + (aa as f64) * a[k];
        }
        let jp1 = j + 1;
        for k in 3..(jp1 + 1) {
            p[k - 1] = p[k - 1] + (aa as f64) * a[k - 2];
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

pub fn calculate_impulse_response(alpha: &[f32; 66], num_coefficients: usize, odd: bool, neg: bool) -> Vec<f32> {
    let mut impulse_response = vec![0.0; 66];
    match (odd, neg) {
        (false, false) => {
            // Impulse response has an EVEN number of samples,
            // REFLECTED horizontally across the midpoint of the middle two points.
            impulse_response[0] = 0.25 * alpha[num_coefficients - 1];
            for j in 2..num_coefficients {
                impulse_response[j - 1] = 0.25 * (alpha[num_coefficients - j] + alpha[num_coefficients + 1 - j]);
            }
            impulse_response[num_coefficients - 1] = 0.5 * alpha[0] + 0.25 * alpha[1];
        },
        (false, true) => {
            // Impulse response has an EVEN number of samples,
            // ROTATED 180 degrees around the midpoint of the middle two points.
            impulse_response[0] = 0.25 * alpha[num_coefficients - 1];
            for j in 2..num_coefficients {
                impulse_response[j - 1] = 0.25 * (alpha[num_coefficients - j] - alpha[num_coefficients + 1 - j]);
            }
            impulse_response[num_coefficients - 1] = 0.5 * alpha[0] - 0.25 * alpha[1];
        },
        (true, false) => {
            // Impulse response has an ODD number of samples,
            // REFLECTED horizontally across the middle point.
            for j in 1..num_coefficients {
                impulse_response[j - 1] = 0.5 * alpha[num_coefficients - j];
            }
            impulse_response[num_coefficients - 1] = alpha[0];
        },
        (true, true) => {
            // Impulse response has an ODD number of samples,
            // ROTATED 180 degrees around the middle point.
            impulse_response[0] = 0.25 * alpha[num_coefficients - 1];
            if num_coefficients > 1 {
                impulse_response[1] = 0.25 * alpha[num_coefficients - 2];
            }
            for j in 3..num_coefficients {
                impulse_response[j - 1] = 0.25 * (alpha[num_coefficients - j] - alpha[num_coefficients + 2 - j]);
            }
            impulse_response[num_coefficients -1] = 0.5 * alpha[0] - 0.25 * alpha[2];
        }
    }
    impulse_response
}