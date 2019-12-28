#[cfg(test)] mod tests;

pub mod dense_grid;
use dense_grid::DenseGrid;

pub mod extremal_frequencies;
use extremal_frequencies::ExtremalFrequencies;

pub mod lagrange_interpolation;
pub mod extremal_frequency_search;
pub mod endpoints_search;

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
    pub filter_type: JType,
    pub impulse_response: Vec<f32>,
    pub lower_band_edges: Vec<f32>,
    pub upper_band_edges: Vec<f32>,
    pub desired_values: Vec<f32>,
    pub weightings: Vec<f32>,
    pub deviations: Vec<f32>,
    pub deviation_dbs: Vec<f32>,
    pub extremal_frequencies: Vec<f32>,
}

fn calculate_impulse_response(alpha: &[f32; 66], num_coefficients: usize, odd: bool, neg: bool) -> Vec<f32> {
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

fn print_filter_outputs(parks_mcclellan_output: &ParksMcclellanOutput) {
    match parks_mcclellan_output.filter_type {
        JType::MultipleBand => println!("BANDPASS FILTER"),
        JType::Differentiator => println!("DIFFERENTIATOR"),
        JType::Hilbert => println!("HILBERT TRANSFORMER"),
    }
    println!("Filter length: {}", parks_mcclellan_output.filter_length);
    println!();
    println!("***** IMPULSE RESPONSE *****");
    for n in parks_mcclellan_output.impulse_response.iter() {
        println!("{:e}", *n);
    }
    println!();
    for k in 0..parks_mcclellan_output.lower_band_edges.len() {
        println!("Band {}:", k);
        println!("    lower edge: {}", parks_mcclellan_output.lower_band_edges[k]);
        println!("    upper edge: {}", parks_mcclellan_output.upper_band_edges[k]);
        match parks_mcclellan_output.filter_type {
            JType::MultipleBand => println!("    desired value: {}", parks_mcclellan_output.desired_values[k]),
            JType::Differentiator => println!("    desired slope: {}", parks_mcclellan_output.desired_values[k]),
            JType::Hilbert => println!("    desired value: {}", parks_mcclellan_output.desired_values[k]),
        }
        println!("    weighting: {}", parks_mcclellan_output.weightings[k]);
        println!("    deviation: {}", parks_mcclellan_output.deviations[k]);
        match parks_mcclellan_output.filter_type {
            JType::MultipleBand => println!("    deviation in dB: {}", parks_mcclellan_output.deviation_dbs[k]),
            _ => {}
        }
    }
    println!("Extremal frequencies (maxima of the error curve):");
    for e in parks_mcclellan_output.extremal_frequencies.iter() {
        println!("{}", *e);
    }
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

    // Call the remez exchange algorithm to do the approximation problem.
    let (alpha, deviation, extremal_frequencies) = remez(num_coefficients, &grid);

    // Calculate the impulse response.
    let impulse_response = calculate_impulse_response(&alpha, num_coefficients, odd, neg);

    // Program output section
    let mut parks_mcclellan_output = ParksMcclellanOutput {
        filter_length: n_filt,
        filter_type: j_type,
        impulse_response: vec![],
        lower_band_edges: vec![],
        upper_band_edges: vec![],
        desired_values: vec![],
        weightings: vec![],
        deviations: vec![],
        deviation_dbs: vec![],
        extremal_frequencies: vec![],
    };
    for j in 0..num_coefficients {
        parks_mcclellan_output.impulse_response.push(impulse_response[j]);
    }
    if neg && odd {
        parks_mcclellan_output.impulse_response.push(0.0);
    }
    for k in 0..bands.len() {
        parks_mcclellan_output.lower_band_edges.push(bands[k].lower_edge);
        parks_mcclellan_output.upper_band_edges.push(bands[k].upper_edge);
        match j_type {
            JType::MultipleBand => parks_mcclellan_output.desired_values.push(bands[k].desired_value),
            JType::Differentiator => parks_mcclellan_output.desired_values.push(bands[k].desired_value),
            JType::Hilbert => parks_mcclellan_output.desired_values.push(bands[k].desired_value),
        }
        parks_mcclellan_output.weightings.push(bands[k].weight);

        let deviation = (deviation /(bands[k].weight as f64)) as f32;
        parks_mcclellan_output.deviations.push(deviation);
        match j_type {
            JType::MultipleBand => {
                let deviation_db: f32 = 20.0 * (deviation + bands[k].desired_value).log10();
                parks_mcclellan_output.deviation_dbs.push(deviation_db);
            },
            _ => {}
        }
    }
    for j in 0..(num_coefficients + 1) {
        let grid_index = extremal_frequencies.get_grid_index(j);
        let temp = grid.get_grid(grid_index as usize - 1);
        parks_mcclellan_output.extremal_frequencies.push(temp);
    }

    print_filter_outputs(&parks_mcclellan_output);

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
    for num_iterations in 0..25 {
        println!("num_iterations: {}", num_iterations);

        extremal_frequencies.set_grid_index(last_coefficient_index, grid.n_grid() as i64 + 1);

        // Lagrange interpolation
        *x = lagrange_interpolation::lagrange_x_coordinates(num_coefficients, grid, extremal_frequencies);
        *ad = lagrange_interpolation::lagrange_interpolation_coefficients(num_coefficients, x);
        deviation = lagrange_interpolation::deviation(num_coefficients, grid, ad, extremal_frequencies);
        println!("DEVIATION: {}", deviation);
        let nu = if deviation > 0.0 { -1 } else { 1 };
        deviation = deviation.abs();
        *y = lagrange_interpolation::lagrange_y_coordinates(num_coefficients, grid, nu, deviation, extremal_frequencies);

        // If deviation this time is less than (or equal to) the deviation last time (devl),
        // then we've failed to converge. Break out of iteration, and calculate filter coefficients.
        if deviation <= devl as f64 {
            println!("***** FAILURE TO CONVERGE *****");
            println!("Number of iterations: {}", num_iterations);
            println!("If the number of iterations is greater than 3,");
            println!("the design might be correct, but should be verified by FFT.");
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

    let mut nut = -nu;
    let mut comp = 0.0;

    // Capture k1 and knz before the extremal frequencies are changed by `find_nth_extremal_frequency`
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
        &mut nut,
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
        nu,
        nut,
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