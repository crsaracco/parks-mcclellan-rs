use super::ExtremalFrequencies;
use super::DenseGrid;
use super::PI2;

pub fn lagrange_x_coordinates(
    num_coefficients: usize,
    grid: &DenseGrid,
    extremal_frequencies: &mut ExtremalFrequencies,
) -> [f64; 66] {
    let mut x = [0.0; 66];
    let last_coefficient_index = num_coefficients + 1;
    for j in 1..(last_coefficient_index+1) {
        let jxt = extremal_frequencies.get_grid_index(j-1);
        let mut dtemp: f64 = grid.get_grid((jxt-1) as usize) as f64;
        dtemp = (dtemp * PI2).cos();
        x[j-1] = dtemp;
    }
    x
}

pub fn lagrange_interpolation_coefficients(
    num_coefficients: usize,
    x: &[f64; 66],
) -> [f64; 66] {
    let mut ad = [0.0; 66];
    let last_coefficient_index = num_coefficients + 1;
    let jet = (num_coefficients - 1) / 15 + 1;
    for j in 1..(last_coefficient_index+1) {
        ad[j-1] = d_func(x, j, last_coefficient_index, jet);
    }
    ad
}

// Calculate the deviation between the *desired* filter function
// and the current Lagrange-interpolated function.
pub fn deviation(
    num_coefficients: usize,
    grid: &DenseGrid,
    ad: &[f64; 66],
    extremal_frequencies: &mut ExtremalFrequencies,
) -> f64 {
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    let mut k = 1;
    let last_coefficient_index = num_coefficients + 1;
    for j in 1..(last_coefficient_index+1) {
        let grid_index = extremal_frequencies.get_grid_index(j-1);
        numerator += ad[j-1] * grid.get_des((grid_index-1) as usize) as f64;
        denominator += (k as f64) * ad[j-1] / grid.get_wt((grid_index-1) as usize) as f64;
        k = -k;
    }
    numerator / denominator
}

pub fn lagrange_y_coordinates(
    num_coefficients: usize,
    grid: &DenseGrid,
    nu: i32,
    deviation: f64,
    extremal_frequencies: &mut ExtremalFrequencies,
) -> [f64; 66] {
    let mut y = [0.0; 66];
    let mut k = nu;
    let last_coefficient_index = num_coefficients + 1;
    for j in 1..(last_coefficient_index+1) {
        let grid_index = extremal_frequencies.get_grid_index(j-1);
        let temp = (k as f64) * deviation / grid.get_wt((grid_index-1) as usize) as f64;
        y[j-1] = grid.get_des((grid_index-1) as usize) as f64 + temp;
        k = -k;
    }
    y
}

// Function to calculate the lagrange interpolation coefficients
// for use in the function `lagrange_interpolation_coefficients`
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