use crate::Band;
use super::FilterType;

const PI: f64 = std::f64::consts::PI;
const PI2: f64 = PI * 2.0;

pub struct DenseGrid {
    grid: Vec<f32>,
    desired_values: [f32; 1045], // TODO: vec with exact length needed?
    weights: [f32; 1045], // TODO: vec with exact length needed?
    n_grid: usize,
}

impl DenseGrid {
    pub fn new(bands: &[Band], j_type: FilterType, grid_density: usize, num_coefficients: usize, rotational_symmetry: bool, odd_length: bool) -> Self {
        let delta_f = 0.5 / ((grid_density as f32) * (num_coefficients as f32));
        let (grid, grid_band) = generate_grid(bands, delta_f, rotational_symmetry);
        let desired_values = generate_desired_values(bands, j_type, &grid, &grid_band);
        let weights = generate_weights(bands, j_type, &grid, &grid_band);

        let mut n_grid = grid.len();
        if rotational_symmetry == odd_length && grid[n_grid-1] > (0.5 - delta_f) {
            n_grid -= 1;
        }

        Self {
            grid,
            desired_values,
            weights,
            n_grid,
        }
    }

    pub fn get_grid(&self, index: usize) -> f32 {
        self.grid[index]
    }

    pub fn get_des(&self, index: usize) -> f32 {
        self.desired_values[index]
    }
    
    pub fn get_wt(&self, index: usize) -> f32 {
        self.weights[index]
    }

    pub fn n_grid(&self) -> usize {
        self.n_grid
    }

    pub fn reformat_grid_for_remez(&mut self, neg: bool, odd: bool) {
        if !neg && odd {
            // Don't have to reformat the grid for this case.
            return;
        }

        let f = if neg { f64::sin } else { f64::cos };
        let constant = if odd { PI2 } else { PI };

        for j in 0..self.n_grid {
            let change = f(constant * self.grid[j] as f64) as f32;
            self.desired_values[j] = self.desired_values[j] / change;
            self.weights[j] = self.weights[j] * change;
        }
    }

    // Find the frequency response using lagrange interpolation
    pub fn frequency_response(&self, zeroth_value_override: Option<f32>, x: &[f64; 66], y: &[f64; 66], ad: &[f64; 66], k: i64, n: usize) -> f64 {
        let mut p = 0.0;
        let mut d = 0.0;

        let mut xf = if k-1 == 0 && zeroth_value_override.is_some() {
            zeroth_value_override.unwrap() as f64
        } else {
            self.grid[(k-1) as usize] as f64
        };

        xf = (PI2 * xf).cos();
        for j in 1..(n+1) {
            let mut c = xf - x[j-1];
            c = ad[j-1] / c;
            d += c;
            p += c*y[j-1];
        }

        p/d
    }

    // Calculate the error between the desired value and the currently-calculated frequency response
    pub fn calculate_error(
        &self,
        zeroth_value_override: Option<f32>,
        x: &[f64; 66],
        y: &[f64; 66],
        ad: &[f64; 66],
        k: i64,
        last_coefficient_index: usize,
    ) -> f32 {
        let frequency_response = self.frequency_response(zeroth_value_override, x, y, ad, k, last_coefficient_index) as f32;

        // Error is the delta between current frequency response and desired frequency response
        let error = frequency_response - self.get_des((k-1) as usize);

        // Weighted error is the error multiplied by the band's weight
        let weighted_error = error * self.get_wt((k-1) as usize);

        weighted_error
    }
}

fn generate_grid(bands: &[Band], delta_f: f32, neg: bool) -> (Vec<f32>, Vec<usize>) {
    let mut grid_buffer = vec![];
    let mut band_buffer = vec![];

    for b in 0..bands.len() {
        let band_coefficients = band_coefficients(&bands, delta_f, b, neg);
        for coeff in band_coefficients {
            grid_buffer.push(coeff);
            band_buffer.push(b);
        }
    }

    (grid_buffer, band_buffer)
}

fn band_coefficients(bands: &[Band], delta_f: f32, l_band: usize, neg: bool) -> Vec<f32> {
    let mut coefficients = vec![];
    if l_band == 0 && neg && bands[0].lower_edge < delta_f {
        coefficients.push(delta_f);
    } else {
        coefficients.push(bands[l_band].lower_edge);
    }

    let mut last_frequency = coefficients[0];

    while last_frequency + delta_f <= bands[l_band].upper_edge {
        let current_frequency = last_frequency + delta_f;
        coefficients.push(current_frequency);
        last_frequency = current_frequency;
    }

    let coeffs_len = coefficients.len();
    coefficients[coeffs_len-1] = bands[l_band].upper_edge;

    coefficients
}

fn generate_desired_values(bands: &[Band], j_type: FilterType, grid: &Vec<f32>, grid_band: &Vec<usize>) -> [f32; 1045] {
    let mut des = [0.0f32; 1045];

    for j in 0..grid.len() {
        // Which band does this grid index correspond to?
        let band_number = grid_band[j];

        let f_up = bands[band_number].upper_edge;
        let new_grid_value = grid[j];

        if new_grid_value < f_up {
            des[j] = desired_value_at(&bands, grid[j], band_number, j_type);
        } else {
            des[j] = desired_value_at(&bands, f_up, band_number, j_type);
        }
    }

    des
}

fn generate_weights(bands: &[Band], j_type: FilterType, grid: &Vec<f32>, grid_band: &Vec<usize>) -> [f32; 1045] {
    let mut wt = [0.0f32; 1045];

    for j in 0..grid.len() {
        // Which band does this grid index correspond to?
        let band_number = grid_band[j];

        let f_up = bands[band_number].upper_edge;
        let new_grid_value = grid[j];

        if new_grid_value < f_up {
            wt[j] = weight_at(&bands, grid[j], band_number, j_type);
        } else {
            wt[j] = weight_at(&bands, f_up, band_number, j_type);
        }
    }

    wt
}

// Used internally to generate `des`
// Function to calculate the desired magnitude response as a function of frequency.
// An arbitrary function of frequency can be approximated if the user replaces this function
// with the appropriate code to evaluate the ideal magnitude.
// Note that the parameter `freq` is the value of **normalized** frequency needed for evaluation.
fn desired_value_at(bands: &[Band], freq: f32, l_band: usize, j_type: FilterType) -> f32 {
    match j_type {
        FilterType::Differentiator => {
            bands[l_band].desired_value * freq
        },
        _ => {
            bands[l_band].desired_value
        },
    }
}

// Used internally to generate `wt`
// Function to calculate the weight function as a function of frequency.
// Similar to the function `eff`, this function can be replaced by a user-written function
// to calculate any desired weighting function.
fn weight_at(bands: &[Band], freq: f32, l_band: usize, j_type: FilterType) -> f32 {
    match j_type {
        FilterType::Differentiator => {
            if bands[l_band].desired_value < 0.0001 {
                bands[l_band].weight
            } else {
                bands[l_band].weight / freq
            }
        },
        _ => {
            bands[l_band].weight
        }
    }
}