use super::Band;
use super::JType;

const PI: f64 = std::f64::consts::PI;
const PI2: f64 = PI * 2.0;

pub struct DenseGrid {
    grid: Vec<f32>,
    des: [f32; 1045], // num elements?
    wt: [f32; 1045], // num elements?
    del_f: f32,
    n_grid: usize,
}

impl DenseGrid {
    pub fn new(bands: &Vec<Band>, j_type: JType, filter_length: usize, grid_density: usize, num_coefficients: usize, neg: i32, n_odd: i32) -> Self {
        let del_f = 0.5 / ((grid_density as f32) * (num_coefficients as f32));
        let grid = generate_grid(bands, del_f, neg);
        
        let mut grid = Self {
            grid,
            des: [0.0; 1045], // Will be filled in after
            wt: [0.0; 1045], // Will be filled in after
            del_f,
            n_grid: 0, // Will be filled in after
        };

        let (des, wt, n_grid) = grid.generate_des_wt(bands, j_type, neg, n_odd);

        grid.des = des;
        grid.wt = wt;
        grid.n_grid = n_grid;

        grid
    }

    pub fn get_grid(&self, index: usize) -> f32 {
        self.grid[index]
    }

    pub fn get_des(&self, index: usize) -> f32 {
        self.des[index]
    }

    pub fn set_des(&mut self, index: usize, value: f32) {
        self.des[index] = value;
    }

    pub fn get_wt(&self, index: usize) -> f32 {
        self.wt[index]
    }

    pub fn set_wt(&mut self, index: usize, value: f32) {
        self.wt[index] = value;
    }

    pub fn del_f(&self) -> f32 {
        self.del_f
    }

    pub fn n_grid(&self) -> usize {
        self.n_grid
    }

    // Used internally to generate `des` and `wt`
    fn generate_des_wt(&self, bands: &Vec<Band>, j_type: JType, neg: i32, n_odd: i32) -> ([f32; 1045], [f32; 1045], usize) {
        let mut des = [0.0f32; 1045];
        let mut wt = [0.0f32; 1045];

        let mut j = 1;
        let mut l_band = 1;
        loop {
            let f_up = bands[l_band-1].upper_edge;
            let new_grid_value = self.get_grid(j-1) + self.del_f();

            if new_grid_value > f_up {
                l_band = l_band + 1;
            }

            if new_grid_value <= f_up {
                des[j-1] = eff(&bands, self.get_grid(j-1), l_band, j_type);
                wt[j-1] = wate(&bands, self.get_grid(j-1), l_band, j_type);
            } else {
                des[j-1] = eff(&bands, f_up, l_band-1, j_type);
                wt[j-1] = wate(&bands, f_up, l_band-1, j_type);
            }

            if new_grid_value > f_up && l_band > bands.len() {
                break;
            }
            j += 1;
        }
        j += 1;

        let mut n_grid = j-1;

        if neg == n_odd {
            if self.get_grid(n_grid-1) > (0.5 - self.del_f()) {
                n_grid = n_grid - 1;
            }
        }

        (des, wt, n_grid)
    }

    // Function to evaluate the frequency response using the lagrange interpolation formula
    // in the barycentric form
    pub fn gee(&self, zeroth_value_override: Option<f32>, x: &[f64; 66], y: &[f64; 66], ad: &[f64; 66], k: i64, n: usize) -> f64 {
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
}

// Used internally to generate the grid.
fn generate_grid(bands: &Vec<Band>, del_f: f32, neg: i32) -> Vec<f32> {
    let mut new_buffer = vec![];
    for b in 0..bands.len() {
        let band_coefficients = band_coefficients(&bands, del_f, b, neg);
        for coeff in band_coefficients {
            new_buffer.push(coeff);
        }
    }

    new_buffer
}

// Used internally to generate the grid.
fn band_coefficients(bands: &Vec<Band>, del_f: f32, l_band: usize, neg: i32) -> Vec<f32> {
    let mut coefficients = vec![];
    if l_band == 0 && neg != 0 && bands[0].lower_edge < del_f {
        coefficients.push(del_f);
    } else {
        coefficients.push(bands[l_band].lower_edge);
    }

    let mut last_frequency = coefficients[0];

    while last_frequency + del_f <= bands[l_band].upper_edge {
        let current_frequency = last_frequency + del_f;
        coefficients.push(current_frequency);
        last_frequency = current_frequency;
    }

    let coeffs_len = coefficients.len();
    coefficients[coeffs_len-1] = bands[l_band].upper_edge;

    coefficients
}

// Used internally to generate `des`
// Function to calculate the desired magnitude response as a function of frequency.
// An arbitrary function of frequency can be approximated if the user replaces this function
// with the appropriate code to evaluate the ideal magnitude.
// Note that the parameter `freq` is the value of **normalized** frequency needed for evaluation.
fn eff(bands: &Vec<Band>, freq: f32, l_band: usize, j_type: JType) -> f32 {
    match j_type {
        JType::Differentiator => {
            bands[l_band-1].desired_value * freq
        },
        _ => {
            bands[l_band-1].desired_value
        },
    }
}

// Used internally to generate `wt`
// Function to calculate the weight function as a function of frequency.
// Similar to the function `eff`, this function can be replaced by a user-written function
// to calculate any desired weighting function.
fn wate(bands: &Vec<Band>, freq: f32, l_band: usize, j_type: JType) -> f32 {
    match j_type {
        JType::Differentiator => {
            if bands[l_band-1].desired_value < 0.0001 {
                bands[l_band-1].weight
            } else {
                bands[l_band-1].weight / freq
            }
        },
        _ => {
            bands[l_band-1].weight
        }
    }
}