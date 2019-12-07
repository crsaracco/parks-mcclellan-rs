use super::Band;

const PI: f64 = std::f64::consts::PI;
const PI2: f64 = PI * 2.0;

pub struct DenseGrid {
    buffer: Vec<f32>,
    del_f: f32,
}

impl DenseGrid {
    pub fn new(bands: &Vec<Band>, filter_length: usize, grid_density: usize, num_coefficients: usize, neg: i32) -> Self {
        let del_f = 0.5 / ((grid_density as f32) * (num_coefficients as f32));
        let buffer = generate_buffer(bands, del_f, neg);

        Self {
            buffer,
            del_f,
        }
    }

    pub fn get(&self, index: usize) -> f32 {
        self.buffer[index]
    }

    pub fn del_f(&self) -> f32 {
        self.del_f
    }

    // Function to evaluate the frequency response using the lagrange interpolation formula
    // in the barycentric form
    pub fn gee(&self, zeroth_value_override: Option<f32>, x: &[f64; 66], y: &[f64; 66], ad: &[f64; 66], k: i64, n: usize) -> f64 {
        let mut p = 0.0;
        let mut d = 0.0;

        let mut xf = if k-1 == 0 && zeroth_value_override.is_some() {
            zeroth_value_override.unwrap() as f64
        } else {
            self.get((k-1) as usize) as f64
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
fn generate_buffer(bands: &Vec<Band>, del_f: f32, neg: i32) -> Vec<f32> {
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