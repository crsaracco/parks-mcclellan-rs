use super::DenseGrid;

pub struct ExtremalFrequencies {
    // TODO: Probably can be rolled into a Vec. Note there are (num_coefficients + 1) frequencies.
    dense_grid_indexes: [i64; 66],
    num_coefficients: usize,
}

impl ExtremalFrequencies {
    pub fn new(num_coefficients: usize) -> Self {
        Self {
            dense_grid_indexes: [0; 66],
            num_coefficients,
        }
    }

    // Initialize the current indexes to be equally spaced among the available grid indexes.
    pub fn initialize_guess(&mut self, grid: &DenseGrid) {
        let temp = ((grid.n_grid()-1) as f32) / (self.num_coefficients as f32);
        for j in 0..self.num_coefficients {
            self.dense_grid_indexes[j] = ((j as f32) * temp + 1.0) as i64;
        }
        // The last frequency is the end of the grid.
        self.dense_grid_indexes[self.num_coefficients] = grid.n_grid() as i64;
    }

    pub fn get_grid_index(&self, grid_index: usize) -> i64 {
        self.dense_grid_indexes[grid_index]
    }

    pub fn set_grid_index(&mut self, grid_index: usize, value: i64) {
        self.dense_grid_indexes[grid_index] = value;
    }

    pub fn shift_grid_indexes_left(&mut self, last_grid_index: i64) {
        for j in 0..=self.num_coefficients {
            let temp_freq = self.get_grid_index(j+1);
            self.set_grid_index(j, temp_freq);
        }
        self.dense_grid_indexes[self.num_coefficients] = last_grid_index;
    }

    pub fn shift_grid_indexes_right(&mut self, first_grid_index: i64) {
        for j in 0..self.num_coefficients {
            let temp_freq = self.get_grid_index(self.num_coefficients - j - 1);
            self.set_grid_index(self.num_coefficients - j, temp_freq);
        }
        self.dense_grid_indexes[0] = first_grid_index;
    }
}