use super::DenseGrid;

pub struct ExtremalFrequencies {
    dense_grid_indexes: [i64; 66],
    // Probably can be rolled into a Vec. Note there are (num_coefficients + 1) frequencies.
    // TODO: actually it looks like Parks-McClellan uses this last index as a scratch area. wtf?
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

    pub fn get_iext(&self) -> [i64; 66] {
        self.dense_grid_indexes
    }

    pub fn get_grid_index(&self, grid_index: usize) -> i64 {
        self.dense_grid_indexes[grid_index]
    }

    pub fn set_grid_index(&mut self, grid_index: usize, value: i64) {
        self.dense_grid_indexes[grid_index] = value;
    }

    pub fn shift_grid_indexes_left(&mut self) {
        let last_coefficient_index = self.num_coefficients + 1;

        for j in 1..last_coefficient_index+1 {
            let temp_freq = self.get_grid_index(j);
            self.set_grid_index(j-1, temp_freq);
        }

        // TODO: there's some weirdness here with the "last_grid_index"
        //       (which seems to be a scratch value)
    }

    pub fn shift_grid_indexes_right(&mut self, last_grid_index: i64) {
        self.debug_print_dense_grid_indexes("beforeshift");

        for j in 1..(self.num_coefficients +1) {
            let nzzmj = self.num_coefficients + 2 - j;
            let nzmj = self.num_coefficients + 1 - j;
            let temp_freq = self.get_grid_index(nzmj-1);
            self.set_grid_index(nzzmj-1, temp_freq);
        }

        // TODO: there's some weirdness here with the "last_grid_index"
        //       (which seems to be a scratch value)
        self.set_grid_index(0, last_grid_index);

        self.debug_print_dense_grid_indexes("aftershift");
        println!();
    }

    fn debug_print_dense_grid_indexes(&self, prefix: &str) {
        print!("{}: ", prefix);
        for x in self.dense_grid_indexes.iter() {
            print!("{}, ", *x);
        }
        println!();
    }
}