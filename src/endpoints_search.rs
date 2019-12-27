use crate::DenseGrid;
use super::ExtremalFrequencies;

pub enum EndpointSearchResult {
    KeepIteratingRemez,
    StopIteratingRemez,
}

pub fn endpoints_search(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    nu: i32,
    k1: i64,
    knz: i64,
    deviation: f64,

    comp: &mut f64,
    nut: &mut i32,
    extremal_frequencies_changed: bool,
    extremal_frequencies: &mut ExtremalFrequencies,
) -> EndpointSearchResult{
    let last_coefficient_index = num_coefficients + 1;
    let coefficient_off_end_index = num_coefficients + 2;

    let mut nut1 = 0;
    let mut k1 = k1;
    let mut knz = knz;
    let mut y1 = deviation;

    // Restrict the range of endpoint frequencies to search in.
    // The highest left-endpoint is min(old_k1, new_k1)
    // The lowest right-endpoint is max(old_knz, new_knz)
    k1 = k1.min(extremal_frequencies.get_grid_index(0));
    knz = knz.max(extremal_frequencies.get_grid_index(last_coefficient_index - 1));

    let err_to_beat = *comp * 1.00001;

    // Some other variables
    nut1 = *nut;
    *nut = -nu;
    let mut ell = 1;
    let kup = k1;
    *comp = err_to_beat;
    let mut err;

    while ell < kup {
        err = calculate_err(grid, x, y, ad, ell, last_coefficient_index);
        if (*nut as f64) * (err as f64) - *comp > 0.0 {
            *comp = (*nut as f64) * (err as f64);

            search_for_lower_endpoint(num_coefficients, grid, x, y, ad, kup, comp, nut, &mut ell, extremal_frequencies);

            if *comp > y1 {
                y1 = *comp;
            }



            ell = grid.n_grid() as i64;
            *nut = -nut1;
            *comp = y1 * 1.00001;

            while ell > knz {
                err = calculate_err(grid, x, y, ad, ell, last_coefficient_index);
                if (*nut as f64) * (err as f64) - *comp > 0.0 {
                    search_for_upper_endpoint(num_coefficients, grid, x, y, ad, knz, comp, nut, &mut ell, extremal_frequencies);
                    return EndpointSearchResult::KeepIteratingRemez;
                }
                ell -= 1;
            }

            k1 = extremal_frequencies.get_grid_index(last_coefficient_index);
            extremal_frequencies.shift_grid_indexes_right(k1);
            return EndpointSearchResult::KeepIteratingRemez;


        }
        ell += 1;
    }



    ell = grid.n_grid() as i64;
    *nut = -nut1;
    *comp = y1 * 1.00001;

    while ell > knz {
        err = calculate_err(grid, x, y, ad, ell, last_coefficient_index);
        if (*nut as f64) * (err as f64) - *comp > 0.0 {
            search_for_upper_endpoint(num_coefficients, grid, x, y, ad, knz, comp, nut, &mut ell, extremal_frequencies);
            return EndpointSearchResult::KeepIteratingRemez;
        }
        ell -= 1;
    }

    if extremal_frequencies_changed {
        return EndpointSearchResult::KeepIteratingRemez;
    }
    return EndpointSearchResult::StopIteratingRemez;


}

fn search_for_lower_endpoint(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    kup: i64,

    comp: &mut f64,
    nut: &mut i32,
    ell: &mut i64,
    extremal_frequencies: &mut ExtremalFrequencies,
) {
    let last_coefficient_index = num_coefficients + 1;

    loop {
        *ell += 1;
        if *ell >= kup {
            break;
        }
        let err = calculate_err(grid, x, y, ad, *ell, last_coefficient_index);
        if (*nut as f64) * (err as f64) - *comp <= 0.0 {
            break;
        }
        *comp = (*nut as f64) * (err as f64);
    }

    extremal_frequencies.set_grid_index(last_coefficient_index,  *ell - 1);
}

fn search_for_upper_endpoint(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    klow: i64,

    comp: &mut f64,
    nut: &mut i32,
    ell: &mut i64,
    extremal_frequencies: &mut ExtremalFrequencies,
) {
    let last_coefficient_index = num_coefficients + 1;

    while *ell > klow {
        let err = calculate_err(grid, x, y, ad, *ell, last_coefficient_index);
        if (*nut as f64) * (err as f64) - *comp <= 0.0 {
            break;
        }
        *comp = (*nut as f64) * (err as f64);
        *ell -= 1;
    }

    extremal_frequencies.set_grid_index(last_coefficient_index, *ell + 1);
    extremal_frequencies.shift_grid_indexes_left();
}

fn calculate_err(
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    ell: i64,
    last_coefficient_index: usize,
) -> f32 {
    let err = grid.gee(None, x, y, ad, ell, last_coefficient_index) as f32;
    (err - grid.get_des((ell - 1) as usize)) * grid.get_wt((ell - 1) as usize)
}