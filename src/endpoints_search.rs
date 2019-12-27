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

    let left_nut = -nu;
    let right_nut = -*nut;

    // Restrict the range of endpoint frequencies to search in.
    // The highest left-endpoint is min(old_k1, new_k1)
    // The lowest right-endpoint is max(old_knz, new_knz)
    let kup = k1.min(extremal_frequencies.get_grid_index(0));
    let knz = knz.max(extremal_frequencies.get_grid_index(last_coefficient_index - 1));

    let (my_index_left, my_weighted_err_left) = my_search_lower(num_coefficients, grid, x, y, ad, kup, deviation*1.00001, left_nut);
    let (my_index_right, my_weighted_err_right) = my_search_upper(num_coefficients, grid, x, y, ad, knz, deviation*1.00001, right_nut);

    match (my_weighted_err_left, my_weighted_err_right) {
        (Some(left), Some(right)) => {
            println!("left: {} | right: {}", left, right);
            if left>right {
                extremal_frequencies.set_grid_index(last_coefficient_index,  my_index_left.unwrap() - 1);
                let k1 = extremal_frequencies.get_grid_index(last_coefficient_index);
                extremal_frequencies.shift_grid_indexes_right(k1);
                return EndpointSearchResult::KeepIteratingRemez;
            } else {
                extremal_frequencies.set_grid_index(last_coefficient_index, my_index_right.unwrap() + 1);
                extremal_frequencies.shift_grid_indexes_left();
                return EndpointSearchResult::KeepIteratingRemez;
            }
        },
        (Some(_), None) => {
            println!("LEFT IS SOME");
            extremal_frequencies.set_grid_index(last_coefficient_index,  my_index_left.unwrap() - 1);
            let k1 = extremal_frequencies.get_grid_index(last_coefficient_index);
            extremal_frequencies.shift_grid_indexes_right(k1);
            return EndpointSearchResult::KeepIteratingRemez;
        },
        (None, Some(_)) => {
            println!("RIGHT IS SOME");
            extremal_frequencies.set_grid_index(last_coefficient_index, my_index_right.unwrap() + 1);
            extremal_frequencies.shift_grid_indexes_left();
            return EndpointSearchResult::KeepIteratingRemez;
        },
        (None, None) => {
            println!("ALL IS NONE");
            if extremal_frequencies_changed {
                return EndpointSearchResult::KeepIteratingRemez;
            } else {
                return EndpointSearchResult::StopIteratingRemez;
            }
        },
    };
}

fn my_search_lower(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    kup: i64,
    deviation: f64,
    nut: i32,
) -> (/* index */ Option<i64>, /* weighted_err */ Option<f64>) {
    let last_coefficient_index = num_coefficients + 1;
    let mut index = 0;
    let mut comp = deviation;
    let mut found = false;

    for k in 1..kup {
        let err = calculate_err(grid, x, y, ad, k, last_coefficient_index);
        let weighted_err = (nut as f64) * (err as f64);
        if !found {
            if weighted_err > comp {
                index = k;
                comp = weighted_err;
                found = true;
            }
        } else {
            if weighted_err > comp {
                index = k;
                comp = weighted_err;
            } else {
                index = k;
                return (Some(index), Some(comp));
            }
        }
    }

    if !found {
        return (None, None);
    } else {
        return (Some(index+1), Some(comp));
    }
}

fn my_search_upper(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    klow: i64,
    comp: f64,
    nut: i32,
) -> (/* ell */ Option<i64>, /* comp */ Option<f64>) {
    let last_coefficient_index = num_coefficients + 1;
    let mut index = 0;
    let mut comp = comp;
    let mut found = false;

    for k in ((klow)..=(grid.n_grid() as i64)).rev() {
        println!("k: {}", k);
        let err = calculate_err(grid, x, y, ad, k, last_coefficient_index);
        let weighted_err = (nut as f64) * (err as f64);
        if !found {
            if weighted_err > comp {
                index = k;
                comp = weighted_err;
                found = true;
            }
        } else {
            if weighted_err > comp {
                index = k;
                comp = weighted_err;
            } else {
                index = k;
                return (Some(index), Some(comp));
            }
        }
    }

    if !found {
        return (None, None)
    } else {
        return (Some(index), Some(comp));
    }
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