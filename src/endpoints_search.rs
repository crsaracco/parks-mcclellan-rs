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
    let left_nut = -nu;
    let right_nut = -*nut;

    let mut left_endpoint_wins = false;
    let mut left_endpoint_index: Option<i64> = None;
    let mut left_weighted_err = 0.0;
    let mut right_endpoint_wins = false;
    let mut right_endpoint_index: Option<i64> = None;
    let mut right_weighted_err = 0.0;

    // Restrict the range of endpoint frequencies to search in.
    // The highest left-endpoint is min(old_k1, new_k1)
    // The lowest right-endpoint is max(old_knz, new_knz)
    let k1 = k1.min(extremal_frequencies.get_grid_index(0));
    let knz = knz.max(extremal_frequencies.get_grid_index(last_coefficient_index - 1));

    // Some other variables
    nut1 = *nut;
    *nut = -nu;
    let mut ell = 1;
    let kup = k1;
    *comp = *comp * 1.00001;

    let (my_index_left, my_weighted_err_left) = my_search_lower(num_coefficients, grid, x, y, ad, kup, deviation, left_nut);
    let (my_index_right, my_weighted_err_right) = my_search_upper(num_coefficients, grid, x, y, ad, knz, deviation, right_nut);

    while ell < kup {
        let err = calculate_err(grid, x, y, ad, ell, last_coefficient_index);
        if (*nut as f64) * (err as f64) - *comp > 0.0 {
            *comp = (*nut as f64) * (err as f64);

            let p = search_for_lower_endpoint(num_coefficients, grid, x, y, ad, kup, *comp, *nut, ell);
            ell = p.0;
            *comp = p.1;

            left_endpoint_wins = true;
            left_endpoint_index = Some(ell);
            left_weighted_err = *comp;

//            assert_eq!(my_index_left, left_endpoint_index.unwrap());
//            assert_eq!(my_weighted_err_left, left_weighted_err);

            if *comp > y1 {
                y1 = *comp;
            }

            break;
        }
        ell += 1;
    }

    //assert_eq!(left_endpoint_index, my_index_left);


    ell = grid.n_grid() as i64;
    *nut = -nut1;
    *comp = y1 * 1.00001;

    let (my_index_right, my_weighted_err_right) = my_search_upper(num_coefficients, grid, x, y, ad, knz, deviation, right_nut);

    while ell > knz {
        let err = calculate_err(grid, x, y, ad, ell, last_coefficient_index);
        if (*nut as f64) * (err as f64) - *comp > 0.0 {
            // TODO: should comp be set here???

            ell = search_for_upper_endpoint(num_coefficients, grid, x, y, ad, knz, *comp, *nut, ell);

            right_endpoint_wins = true;
            left_endpoint_wins = false;
            right_endpoint_index = Some(ell);

            assert_eq!(my_index_right, right_endpoint_index.unwrap());

            break;
        }
        ell -= 1;
    }


    // If one of the endpoints "won", save it to the extremal frequencies.
    if left_endpoint_wins {
        extremal_frequencies.set_grid_index(last_coefficient_index,  left_endpoint_index.unwrap() - 1);
        let k1 = extremal_frequencies.get_grid_index(last_coefficient_index);
        extremal_frequencies.shift_grid_indexes_right(k1);
        return EndpointSearchResult::KeepIteratingRemez;
    } else if right_endpoint_wins {
        extremal_frequencies.set_grid_index(last_coefficient_index, right_endpoint_index.unwrap() + 1);
        extremal_frequencies.shift_grid_indexes_left();
        return EndpointSearchResult::KeepIteratingRemez;
    }

    // If none of the endpoints "won", figure out if we should keep iterating.
    if extremal_frequencies_changed {
        return EndpointSearchResult::KeepIteratingRemez;
    } else {
        return EndpointSearchResult::StopIteratingRemez;
    }
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
) -> (/* index */ i64, /* weighted_err */ f64) {
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
                return (index, comp);
            }
        }
    }

    (index+1, comp)
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
) -> (/* ell */ i64, /* comp */ f64) {
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
                return (index, comp);
            }
        }
    }

    (index, comp)
}

fn search_for_lower_endpoint(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    kup: i64,
    comp: f64,
    nut: i32,
    ell: i64,
) -> (/* ell */ i64, /* comp */ f64) {
    let last_coefficient_index = num_coefficients + 1;
    let mut comp = comp;
    let mut ell = ell;

    loop {
        ell += 1;
        if ell >= kup {
            break;
        }
        let err = calculate_err(grid, x, y, ad, ell, last_coefficient_index);
        if (nut as f64) * (err as f64) - comp <= 0.0 {
            break;
        }
        comp = (nut as f64) * (err as f64);
    }

    (ell, comp)
}

fn search_for_upper_endpoint(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    klow: i64,
    comp: f64,
    nut: i32,
    ell: i64,
) -> /* ell */ i64 {
    let last_coefficient_index = num_coefficients + 1;
    let mut comp = comp;
    let mut ell = ell;

    while ell > klow {
        let err = calculate_err(grid, x, y, ad, ell, last_coefficient_index);
        if (nut as f64) * (err as f64) - comp <= 0.0 {
            break;
        }
        comp = (nut as f64) * (err as f64);
        ell -= 1;
    }

    ell
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