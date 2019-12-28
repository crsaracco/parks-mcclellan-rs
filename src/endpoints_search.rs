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
    left_nu: i32,
    right_nu: i32,
    k1: i64,
    knz: i64,
    deviation: f64,
    extremal_frequencies_changed: bool,
    extremal_frequencies: &mut ExtremalFrequencies,
) -> EndpointSearchResult{
    let last_coefficient_index = num_coefficients + 1;
    let deviation = deviation * 1.00001;

    let left_sign = -left_nu;
    let right_sign = -right_nu;

    // Restrict the range of endpoint frequencies to search in.
    // The highest left-endpoint is min(old_k1, new_k1)
    // The lowest right-endpoint is max(old_knz, new_knz)
    let left_k_max = k1.min(extremal_frequencies.get_grid_index(0));
    let right_k_min = knz.max(extremal_frequencies.get_grid_index(last_coefficient_index - 1));

    let left = search_for_lower_endpoint(num_coefficients, grid, x, y, ad, left_k_max, deviation, left_sign);
    let right = search_for_upper_endpoint(num_coefficients, grid, x, y, ad, right_k_min, deviation, right_sign);

    match (left, right) {
        (Some(left), Some(right)) => {
            let left_error = left.1;
            let right_error = right.1;
            if left_error > right_error {
                // Left endpoint wins.
                let left_index = left.0;
                extremal_frequencies.shift_grid_indexes_right(left_index - 1);
                return EndpointSearchResult::KeepIteratingRemez;
            } else {
                // Right endpoint wins.
                let right_index = right.0;
                extremal_frequencies.shift_grid_indexes_left(right_index + 1);
                return EndpointSearchResult::KeepIteratingRemez;
            }
        },
        (Some(left), None) => {
            // Left endpoint wins.
            let left_index = left.0;
            extremal_frequencies.shift_grid_indexes_right(left_index - 1);
            return EndpointSearchResult::KeepIteratingRemez;
        },
        (None, Some(right)) => {
            // Right endpoint wins.
            let right_index = right.0;
            extremal_frequencies.shift_grid_indexes_left(right_index + 1);
            return EndpointSearchResult::KeepIteratingRemez;
        },
        (None, None) => {
            // Neither endpoints were sufficient. If we've changed some other extremal frequency
            // already, just keep iterating the Remez exchange algorithm. Otherwise, stop.
            if extremal_frequencies_changed {
                return EndpointSearchResult::KeepIteratingRemez;
            } else {
                return EndpointSearchResult::StopIteratingRemez;
            }
        },
    };
}

// TODO: You can probably combine these two search functions.
//       They're mostly doing the same thing, except in opposite directions.
fn search_for_lower_endpoint(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    k_max: i64,
    deviation: f64,
    sign: i32,
) -> Option<(/* k */ i64, /* error */ f64)> {
    let last_coefficient_index = num_coefficients + 1;
    let mut largest_error = deviation;
    let mut found = false;

    // First, try to find a weighted error that's above `deviation`.
    // If we find such an error, then there is a local max somewhere that we can use
    // as an endpoint frequency.
    let mut k = 1;
    while k <= k_max {
        let error = calculate_err(grid, x, y, ad, k, last_coefficient_index);
        let sign_corrected_error = (sign as f64) * (error as f64);
        if sign_corrected_error > deviation {
            largest_error = sign_corrected_error;
            found = true;
            break;
        }
        k += 1;
    }

    // If we didn't find such an error, then there's no possible endpoint to use here.
    if !found {
        return None;
    }

    // If we did find such an error, continue searching for the first local max.
    // As soon as the weighted error starts decreasing, we've found our local max.
    k += 1;
    while k <= k_max {
        let error = calculate_err(grid, x, y, ad, k, last_coefficient_index);
        let sign_corrected_error = (sign as f64) * (error as f64);
        if sign_corrected_error <= largest_error {
            return Some((k, largest_error));
        } else {
            largest_error = sign_corrected_error;
        }
        k += 1;
    }

    // We didn't find a local max (the error kept increasing all the way to k_max)
    // Just return k_max as the new endpoint.
    return Some((k_max, largest_error));
}

fn search_for_upper_endpoint(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    k_min: i64,
    deviation: f64,
    sign: i32,
) -> Option<(/* k */ i64, /* error */ f64)> {
    let last_coefficient_index = num_coefficients + 1;
    let mut largest_error = deviation;
    let mut found = false;

    // First, try to find a weighted error that's above `deviation`.
    // If we find such an error, then there is a local max somewhere that we can use
    // as an endpoint frequency.
    let mut k = grid.n_grid() as i64;
    while k >= k_min {
        let error = calculate_err(grid, x, y, ad, k, last_coefficient_index);
        let sign_corrected_error = (sign as f64) * (error as f64);
        if sign_corrected_error > deviation {
            largest_error = sign_corrected_error;
            found = true;
            break;
        }
        k -= 1;
    }

    // If we didn't find such an error, then there's no possible endpoint to use here.
    if !found {
        return None;
    }

    // If we did find such an error, continue searching for the first local max.
    // As soon as the weighted error starts decreasing, we've found our local max.
    k -= 1;
    while k >= k_min {
        let error = calculate_err(grid, x, y, ad, k, last_coefficient_index);
        let sign_corrected_error = (sign as f64) * (error as f64);
        if sign_corrected_error <= largest_error {
            return Some((k, largest_error));
        } else {
            largest_error = sign_corrected_error;
        }
        k -= 1;
    }

    // We didn't find a local max (the error kept increasing all the way to k_min)
    // Just return k_min as the new endpoint.
    return Some((k_min, largest_error));
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