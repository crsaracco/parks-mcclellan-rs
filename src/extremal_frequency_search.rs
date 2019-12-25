use super::DenseGrid;
use super::ExtremalFrequencies;

pub fn find_extremal_frequencies(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    deviation: f64,

    nut: &mut i32,
    comp: &mut f64,
    extremal_frequencies: &mut ExtremalFrequencies,
) -> bool /* extremal_frequencies_changed */ {
    let coefficient_off_end_index = num_coefficients + 2;
    let mut klow = 0;

    let mut extremal_frequencies_changed = false;

    for j in 1..coefficient_off_end_index {
        *nut = -*nut;
        *comp = deviation;

        search_for_frequency(
            num_coefficients,
            grid,
            x,
            y,
            ad,
            deviation,
            j,
            *nut,
            comp,
            &mut klow,
            &mut extremal_frequencies_changed,
            extremal_frequencies
        );
    }

    extremal_frequencies_changed
}

fn search_for_frequency(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    deviation: f64,
    j: usize,
    nut: i32,

    comp: &mut f64,
    klow: &mut i64,
    extremal_frequencies_changed: &mut bool,
    extremal_frequencies: &mut ExtremalFrequencies,
) {
    let last_coefficient_index = num_coefficients + 1;
    let kup = extremal_frequencies.get_grid_index(j);

    // Start the search at k+1. If this frequency is "suitable", keep searching upwards (towards
    // kup) for a local maximum.
    let mut ell = extremal_frequencies.get_grid_index(j - 1) + 1;
    if ell < kup { // Careful not to overflow...
        let err = calculate_err(grid, &x, &y, &ad, ell, last_coefficient_index);
        if deviation < (nut as f64) * (err as f64) {
            search_upwards(num_coefficients, grid, x, y, ad, nut, err, ell, kup, j,
                           klow, comp, extremal_frequencies_changed, extremal_frequencies);
            return; // ON TO THE NEXT EXTREMA
        }
    }

    // If that frequency wasn't "suitable", check k-1. If that is "suitable", keep searching
    // downwards (towards klow) for a local maximum.
    ell = extremal_frequencies.get_grid_index(j - 1) - 1;
    if ell > *klow { // Careful not to underflow...
        let err = calculate_err(grid, &x, &y, &ad, ell, last_coefficient_index);
        if (nut as f64) * (err as f64) - *comp > 0.0 {
            // There is a local max of error curve. Keep searching k-2, k-3, ..., klow
            // until the local max is found.
            search_downwards(num_coefficients, grid, x, y, ad, nut, err, ell, j,
                             klow, comp, extremal_frequencies_changed, extremal_frequencies);
            return; // ON TO THE NEXT EXTREMA
        }
    }

    // We didn't find a "suitable" frequency at either k+1 or k-1.
    // If we already changed an extrema this round, don't worry about it. Just keep this extremal
    // frequency as-is, and move on to the next extrema.
    if *extremal_frequencies_changed {
        *klow = extremal_frequencies.get_grid_index(j - 1);
        return; // ON TO THE NEXT EXTREMA
    }

    // We're getting a little desperate to change some frequency this round,
    // so search **ALL** the frequencies!

    // Start by searching downwards (starting at k-2, going towards klow)
    ell = extremal_frequencies.get_grid_index(j - 1) - 2;
    while ell > *klow {
        let err = calculate_err(grid, &x, &y, &ad, ell, last_coefficient_index);
        if (nut as f64) * (err as f64) - *comp > 0.0 {
            // There is a local max of error curve. Keep searching k-2, k-3, ..., klow
            // until the local max is found.
            search_downwards(num_coefficients, grid, x, y, ad, nut, err, ell, j,
                             klow, comp, extremal_frequencies_changed, extremal_frequencies);
            return; // ON TO THE NEXT EXTREMA
        }
        ell -= 1;
    }

    // Still didn't find anything. Search upwards (starting at k+2, going towards kup)
    ell = extremal_frequencies.get_grid_index(j - 1) + 2;
    while ell < kup {
        let err = calculate_err(grid, &x, &y, &ad, ell, last_coefficient_index);
        if (nut as f64) * (err as f64) - *comp > 0.0 {
            // Found a local maximum!
            search_upwards(num_coefficients, grid, x, y, ad, nut, err, ell, kup, j,
                           klow, comp, extremal_frequencies_changed, extremal_frequencies);
            return; // ON TO THE NEXT EXTREMA
        }
        ell += 1;
    }

    // We didn't find **ANY** local max in this entire range. We can't do anything else, so
    // continue on to the next frequency.
    *klow = extremal_frequencies.get_grid_index(j - 1);
}

// There is a local maximum of the error curve where the signed error is greater
// than the present deviation. Continue the search at k-2, k-3, ... klow
// until this local maximum is found.
fn search_downwards(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    nut: i32,
    err: f32,
    ell: i64,
    j: usize,
    klow: &mut i64,
    comp: &mut f64,
    extremal_frequencies_changed: &mut bool,
    extremal_frequencies: &mut ExtremalFrequencies,
) {
    let last_coefficient_index = num_coefficients + 1;
    let mut local_max_found = false;
    *comp = (nut as f64) * (err as f64);

    for k in ((*klow+1)..(ell)).rev() {
        let err = calculate_err(grid, &x, &y, &ad, k, last_coefficient_index);
        let dtemp = (nut as f64) * (err as f64) - *comp;
        if dtemp <= 0.0 {
            // Local max found!
            local_max_found = true;
            update_loop_variables_downwards_loop(j-1, k+1, klow, extremal_frequencies_changed, extremal_frequencies);
            *extremal_frequencies_changed = true;
            break;
        }
        *comp = (nut as f64) * (err as f64);
    }

    if !local_max_found {
        update_loop_variables_downwards_loop(j-1, *klow+1, klow, extremal_frequencies_changed, extremal_frequencies);
    }
}

// There is a local maximum of the error curve where the signed error is greater
// than the present deviation. Continue the search at k+2, k+3, ... kup
// until this local maximum is found.
fn search_upwards(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    nut: i32,
    err: f32,
    ell: i64,
    kup: i64,
    j: usize,
    klow: &mut i64,
    comp: &mut f64,
    extremal_frequencies_changed: &mut bool,
    extremal_frequencies: &mut ExtremalFrequencies,
) {
    let last_coefficient_index = num_coefficients + 1;
    let mut local_max_found = false;
    *comp = (nut as f64) * (err as f64);

    for k in (ell+1)..kup {
        let err = calculate_err(grid, &x, &y, &ad, k, last_coefficient_index);
        if (nut as f64) * (err as f64) - *comp <= 0.0 {
            // Error is starting to fall. Record the current best.
            // Go to next extremal frequency.
            update_loop_variables_upwards_loop(j-1, k-1, klow, extremal_frequencies_changed, extremal_frequencies);
            local_max_found = true;
            break; // ON TO THE NEXT EXTREMA
        }
        *comp = (nut as f64) * (err as f64);
    }

    if !local_max_found {
        update_loop_variables_upwards_loop(j-1, kup-1, klow, extremal_frequencies_changed, extremal_frequencies);
    }
}

// Change the extremal frequency, update klow and kup for the next iteration, and note
// that an extremal frequency has been changed.
fn update_loop_variables_upwards_loop(
    grid_index: usize,
    grid_value: i64,
    klow: &mut i64,
    extremal_frequencies_changed: &mut bool,
    extremal_frequencies: &mut ExtremalFrequencies,
) {
    extremal_frequencies.set_grid_index(grid_index, grid_value);
    *klow = grid_value;
    *extremal_frequencies_changed = true;
}

fn update_loop_variables_downwards_loop(
    grid_index: usize,
    grid_value: i64,
    klow: &mut i64,
    extremal_frequencies_changed: &mut bool,
    extremal_frequencies: &mut ExtremalFrequencies,
) {
    let new_klow = extremal_frequencies.get_grid_index(grid_index);
    extremal_frequencies.set_grid_index(grid_index, grid_value);
    *klow = new_klow;
    *extremal_frequencies_changed = true;
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