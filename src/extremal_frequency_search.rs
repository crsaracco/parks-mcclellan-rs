use super::DenseGrid;
use super::ExtremalFrequencies;

pub fn find_nth_extremal_frequency(
    num_coefficients: usize,
    grid: &DenseGrid,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    deviation: f64,

    nut: &mut i32,
    comp: &mut f64,
    y1: &mut f64,
    klow: &mut i64,
    extremal_frequencies_changed: &mut bool,
    extremal_frequencies: &mut ExtremalFrequencies,
) {
    let last_coefficient_index = num_coefficients + 1;
    let coefficient_off_end_index = num_coefficients + 2;

    'find_extremal: for j in 1..coefficient_off_end_index {
        *nut = -*nut;
        if j == 2 {
            *y1 = *comp;
        }
        *comp = deviation;

        let kup = extremal_frequencies.get_grid_index(j);
        let mut ell = extremal_frequencies.get_grid_index(j - 1) + 1;

        if ell < kup {
            // NORMAL OPERATION
            let mut err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
            if deviation < (*nut as f64) * (err as f64) {
                // There is a local maximum of the error curve where the signed error is greater
                // than the present deviation. Continue the search at k+2, k+3, ... kup
                // until this local maximum is found.
                let mut local_max_found = false;
                *comp = (*nut as f64) * (err as f64);

                for k in (ell+1)..kup {
                    err = calculate_err(grid, None, &x, &y, &ad, k, last_coefficient_index);
                    if (*nut as f64) * (err as f64) - *comp <= 0.0 {
                        // found local maximum!
                        // Change the extremal frequency, update klow and kup for the next iteration, and note
                        // that an extremal frequency has been changed.
                        // (on to the next frequency)
                        extremal_frequencies.set_grid_index(j - 1, k - 1);
                        *klow = k - 1;
                        *extremal_frequencies_changed = true;
                        local_max_found = true;
                        break;
                    } else {
                        *comp = (*nut as f64) * (err as f64);
                    }
                }

                if !local_max_found {
                    // Change the extremal frequency, update klow and kup for the next iteration, and note
                    // that an extremal frequency has been changed.
                    // (on to the next frequency)
                    extremal_frequencies.set_grid_index(j - 1, kup - 1);
                    *klow = kup - 1;
                    *extremal_frequencies_changed = true;
                }
            } else {
                // Calculate frequency response and weighted error at index k-1
                ell = ell - 1;
                'loop_13: loop {
                    // Continue searching k-2, k-3, ..., klow
                    // for a local max where signed error is >= dev
                    ell = ell - 1;
                    if ell > *klow {
                        // Searching...
                        err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                        if (*nut as f64) * (err as f64) - *comp > 0.0 {
                            // There is a local max of error curve. Keep searching k-2, k-3, ..., klow
                            // until the local max is found.
                            let mut local_max_found = false;
                            *comp = (*nut as f64) * (err as f64);

                            for k in ((*klow)..=(ell-1)).rev() {
                                ell -= 1;
                                if ell <= *klow {
                                    break;
                                }
                                assert_eq!(k, ell);
                                err = calculate_err(grid, None, &x, &y, &ad, k, last_coefficient_index);
                                let dtemp = (*nut as f64) * (err as f64) - *comp;
                                if dtemp <= 0.0 {
                                    // Local max found!
                                    local_max_found = true;
                                    let new_klow = extremal_frequencies.get_grid_index(j - 1);
                                    extremal_frequencies.set_grid_index(j - 1, k + 1);
                                    *klow = new_klow;
                                    *extremal_frequencies_changed = true;
                                    break;
                                }
                                *comp = (*nut as f64) * (err as f64);
                            }

                            if !local_max_found {
                                let new_klow = extremal_frequencies.get_grid_index(j - 1);
                                extremal_frequencies.set_grid_index(j - 1, *klow + 1);
                                *klow = new_klow;
                                *extremal_frequencies_changed = true;
                            }
                        } else {
                            if *extremal_frequencies_changed {
                                *klow = extremal_frequencies.get_grid_index(j - 1);
                            } else {
                                // No extremal frequencies changed yet, keep searching.
                                continue 'loop_13;
                            }
                        }
                        break 'loop_13;
                    } else {
                        // k is less than klow, so we didn't find any such local max.
                        // If any extremal frequencies have changed by now, don't worry about it.
                        // (???)
                        ell = extremal_frequencies.get_grid_index(j - 1) + 1;
                        if *extremal_frequencies_changed {
                            extremal_frequencies.set_grid_index(j - 1, ell - 1);
                            *klow = ell - 1;
                            continue 'find_extremal;
                        }
                        else {
                            // We haven't changed any extremal frequencies yet.
                            // Start searching k+2, j+3, ..., kup for a local max.
                            'loop_15: loop {
                                ell += 1;
                                if ell >= kup {
                                    // We didn't find **ANY** local max in this entire range.
                                    // I guess we just give up here. Keep this extremal frequency
                                    // the same. (???)
                                    *klow = extremal_frequencies.get_grid_index(j - 1);
                                    continue 'find_extremal;
                                }
                                err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                                if (*nut as f64) * (err as f64) - *comp > 0.0 {
                                    // Found a local maximum!
                                    *comp = (*nut as f64) * (err as f64);
                                    loop {
                                        ell = ell + 1;
                                        if ell >= kup {
                                            // Reached kup; stop iterating.
                                            // Go to next extremal frequency.
                                            extremal_frequencies.set_grid_index(j - 1, ell - 1);
                                            *klow = ell - 1;
                                            *extremal_frequencies_changed = true;
                                            continue 'find_extremal;
                                        }
                                        err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                                        if (*nut as f64) * (err as f64) - *comp <= 0.0 {
                                            // Error is starting to fall. Record the current best.
                                            // Go to next extremal frequency.
                                            extremal_frequencies.set_grid_index(j - 1, ell - 1);
                                            *klow = ell - 1;
                                            *extremal_frequencies_changed = true;
                                            continue 'find_extremal;
                                        }
                                        *comp = (*nut as f64) * (err as f64);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            // Edge case: if ell is greater than or equal to kup for this extrema
            ell = ell - 1;
            'loop_03: loop {
                ell = ell - 1;
                if ell <= *klow {
                    ell = extremal_frequencies.get_grid_index(j -1) + 1;
                    if *extremal_frequencies_changed {
                        extremal_frequencies.set_grid_index(j -1, ell-1);
                        *klow = ell - 1;
                        continue 'find_extremal;
                    }
                    'loop_05: loop {
                        ell += 1;
                        if ell >= kup {
                            *klow = extremal_frequencies.get_grid_index(j -1);
                            continue 'find_extremal;
                        }
                        let mut err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                        let dtemp = (*nut as f64) * (err as f64) - *comp;
                        if dtemp <= 0.0 {
                            continue 'loop_05;
                        }
                        *comp = (*nut as f64) * (err as f64);

                        loop {
                            ell = ell + 1;
                            if ell >= kup {
                                break;
                            }
                            err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                            let dtemp = (*nut as f64) * (err as f64) - *comp;
                            if dtemp <= 0.0 {
                                break;
                            }
                            *comp = (*nut as f64) * (err as f64);
                        }

                        extremal_frequencies.set_grid_index(j -1, ell-1);
                        *klow = ell - 1;
                        *extremal_frequencies_changed = true;
                        continue 'find_extremal;
                    }
                }
                let mut err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                let dtemp = (*nut as f64) * (err as f64) - *comp;
                if dtemp > 0.0 {
                    *comp = (*nut as f64) * (err as f64);

                    loop {
                        ell -= 1;
                        if ell <= *klow {
                            break;
                        }
                        err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                        let dtemp = (*nut as f64) * (err as f64) - *comp;
                        if dtemp <= 0.0 {
                            break;
                        }
                        *comp = (*nut as f64) * (err as f64);
                    }

                    *klow = extremal_frequencies.get_grid_index(j -1);
                    extremal_frequencies.set_grid_index(j -1, ell+1);
                    *extremal_frequencies_changed = true;
                    continue 'find_extremal;
                }
                if !*extremal_frequencies_changed {
                    continue 'loop_03;
                }
                *klow = extremal_frequencies.get_grid_index(j -1);
                continue 'find_extremal;
            }
        }
    }
}

fn calculate_err(
    grid: &DenseGrid,
    zeroth_value_override: Option<f32>,
    x: &[f64; 66],
    y: &[f64; 66],
    ad: &[f64; 66],
    ell: i64,
    last_coefficient_index: usize,
) -> f32 {
    let err = grid.gee(zeroth_value_override, x, y, ad, ell, last_coefficient_index) as f32;
    (err - grid.get_des((ell - 1) as usize)) * grid.get_wt((ell - 1) as usize)
}