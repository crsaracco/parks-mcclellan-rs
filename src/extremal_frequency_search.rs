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
    jchnge: &mut i32,
    extremal_frequencies: &mut ExtremalFrequencies,
) {
    let last_coefficient_index = num_coefficients + 1;
    let coefficient_off_end_index = num_coefficients + 2;

    let mut non_loop_j = 1;

    'find_extremal: while non_loop_j < coefficient_off_end_index {
        let kup = extremal_frequencies.get_grid_index(non_loop_j);
        let mut ell = extremal_frequencies.get_grid_index(non_loop_j-1) + 1;
        *nut = -*nut;
        if non_loop_j == 2 {
            *y1 = *comp;
        }
        *comp = deviation;
        if ell >= kup {
            ell = ell - 1;
            'loop_03: loop {
                ell = ell - 1;
                if ell <= *klow {
                    ell = extremal_frequencies.get_grid_index(non_loop_j-1) + 1;
                    if *jchnge > 0 {
                        extremal_frequencies.set_grid_index(non_loop_j-1, ell-1);
                        non_loop_j += 1;
                        *klow = ell - 1;
                        *jchnge += 1;
                        continue 'find_extremal;
                    }
                    'loop_05: loop {
                        ell += 1;
                        if ell >= kup {
                            *klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                            non_loop_j += 1;
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

                        extremal_frequencies.set_grid_index(non_loop_j-1, ell-1);
                        non_loop_j += 1;
                        *klow = ell - 1;
                        *jchnge += 1;
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

                    *klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                    extremal_frequencies.set_grid_index(non_loop_j-1, ell+1);
                    non_loop_j += 1;
                    *jchnge += 1;
                    continue 'find_extremal;
                }
                if *jchnge <= 0 {
                    continue 'loop_03;
                }
                *klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                non_loop_j += 1;
                continue 'find_extremal;
            }
        }
        let mut err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
        let dtemp = (*nut as f64) * (err as f64) - *comp;
        if dtemp <= 0.0 {
            ell = ell - 1;
            'loop_13: loop {
                ell = ell - 1;
                if ell <= *klow {
                    ell = extremal_frequencies.get_grid_index(non_loop_j-1) + 1;
                    if *jchnge > 0 {
                        extremal_frequencies.set_grid_index(non_loop_j-1, ell-1);
                        non_loop_j += 1;
                        *klow = ell - 1;
                        *jchnge += 1;
                        continue 'find_extremal;
                    }
                    'loop_15: loop {
                        ell += 1;
                        if ell >= kup {
                            *klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                            non_loop_j += 1;
                            continue 'find_extremal;
                        }
                        err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
                        let dtemp = (*nut as f64) * (err as f64) - *comp;
                        if dtemp <= 0.0 {
                            continue 'loop_15;
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

                        extremal_frequencies.set_grid_index(non_loop_j-1, ell-1);
                        non_loop_j += 1;
                        *klow = ell - 1;
                        *jchnge += 1;
                        continue 'find_extremal;
                    }
                }
                err = calculate_err(grid, None, &x, &y, &ad, ell, last_coefficient_index);
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

                    *klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                    extremal_frequencies.set_grid_index(non_loop_j-1, ell+1);
                    non_loop_j += 1;
                    *jchnge += 1;
                    continue 'find_extremal;
                }
                if *jchnge <= 0 {
                    continue 'loop_13;
                }
                *klow = extremal_frequencies.get_grid_index(non_loop_j-1);
                non_loop_j += 1;
                continue 'find_extremal;
            }
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

        extremal_frequencies.set_grid_index(non_loop_j-1, ell-1);
        non_loop_j += 1;
        *klow = ell - 1;
        *jchnge += 1;
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