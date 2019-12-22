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

    comp: &mut f64,
    y1: &mut f64,
    klow: &mut i64,
    nut: &mut i32,
    jchnge: &mut i32,
    extremal_frequencies: &mut ExtremalFrequencies,
) -> EndpointSearchResult{
    let last_coefficient_index = num_coefficients + 1;
    let coefficient_off_end_index = num_coefficients + 2;

    let mut nut1 = 0;
    let mut k1 = k1;
    let mut knz = knz;

    'search_endpoints: loop {
        let ynz = *comp;
        let zeroth_grid_index = extremal_frequencies.get_grid_index(0);
        if k1 > zeroth_grid_index {
            k1 = zeroth_grid_index;
        }
        let almost_last_grid_index = extremal_frequencies.get_grid_index(last_coefficient_index - 1);
        if knz < almost_last_grid_index {
            knz = almost_last_grid_index;
        }
        nut1 = *nut;
        *nut = -nu;
        let mut ell = 0;
        let kup = k1;
        *comp = ynz * 1.00001;
        'loop_06: loop {
            ell += 1;
            if ell >= kup {
                ell = grid.n_grid() as i64 + 1;
                *klow = knz;
                *nut = -nut1;
                *comp = *y1 * 1.00001;

                ell = ell - 1;
                if ell <= *klow {
                    if *jchnge > 0 {
                        return EndpointSearchResult::KeepIteratingRemez;
                    }
                    return EndpointSearchResult::StopIteratingRemez;
                }
                let mut err = calculate_err(grid, None, x, y, ad, ell, last_coefficient_index);

                'loop_07: while (*nut as f64) * (err as f64) - *comp <= 0.0 {
                    ell = ell - 1;
                    if ell <= *klow {
                        if *jchnge > 0 {
                            return EndpointSearchResult::KeepIteratingRemez;
                        }
                        return EndpointSearchResult::StopIteratingRemez;
                    }
                    err = calculate_err(grid, None, x, y, ad, ell, last_coefficient_index);
                }

                *comp = (*nut as f64) * (err as f64);

                loop {
                    ell -= 1;
                    if ell <= *klow {
                        break;
                    }
                    err = calculate_err(grid, None, x, y, ad, ell, last_coefficient_index);
                    let dtemp = (*nut as f64) * (err as f64) - *comp;
                    if dtemp <= 0.0 {
                        break;
                    }
                    *comp = (*nut as f64) * (err as f64);
                }

                *klow = extremal_frequencies.get_grid_index(coefficient_off_end_index - 1);
                extremal_frequencies.set_grid_index(coefficient_off_end_index - 1, ell + 1);
                *jchnge += 1;
                extremal_frequencies.shift_grid_indexes_left();
                return EndpointSearchResult::KeepIteratingRemez;
            }
            let mut err = calculate_err(grid, None, x, y, ad, ell, last_coefficient_index);
            let dtemp = (*nut as f64) * (err as f64) - *comp;
            if dtemp <= 0.0 {
                continue 'loop_06;
            }

            *comp = (*nut as f64) * (err as f64);

            loop {
                ell = ell + 1;
                if ell >= kup {
                    break;
                }
                err = calculate_err(grid, None, x, y, ad, ell, last_coefficient_index);
                let dtemp = (*nut as f64) * (err as f64) - *comp;
                if dtemp <= 0.0 {
                    break;
                }
                *comp = (*nut as f64) * (err as f64);
            }

            extremal_frequencies.set_grid_index(coefficient_off_end_index - 1, ell - 1);
            *klow = ell - 1;
            *jchnge += 1;
            break 'search_endpoints;
        }
    }

    if *comp > *y1 {
        *y1 = *comp;
    }
    k1 = extremal_frequencies.get_grid_index(last_coefficient_index);
    let mut local_ell = grid.n_grid() as i64 + 1;
    *klow = knz;
    *nut = -nut1;
    *comp = *y1 * 1.00001;
    'loop_11: loop {
        local_ell -= 1;
        if local_ell <= *klow {
            extremal_frequencies.shift_grid_indexes_right(k1);
            return EndpointSearchResult::KeepIteratingRemez;
        }
        let err = calculate_err(grid, None, x, y, ad, local_ell, last_coefficient_index);
        let dtemp = (*nut as f64) * (err as f64) - *comp;
        if dtemp <= 0.0 {
            continue 'loop_11;
        }
        break 'loop_11;
    }

    loop {
        let err = calculate_err(grid, None, x, y, ad, local_ell, last_coefficient_index);
        *comp = (*nut as f64) * (err as f64);
        local_ell -= 1;
        if local_ell <= *klow {
            break;
        }
        if (*nut as f64) * (err as f64) - *comp <= 0.0 {
            break;
        }
    }
    //*klow = extremal_frequencies.get_grid_index(coefficient_off_end_index-1);
    extremal_frequencies.set_grid_index(coefficient_off_end_index-1, local_ell+1);
    //*jchnge += 1;
    extremal_frequencies.shift_grid_indexes_left();
    return EndpointSearchResult::KeepIteratingRemez;
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