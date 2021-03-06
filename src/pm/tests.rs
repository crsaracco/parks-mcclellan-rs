use super::*;
use crate::Band;

const FLOATING_POINT_TOLERANCE: f64 = std::f32::EPSILON as f64;

fn test_f64_error_within_tolerance(measured: f64, expected: f64) {
    // Don't try to divide by zero.
    if expected == 0.0 {
        assert_eq!(measured.abs(), 0.0);
        return;
    }

    let error = (measured - expected) / expected;
    assert!(
        error.abs() <= FLOATING_POINT_TOLERANCE,
        format!("measured: {}, expected: {}, error: {}", measured, expected, error)
    );
}

fn test_inputs(filter_output: &ParksMcclellanOutput, expected_bands: &Vec<Band>) {
    assert_eq!(filter_output.lower_band_edges.len(), expected_bands.len());
    assert_eq!(filter_output.upper_band_edges.len(), expected_bands.len());
    assert_eq!(filter_output.desired_values.len(), expected_bands.len());
    assert_eq!(filter_output.weightings.len(), expected_bands.len());
    for (b, band) in expected_bands.iter().enumerate() {
        assert_eq!(filter_output.lower_band_edges[b], band.lower_edge);
        assert_eq!(filter_output.upper_band_edges[b], band.upper_edge);
        assert_eq!(filter_output.desired_values[b], band.desired_value);
        assert_eq!(filter_output.weightings[b], band.weight);
    }
}

fn test_impulse_response(filter_output: &ParksMcclellanOutput, expected: &Vec<f64>) {
    println!("Impulse response: {:?}", filter_output.impulse_response);
    assert_eq!(filter_output.impulse_response.len(), expected.len());
    for h in 0..filter_output.impulse_response.len() {
        test_f64_error_within_tolerance(
            filter_output.impulse_response[h] as f64,
            expected[h],
        );
    }
}

fn test_deviations(filter_output: &ParksMcclellanOutput, expected: &Vec<f64>) {
    assert_eq!(filter_output.deviations.len(), expected.len());
    for dev in 0..filter_output.deviations.len() {
        test_f64_error_within_tolerance(
            filter_output.deviations[dev] as f64,
            expected[dev],
        );
    }
}

fn test_deviation_dbs(filter_output: &ParksMcclellanOutput, expected: &Vec<f64>) {
    assert_eq!(filter_output.deviation_dbs.len(), expected.len());
    for db in 0..filter_output.deviation_dbs.len() {
        test_f64_error_within_tolerance(
            filter_output.deviation_dbs[db] as f64,
            expected[db],
        );
    }
}

fn test_extremal_frequencies(filter_output: &ParksMcclellanOutput, expected: &Vec<f64>) {
    assert_eq!(filter_output.extremal_frequencies.len(), expected.len());
    for f in 0..filter_output.extremal_frequencies.len() {
        test_f64_error_within_tolerance(
            filter_output.extremal_frequencies[f] as f64,
            expected[f],
        );
    }
}

/////////////////////
// EXAMPLE FILTERS //
/////////////////////

// This is the first example in the documentation comments -- a multi-band filter.
/*
32,1,3,0
0.0,0.1,0.2,0.35
0.425,0.5
0.0,1.0,0.0
10.0,1.0,10.0
*/
#[test]
fn example_multiband_filter() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0.0,
        upper_edge: 0.1,
        desired_value: 0.0,
        weight: 10.0,
    });
    bands.push(Band {
        lower_edge: 0.2,
        upper_edge: 0.35,
        desired_value: 1.0,
        weight: 1.0,
    });
    bands.push(Band {
        lower_edge: 0.425,
        upper_edge: 0.5,
        desired_value: 0.0,
        weight: 10.0,
    });

    let pm_output = design_filter(32, FilterType::MultipleBand, &bands, 16);

    assert_eq!(pm_output.filter_length, 32);

    let expected_impulse_response = vec![
        -0.5753406789E-02,  0.9902641177E-03,  0.7573357783E-02, -0.6514116656E-02,
         0.1396050397E-01,  0.2295166254E-02, -0.1999405399E-01,  0.7136960514E-02,
        -0.3965735435E-01,  0.1126007363E-01,  0.6623364240E-01, -0.1049720775E-01,
         0.8513614535E-01, -0.1202498823E+00, -0.2967857718E+00,  0.3041091561E+00,
    ];
    let expected_deviations = vec![
        0.1513117226E-02, 0.1513117179E-01, 0.1513117226E-02
    ];
    let expected_deviation_dbs = vec![
        -0.5640254974E+02, 0.1304427981E+00, -0.5640254974E+02
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000E+00, 0.2734375000E-01, 0.5273437500E-01, 0.7617187500E-01, 0.9375000000E-01,
        0.1000000015E+00, 0.2000000030E+00, 0.2195312530E+00, 0.2527343631E+00, 0.2839843631E+00,
        0.3132812381E+00, 0.3386718631E+00, 0.3499999940E+00, 0.4250000119E+00, 0.4328125119E+00,
        0.4503906369E+00, 0.4796875119E+00,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

// This is the second example in the documentation comments -- a differentiator filter.
/*
32,2,1,20
0,0.5
1.0
1.0
*/
#[test]
fn example_differentiator() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0.0,
        upper_edge: 0.5,
        desired_value: 1.0,
        weight: 1.0,
    });

    let pm_output = design_filter(32, FilterType::Differentiator, &bands, 20);

    assert_eq!(pm_output.filter_length, 32);

    let expected_impulse_response = vec![
        -0.6257381756E-03,  0.8548910846E-03, -0.4238411202E-03,  0.3989883116E-03,
        -0.4351424286E-03,  0.5008199369E-03, -0.5971363280E-03,  0.7333782269E-03,
        -0.9301801911E-03,  0.1226891065E-02, -0.1701220521E-02,  0.2527206670E-02,
        -0.4160370212E-02,  0.8130028844E-02, -0.2253980562E-01,  0.2026662678E+00,
    ];
    let expected_deviations = vec![0.6192352157E-02];
    let expected_extremal_frequencies = vec![
        0.1562500023E-02, 0.3281250224E-01, 0.6562497467E-01, 0.9843744338E-01, 0.1312499344E+00,
        0.1656250656E+00, 0.1984376907E+00, 0.2312503159E+00, 0.2640629411E+00, 0.2968755662E+00,
        0.3296881914E+00, 0.3625008166E+00, 0.3953134418E+00, 0.4281260669E+00, 0.4593761861E+00,
        0.4859387875E+00, 0.5000000000E+00,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/////////////////////////
// GENERATED BY SCRIPT //
/////////////////////////

/*
3,1,2,3
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n3_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::MultipleBand, &bands, 3);

    let expected_impulse_response = vec![
        3.8196599480e-1,    5.0000000000e-1,
    ];
    let expected_deviations = vec![
        2.6393201950e-1,    2.6393201950e-1,
    ];
    let expected_deviation_dbs = vec![
        2.0344741340e0,    -1.1570158960e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,1,2,4
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n3_grid4() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::MultipleBand, &bands, 4);

    let expected_impulse_response = vec![
        3.8196599480e-1,    5.0000000000e-1,
    ];
    let expected_deviations = vec![
        2.6393201950e-1,    2.6393201950e-1,
    ];
    let expected_deviation_dbs = vec![
        2.0344741340e0,    -1.1570158960e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,1,2,16
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n3_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        3.8196599480e-1,    5.0000000000e-1,
    ];
    let expected_deviations = vec![
        2.6393201950e-1,    2.6393201950e-1,
    ];
    let expected_deviation_dbs = vec![
        2.0344741340e0,    -1.1570158960e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,1,2,17
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n3_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        3.8196599480e-1,    5.0000000000e-1,
    ];
    let expected_deviations = vec![
        2.6393201950e-1,    2.6393201950e-1,
    ];
    let expected_deviation_dbs = vec![
        2.0344741340e0,    -1.1570158960e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,1,2,61
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n3_grid61() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::MultipleBand, &bands, 61);

    let expected_impulse_response = vec![
        3.8196599480e-1,    5.0000000000e-1,
    ];
    let expected_deviations = vec![
        2.6393201950e-1,    2.6393201950e-1,
    ];
    let expected_deviation_dbs = vec![
        2.0344741340e0,    -1.1570158960e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,1,2,64
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n3_grid64() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::MultipleBand, &bands, 64);

    let expected_impulse_response = vec![
        3.8196599480e-1,    5.0000000000e-1,
    ];
    let expected_deviations = vec![
        2.6393201950e-1,    2.6393201950e-1,
    ];
    let expected_deviation_dbs = vec![
        2.0344741340e0,    -1.1570158960e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,1,2,3
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n4_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::MultipleBand, &bands, 3);

    let expected_impulse_response = vec![
        1.5172423420e-1,    4.9483287330e-1,
    ];
    let expected_deviations = vec![
        2.9311430450e-1,    2.9311430450e-1,
    ];
    let expected_deviation_dbs = vec![
        2.2327382560e0,    -1.0659259800e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,1,2,4
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n4_grid4() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::MultipleBand, &bands, 4);

    let expected_impulse_response = vec![
        1.5172423420e-1,    4.9483287330e-1,
    ];
    let expected_deviations = vec![
        2.9311430450e-1,    2.9311430450e-1,
    ];
    let expected_deviation_dbs = vec![
        2.2327382560e0,    -1.0659259800e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,1,2,16
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n4_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        1.5172423420e-1,    4.9483287330e-1,
    ];
    let expected_deviations = vec![
        2.9311430450e-1,    2.9311430450e-1,
    ];
    let expected_deviation_dbs = vec![
        2.2327382560e0,    -1.0659259800e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,1,2,17
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n4_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        1.5172423420e-1,    4.9483287330e-1,
    ];
    let expected_deviations = vec![
        2.9311430450e-1,    2.9311430450e-1,
    ];
    let expected_deviation_dbs = vec![
        2.2327382560e0,    -1.0659259800e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,1,2,61
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n4_grid61() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::MultipleBand, &bands, 61);

    let expected_impulse_response = vec![
        1.5172423420e-1,    4.9483287330e-1,
    ];
    let expected_deviations = vec![
        2.9311430450e-1,    2.9311430450e-1,
    ];
    let expected_deviation_dbs = vec![
        2.2327382560e0,    -1.0659259800e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,1,2,64
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n4_grid64() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::MultipleBand, &bands, 64);

    let expected_impulse_response = vec![
        1.5172423420e-1,    4.9483287330e-1,
    ];
    let expected_deviations = vec![
        2.9311430450e-1,    2.9311430450e-1,
    ];
    let expected_deviation_dbs = vec![
        2.2327382560e0,    -1.0659259800e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.0000000300e-1,    3.0000001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,3
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n16_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 3);

    let expected_impulse_response = vec![
        -1.1928418650e-2,   -2.2930486130e-2,    2.9667090620e-2,    3.5365093500e-2,
        -5.4441429670e-2,   -8.2288876180e-2,    1.4648470280e-1,    4.4791793820e-1,
    ];
    let expected_deviations = vec![
        2.4308655410e-2,    2.4308655410e-2,
    ];
    let expected_deviation_dbs = vec![
        2.0861707630e-1,    -3.2284782410e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.2500000000e-2,    1.2500000000e-1,    1.6666665670e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.2083335520e-1,    3.8333338500e-1,
        4.6666675810e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,4
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n16_grid4() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 4);

    let expected_impulse_response = vec![
        -1.1470410040e-2,   -2.3657456040e-2,    2.8824593870e-2,    3.5916417840e-2,
        -5.4141830650e-2,   -8.2091890280e-2,    1.4592386780e-1,    4.4812351470e-1,
    ];
    let expected_deviations = vec![
        2.5146465750e-2,    2.5146465750e-2,
    ];
    let expected_deviation_dbs = vec![
        2.1571853760e-1,    -3.1990461350e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.2500000000e-2,    1.2500000000e-1,    1.7187500000e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.3125001190e-1,    3.7812501190e-1,
        4.5625001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,16
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n16_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -1.1686353010e-2,   -2.3650348190e-2,    2.8849342840e-2,    3.5625636580e-2,
        -5.3930357100e-2,   -8.2009784880e-2,    1.4596554640e-1,    4.4800698760e-1,
    ];
    let expected_deviations = vec![
        2.5658687580e-2,    2.5658687580e-2,
    ];
    let expected_deviation_dbs = vec![
        2.2005759180e-1,    -3.1815311430e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.2500000000e-2,    1.2500000000e-1,    1.7578125000e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.2734376190e-1,    3.8203126190e-1,
        4.5625001190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,17
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n16_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        -1.1721549560e-2,   -2.3628782480e-2,    2.8876205910e-2,    3.5589434210e-2,
        -5.3939655420e-2,   -8.1979833540e-2,    1.4596402650e-1,    4.4799768920e-1,
    ];
    let expected_deviations = vec![
        2.5685003030e-2,    2.5685003030e-2,
    ];
    let expected_deviation_dbs = vec![
        2.2027969360e-1,    -3.1806407930e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.2500000000e-2,    1.2500005960e-1,    1.7647069690e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.2573533060e-1,    3.8455891610e-1,
        4.5808839800e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,61
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n16_grid61() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 61);

    let expected_impulse_response = vec![
        -1.1719172820e-2,   -2.3625576870e-2,    2.8871141370e-2,    3.5584531720e-2,
        -5.3938426080e-2,   -8.1987418230e-2,    1.4596895870e-1,    4.4799792770e-1,
    ];
    let expected_deviations = vec![
        2.5696048510e-2,    2.5696048510e-2,
    ];
    let expected_deviation_dbs = vec![
        2.2037357090e-1,    -3.1802673340e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.3524559140e-2,    1.2397530670e-1,    1.7622934280e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.2561510800e-1,    3.8299292330e-1,
        4.5778900380e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,64
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n16_grid64() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 64);

    let expected_impulse_response = vec![
        -1.1719881560e-2,   -2.3624714460e-2,    2.8872091320e-2,    3.5583749410e-2,
        -5.3938880560e-2,   -8.1987082960e-2,    1.4596913750e-1,    4.4799762960e-1,
    ];
    let expected_deviations = vec![
        2.5695905090e-2,    2.5695905090e-2,
    ];
    let expected_deviation_dbs = vec![
        2.2037255760e-1,    -3.1802721020e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.3476562500e-2,    1.2402343750e-1,    1.7578125000e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.2539063690e-1,    3.8300782440e-1,
        4.5820313690e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
33,1,2,3
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n33_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(33, FilterType::MultipleBand, &bands, 3);

    let expected_impulse_response = vec![
        -9.7805095720e-5,   -1.9461156330e-3,    2.7027918260e-4,    4.5151663940e-3,
        -4.4772509140e-4,   -9.2990165580e-3,    6.6755095030e-4,    1.7102934420e-2,
        -8.7346771030e-4,   -2.9697250570e-2,    1.0580382080e-3,    5.1456972960e-2,
        -1.2141985350e-3,   -9.8371602590e-2,    1.3149621660e-3,    3.1566381450e-1,
        4.9864473940e-1,
    ];
    let expected_deviations = vec![
        1.1501644040e-3,    1.1501644040e-3,
    ];
    let expected_deviation_dbs = vec![
        9.9841728810e-3,    -5.8784797670e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.9411766680e-2,    5.8823529630e-2,    8.8235296310e-2,
        1.1764705930e-1,    1.4705882970e-1,    1.7647059260e-1,    1.8627451360e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0980393290e-1,    3.2941177490e-1,
        3.4901961680e-1,    3.7843137980e-1,    4.0784314270e-1,    4.3725490570e-1,
        4.6666666870e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
33,1,2,4
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n33_grid4() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(33, FilterType::MultipleBand, &bands, 4);

    let expected_impulse_response = vec![
        -1.4601272600e-5,   -2.1017156540e-3,    1.8827044190e-5,    4.6686795540e-3,
        -2.5761864890e-5,   -9.4365095720e-3,    5.3748757640e-5,    1.7230607570e-2,
        -7.4639872760e-5,   -2.9805587600e-2,    9.0192217610e-5,    5.1532682030e-2,
        -9.8357420940e-5,   -9.8417893050e-2,    1.0190706230e-4,    3.1567656990e-1,
        4.9989736080e-1,
    ];
    let expected_deviations = vec![
        1.3063276420e-3,    1.3063276420e-3,
    ];
    let expected_deviation_dbs = vec![
        1.1338932440e-2,    -5.7678955080e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.9411764820e-2,    5.8823529630e-2,    8.8235296310e-2,
        1.2500000000e-1,    1.4705884460e-1,    1.7647063730e-1,    1.9117653370e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0735296010e-1,    3.2941180470e-1,
        3.5147064920e-1,    3.8088244200e-1,    4.1029423480e-1,    4.3970602750e-1,
        4.6911782030e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
33,1,2,16
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n33_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(33, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -2.0668012720e-7,   -2.1073303650e-3,    4.6326920260e-7,    4.6487320210e-3,
        -1.5791902110e-6,   -9.4253476710e-3,    2.1741527690e-6,    1.7206747080e-2,
        -2.1118710270e-6,   -2.9781030490e-2,    2.8965414460e-6,    5.1522340630e-2,
        -3.1998238230e-6,   -9.8418161270e-2,    3.7095051080e-6,    3.1567791100e-1,
        4.9999570850e-1,
    ];
    let expected_deviations = vec![
        1.3522884110e-3,    1.3522884110e-3,
    ];
    let expected_deviation_dbs = vec![
        1.1738082390e-2,    -5.7378612520e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    3.1250000000e-2,    6.2500029800e-2,    9.1911822560e-2,
        1.2132361530e-1,    1.4889717100e-1,    1.7463248970e-1,    1.9301486020e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0735296010e-1,    3.2573533060e-1,
        3.5147064920e-1,    3.7904420500e-1,    4.0845599770e-1,    4.3786779050e-1,
        4.6911782030e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
33,1,2,17
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n33_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(33, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        -1.5266058430e-6,   -2.1051694640e-3,    4.7413591350e-6,    4.6467911450e-3,
        -8.0761365100e-6,   -9.4230296090e-3,    1.1134432500e-5,    1.7204781990e-2,
        -1.4187123270e-5,   -2.9780218380e-2,    1.8084383560e-5,    5.1521223040e-2,
        -2.0547233360e-5,   -9.8417185250e-2,    2.1976667990e-5,    3.1567797060e-1,
        4.9997681380e-1,
    ];
    let expected_deviations = vec![
        1.3496981700e-3,    1.3496981700e-3,
    ];
    let expected_deviation_dbs = vec![
        1.1715333910e-2,    -5.7395267490e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    3.1141862270e-2,    6.2283717100e-2,    9.1695532200e-2,
        1.2110734730e-1,    1.4878895880e-1,    1.7474044860e-1,    1.9204144180e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0692046880e-1,    3.2595172520e-1,
        3.5017332430e-1,    3.7785515190e-1,    4.0726709370e-1,    4.3840914960e-1,
        4.6955120560e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
33,1,2,61
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n33_grid61() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(33, FilterType::MultipleBand, &bands, 61);

    let expected_impulse_response = vec![
        -2.8809511220e-8,   -2.1071133670e-3,    8.6646210210e-8,    4.6479739250e-3,
        -1.2441832100e-7,   -9.4247106460e-3,    2.3915569610e-7,    1.7205873500e-2,
        -3.3145545330e-7,   -2.9780976470e-2,    4.0247937250e-7,    5.1521617920e-2,
        -4.6643160090e-7,   -9.8417378960e-2,    4.5881864710e-7,    3.1567794080e-1,
        4.9999952320e-1,
    ];
    let expected_deviations = vec![
        1.3535725880e-3,    1.3535725880e-3,
    ];
    let expected_deviation_dbs = vec![
        1.1749455710e-2,    -5.7370368960e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    3.0858255920e-2,    6.1716534200e-2,    9.2092424630e-2,
        1.2150399390e-1,    1.4946909250e-1,    1.7405909300e-1,    1.9286321100e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0723258850e-1,    3.2603728770e-1,
        3.5062804820e-1,    3.7859401110e-1,    4.0800648930e-1,    4.3838331100e-1,
        4.6924230460e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
33,1,2,64
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n33_grid64() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(33, FilterType::MultipleBand, &bands, 64);

    let expected_impulse_response = vec![
        -9.3449735060e-8,   -2.1069769280e-3,    2.4235069420e-7,    4.6478137370e-3,
        -4.3768304180e-7,   -9.4245206560e-3,    6.6265664600e-7,    1.7205690970e-2,
        -8.8533198550e-7,   -2.9780842360e-2,    1.0756434680e-6,    5.1521547140e-2,
        -1.2065131610e-6,   -9.8417341710e-2,    1.3167365300e-6,    3.1567791100e-1,
        4.9999865890e-1,
    ];
    let expected_deviations = vec![
        1.3534586180e-3,    1.3534586180e-3,
    ];
    let expected_deviation_dbs = vec![
        1.1748421940e-2,    -5.7371101380e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    3.0790463090e-2,    6.1580933630e-2,    9.1911844910e-2,
        1.2132363770e-1,    1.4935635030e-1,    1.7417214810e-1,    1.9255422060e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0735284090e-1,    3.2573491330e-1,
        3.5055071120e-1,    3.7858337160e-1,    4.0799468760e-1,    4.3832510710e-1,
        4.6911507840e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
34,1,2,3
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n34_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(34, FilterType::MultipleBand, &bands, 3);

    let expected_impulse_response = vec![
        5.7205540360e-4,   -1.3089210260e-3,   -1.9107917320e-3,    2.9390624260e-3,
        4.1029313580e-3,   -6.0415607880e-3,   -7.9041766000e-3,    1.1051053180e-2,
        1.4058675620e-2,   -1.9002478570e-2,   -2.4063944820e-2,    3.2273974270e-2,
        4.2025640610e-2,   -5.8625299480e-2,   -8.5191436110e-2,    1.4777301250e-1,
        4.4883292910e-1,
    ];
    let expected_deviations = vec![
        8.3852757230e-4,    8.3852757230e-4,
    ];
    let expected_deviation_dbs = vec![
        7.2802240030e-3,    -6.1529655460e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.9411766680e-2,    5.8823529630e-2,    9.8039217290e-2,
        1.2745098770e-1,    1.4705882970e-1,    1.7647059260e-1,    1.8627451360e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0980393290e-1,    3.1960785390e-1,
        3.4901961680e-1,    3.6862745880e-1,    3.9803922180e-1,    4.2745098470e-1,
        4.5686274770e-1,    4.8627451060e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
34,1,2,4
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n34_grid4() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(34, FilterType::MultipleBand, &bands, 4);

    let expected_impulse_response = vec![
        6.8972900040e-4,   -1.3753577370e-3,   -2.1243463270e-3,    2.9001233630e-3,
        4.3872031380e-3,   -5.9462627400e-3,   -8.3188908170e-3,    1.0838104410e-2,
        1.4540242030e-2,   -1.8672540780e-2,   -2.4618841710e-2,    3.1838506460e-2,
        4.2643304910e-2,   -5.8111682530e-2,   -8.5827760400e-2,    1.4720779660e-1,
        4.4945353270e-1,
    ];
    let expected_deviations = vec![
        9.9426682570e-4,    9.9426682570e-4,
    ];
    let expected_deviation_dbs = vec![
        8.6323032160e-3,    -6.0049942020e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.9411764820e-2,    6.6176474090e-2,    9.5588237050e-2,
        1.2500000000e-1,    1.5441179280e-1,    1.7647063730e-1,    1.9117653370e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0735296010e-1,    3.2205885650e-1,
        3.4411770110e-1,    3.7352949380e-1,    4.0294128660e-1,    4.2500013110e-1,
        4.5441192390e-1,    4.8382371660e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
34,1,2,16
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n34_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(34, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        6.9417728810e-4,   -1.3655875810e-3,   -2.1471329960e-3,    2.8724223380e-3,
        4.4264849280e-3,   -5.9042852370e-3,   -8.3666956050e-3,    1.0774061080e-2,
        1.4614485200e-2,   -1.8589740620e-2,   -2.4714220320e-2,    3.1736750160e-2,
        4.2743787170e-2,   -5.7998478410e-2,   -8.5939452050e-2,    1.4708261190e-1,
        4.4956949350e-1,
    ];
    let expected_deviations = vec![
        1.0226215240e-3,    1.0226215240e-3,
    ];
    let expected_deviation_dbs = vec![
        8.8774552570e-3,    -5.9805702210e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    3.3088237050e-2,    6.4338266850e-2,    9.5588296650e-2,
        1.2500008940e-1,    1.5257364510e-1,    1.7647072670e-1,    1.9301486020e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0735296010e-1,    3.2389709350e-1,
        3.4595593810e-1,    3.7169125680e-1,    3.9926481250e-1,    4.2867660520e-1,
        4.5625016090e-1,    4.8566195370e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
34,1,2,17
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n34_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(34, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        6.9215276740e-4,   -1.3687446480e-3,   -2.1421585700e-3,    2.8776898980e-3,
        4.4190855700e-3,   -5.9134559710e-3,   -8.3549153060e-3,    1.0787625800e-2,
        1.4598725360e-2,   -1.8606901170e-2,   -2.4694392460e-2,    3.1757265330e-2,
        4.2720634490e-2,   -5.8021150530e-2,   -8.5914939640e-2,    1.4710736270e-1,
        4.4954389330e-1,
    ];
    let expected_deviations = vec![
        1.0244029110e-3,    1.0244029110e-3,
    ];
    let expected_deviation_dbs = vec![
        8.8929710910e-3,    -5.9790580750e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    3.2871965320e-2,    6.4013823870e-2,    9.5155745740e-2,
        1.2456756090e-1,    1.5224915740e-1,    1.7647054790e-1,    1.9377154110e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0692046880e-1,    3.2422161100e-1,
        3.4671309590e-1,    3.7266480920e-1,    4.0034663680e-1,    4.2802846430e-1,
        4.5744040610e-1,    4.8512223360e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
34,1,2,61
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n34_grid61() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(34, FilterType::MultipleBand, &bands, 61);

    let expected_impulse_response = vec![
        6.9337076270e-4,   -1.3688632750e-3,   -2.1450170320e-3,    2.8764603190e-3,
        4.4223614970e-3,   -5.9110261500e-3,   -8.3593660970e-3,    1.0783646260e-2,
        1.4604222030e-2,   -1.8601924180e-2,   -2.4701211600e-2,    3.1750779600e-2,
        4.2728684840e-2,   -5.8014158160e-2,   -8.5923075680e-2,    1.4709931610e-1,
        4.4955244660e-1,
    ];
    let expected_deviations = vec![
        1.0267206230e-3,    1.0267206230e-3,
    ];
    let expected_deviation_dbs = vec![
        8.9136585590e-3,    -5.9770957950e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    3.2786898310e-2,    6.4609482880e-2,    9.5467522740e-2,
        1.2487909200e-1,    1.5187987690e-1,    1.7598772050e-1,    1.9334536790e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0675041680e-1,    3.2362642880e-1,
        3.4677067400e-1,    3.7232577800e-1,    3.9980956910e-1,    4.2825770380e-1,
        4.5670583840e-1,    4.8563614490e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
34,1,2,64
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n34_grid64() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(34, FilterType::MultipleBand, &bands, 64);

    let expected_impulse_response = vec![
        6.9321598860e-4,   -1.3688101900e-3,   -2.1447176110e-3,    2.8765536840e-3,
        4.4219531120e-3,   -5.9112221930e-3,   -8.3587421100e-3,    1.0784052310e-2,
        1.4603463930e-2,   -1.8602544440e-2,   -2.4700300770e-2,    3.1751561910e-2,
        4.2727679010e-2,   -5.8015033600e-2,   -8.5922054950e-2,    1.4710026980e-1,
        4.4955140350e-1,
    ];
    let expected_deviations = vec![
        1.0265249290e-3,    1.0265249290e-3,
    ];
    let expected_deviation_dbs = vec![
        8.9115891610e-3,    -5.9772609710e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    3.2628700140e-2,    6.4797848460e-2,    9.5588319000e-2,
        1.2454055250e-1,    1.5211366120e-1,    1.7555080350e-1,    1.9347332420e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0643373730e-1,    3.2343715430e-1,
        3.4641474490e-1,    3.7260919810e-1,    3.9972275500e-1,    4.2821496730e-1,
        4.5670717950e-1,    4.8565894370e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
103,1,2,3
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n103_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(103, FilterType::MultipleBand, &bands, 3);

    let expected_impulse_response = vec![
        -3.7950332650e-8,    2.3169121730e-8,    2.3372068370e-7,   -9.3711356900e-8,
        -9.1653976140e-7,    3.0364284950e-7,    2.8392903460e-6,   -8.1638785330e-7,
        -7.5155890040e-6,    1.9159524530e-6,    1.7724798450e-5,   -4.0614891080e-6,
        -3.8220117860e-5,    7.9305063990e-6,    7.6629767140e-5,   -1.4459224990e-5,
        -1.4455559720e-4,    2.4865878000e-5,    2.5884242500e-4,   -4.0620176150e-5,
        -4.4293093380e-4,    6.3388208220e-5,    7.2828284460e-4,   -9.4904440630e-5,
        -1.1558139230e-3,    1.3679912080e-4,    1.7774141160e-3,   -1.9038065510e-4,
        -2.6577806570e-3,    2.5639295930e-4,    3.8771107790e-3,   -3.3478200200e-4,
        -5.5357995440e-3,    4.2450116600e-4,    7.7634216290e-3,   -5.2339321700e-4,
        -1.0736690830e-2,    6.2818860170e-4,    1.4716692270e-2,   -7.3461222930e-4,
        -2.0130241290e-2,    8.3764363080e-4,    2.7763916180e-2,   -9.3187362650e-4,
        -3.9295077320e-2,    1.0119489160e-3,    5.9100992980e-2,   -1.0730552020e-3,
        -1.0330528020e-1,    1.1113656220e-3,    3.1736677890e-1,    4.9887558820e-1,
    ];
    let expected_deviations = vec![
        3.5561236180e-9,    3.5561236180e-9,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.6898046880e2,
    ];
    let expected_extremal_frequencies = vec![
        6.4102564940e-3,    1.2820512990e-2,    2.2435897960e-2,    2.8846153990e-2,
        3.8461539890e-2,    4.4871795920e-2,    5.1282051950e-2,    6.7307695750e-2,
        7.3717951770e-2,    8.6538463830e-2,    9.6153847870e-2,    1.0576923190e-1,
        1.1217948790e-1,    1.2500000000e-1,    1.3141027090e-1,    1.4423081280e-1,
        1.5064108370e-1,    1.6025649010e-1,    1.6666676100e-1,    1.7628216740e-1,
        1.7948730290e-1,    1.8589757380e-1,    1.8910270930e-1,    1.9230784480e-1,
        1.9551298020e-1,    2.0000000300e-1,    3.0000001190e-1,    3.0320513250e-1,
        3.0641025300e-1,    3.0961537360e-1,    3.1282049420e-1,    3.1923073530e-1,
        3.2564097640e-1,    3.3205121760e-1,    3.4166657920e-1,    3.4807682040e-1,
        3.5769218210e-1,    3.6410242320e-1,    3.7371778490e-1,    3.8333314660e-1,
        3.8974338770e-1,    3.9935874940e-1,    4.0897411110e-1,    4.1858947280e-1,
        4.2499971390e-1,    4.3461507560e-1,    4.4423043730e-1,    4.5384579900e-1,
        4.6346116070e-1,    4.7307652240e-1,    4.8269188400e-1,    4.9230724570e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
103,1,2,4
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n103_grid4() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(103, FilterType::MultipleBand, &bands, 4);

    let expected_impulse_response = vec![
        -3.5749888380e-8,    4.4603616800e-8,    2.4496910100e-7,   -1.8862969850e-7,
        -9.9714202410e-7,    5.7697593550e-7,    3.1143924840e-6,   -1.4724787430e-6,
        -8.2455862870e-6,    3.2940497480e-6,    1.9374483600e-5,   -6.6743486970e-6,
        -4.1519975640e-5,    1.2527090800e-5,    8.2663413200e-5,   -2.2047968740e-5,
        -1.5480071310e-4,    3.6730354620e-5,    2.7516484260e-4,   -5.8321111280e-5,
        -4.6752949130e-4,    8.8719905760e-5,    7.6354271730e-4,   -1.2982690530e-4,
        -1.2040785510e-3,    1.8334471680e-4,    1.8406949240e-3,   -2.5052434650e-4,
        -2.7374033820e-3,    3.3192845880e-4,    3.9733783340e-3,   -4.2717825270e-4,
        -5.6476849130e-3,    5.3478067280e-4,    7.8883674000e-3,   -6.5203628040e-4,
        -1.0870564730e-2,    7.7505927770e-4,    1.4853915200e-2,   -8.9893123370e-4,
        -2.0264107730e-2,    1.0179983220e-3,    2.7887109670e-2,   -1.1262479240e-3,
        -3.9400275800e-2,    1.2178077600e-3,    5.9181556110e-2,   -1.2874382080e-3,
        -1.0335590690e-1,    1.3309903440e-3,    3.1738403440e-1,    4.9865418670e-1,
    ];
    let expected_deviations = vec![
        4.5364307870e-9,    4.5364307870e-9,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.6686570740e2,
    ];
    let expected_extremal_frequencies = vec![
        7.2115389630e-3,    1.4423076990e-2,    2.1634615960e-2,    2.8846153990e-2,
        3.8461543620e-2,    4.5673087240e-2,    5.2884630860e-2,    6.0096174480e-2,
        6.7307718100e-2,    8.4134653210e-2,    8.8942348960e-2,    1.0336543620e-1,
        1.1057697980e-1,    1.2500005960e-1,    1.3221158090e-1,    1.4423078300e-1,
        1.5144230430e-1,    1.6105766590e-1,    1.7067302760e-1,    1.7788454890e-1,
        1.8509607020e-1,    1.8990375100e-1,    1.9230759140e-1,    1.9471143190e-1,
        1.9711527230e-1,    2.0000000300e-1,    3.0000001190e-1,    3.0240386720e-1,
        3.0480772260e-1,    3.0721157790e-1,    3.1201928850e-1,    3.1923085450e-1,
        3.2403856520e-1,    3.3125013110e-1,    3.4086555240e-1,    3.4807711840e-1,
        3.5528868440e-1,    3.6490410570e-1,    3.7211567160e-1,    3.8173109290e-1,
        3.9134651420e-1,    3.9855808020e-1,    4.0817350150e-1,    4.1778892280e-1,
        4.2740434410e-1,    4.3461591010e-1,    4.4423133130e-1,    4.5384675260e-1,
        4.6346217390e-1,    4.7307759520e-1,    4.8028916120e-1,    4.8990458250e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
103,1,2,16
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n103_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(103, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -4.4458033700e-8,    5.8217999310e-8,    3.1606907670e-7,   -2.2172622490e-7,
        -1.2673950780e-6,    5.9636391820e-7,    3.8367170420e-6,   -1.3833308690e-6,
        -9.8493555920e-6,    2.9119917140e-6,    2.2576394260e-5,   -5.6246594800e-6,
        -4.7400699260e-5,    1.0121000740e-5,    9.2700873210e-5,   -1.7174708770e-5,
        -1.7088801540e-4,    2.7746318660e-5,    2.9960155370e-4,   -4.2911764470e-5,
        -5.0292070960e-4,    6.3785337260e-5,    8.1259309080e-4,   -9.1435918880e-5,
        -1.2692926680e-3,    1.2680895450e-4,    1.9240583060e-3,   -1.7054704950e-4,
        -2.8400404840e-3,    2.2282543070e-4,    4.0951729750e-3,   -2.8323548030e-4,
        -5.7869576850e-3,    3.5074795600e-4,    8.0417217690e-3,   -4.2364816180e-4,
        -1.1032909150e-2,    4.9953069540e-4,    1.5018650330e-2,   -5.7539332190e-4,
        -2.0423481240e-2,    6.4787024170e-4,    2.8032759200e-2,   -7.1345688780e-4,
        -3.9523962890e-2,    7.6874409570e-4,    5.9275906530e-2,   -8.1066513670e-4,
        -1.0341505710e-1,    8.3681783870e-4,    3.1740418080e-1,    4.9915429950e-1,
    ];
    let expected_deviations = vec![
        6.4661058730e-9,    6.4661058730e-9,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.6378713990e2,
    ];
    let expected_extremal_frequencies = vec![
        7.2115384970e-3,    1.5024043620e-2,    2.2836549210e-2,    3.0649054800e-2,
        3.7860576060e-2,    4.5673057440e-2,    5.3485538810e-2,    6.1298020180e-2,
        6.9110542540e-2,    7.9326927660e-2,    8.8341385130e-2,    9.9759697910e-2,
        1.0757222770e-1,    1.1959150430e-1,    1.2259632350e-1,    1.4543294910e-1,
        1.5024065970e-1,    1.6045704480e-1,    1.6766861080e-1,    1.7668306830e-1,
        1.8209174280e-1,    1.8930330870e-1,    1.9351005550e-1,    1.9771680240e-1,
        1.9891873000e-1,    2.0000000300e-1,    3.0000001190e-1,    3.0060097580e-1,
        3.0300483110e-1,    3.0721157790e-1,    3.1201928850e-1,    3.1802892680e-1,
        3.2463952900e-1,    3.3185109500e-1,    3.3906266090e-1,    3.4687519070e-1,
        3.5528868440e-1,    3.6370217800e-1,    3.7211567160e-1,    3.8113012910e-1,
        3.9014458660e-1,    3.9915904400e-1,    4.0817350150e-1,    4.1718795900e-1,
        4.2620241640e-1,    4.3521687390e-1,    4.4423133130e-1,    4.5384675260e-1,
        4.6286121010e-1,    4.7187566760e-1,    4.8149108890e-1,    4.9050554630e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
103,1,2,17
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n103_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(103, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        -4.1608004150e-8,    6.1180983830e-8,    3.0240931890e-7,   -2.4101981920e-7,
        -1.2318109840e-6,    6.5858966990e-7,    3.7593708840e-6,   -1.5366970270e-6,
        -9.6880085040e-6,    3.2543086950e-6,    2.2267347960e-5,   -6.3307029450e-6,
        -4.6869081420e-5,    1.1453049410e-5,    9.1849018650e-5,   -1.9507146130e-5,
        -1.6957646590e-4,    3.1618939830e-5,    2.9766702210e-4,   -4.9060345190e-5,
        -5.0020904750e-4,    7.3127761420e-5,    8.0895860450e-4,   -1.0506063700e-4,
        -1.2645865790e-3,    1.4598306730e-4,    1.9181631510e-3,   -1.9668265300e-4,
        -2.8329233170e-3,    2.5737608670e-4,    4.0868977080e-3,   -3.2757860030e-4,
        -5.7776579630e-3,    4.0610437280e-4,    8.0316131930e-3,   -4.9099023450e-4,
        -1.1022324670e-2,    5.7943403950e-4,    1.5008038840e-2,   -6.6789961420e-4,
        -2.0413324240e-2,    7.5244333130e-4,    2.8023539110e-2,   -8.2899280820e-4,
        -3.9516158400e-2,    8.9356512760e-4,    5.9269979600e-2,   -9.4253651330e-4,
        -1.0341136900e-1,    9.7307201940e-4,    3.1740292910e-1,    4.9901655320e-1,
    ];
    let expected_deviations = vec![
        6.2971832190e-9,    6.2971832190e-9,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.6401705930e2,
    ];
    let expected_extremal_frequencies = vec![
        7.3529412040e-3,    1.5271493230e-2,    2.2624434900e-2,    3.0542986470e-2,
        3.7895929070e-2,    4.5814480630e-2,    5.3733032200e-2,    6.1085972930e-2,
        6.9004528220e-2,    7.9751133920e-2,    8.8800907140e-2,    1.0124434530e-1,
        1.0690045360e-1,    1.1368778350e-1,    1.2104072420e-1,    1.4705911280e-1,
        1.5101844070e-1,    1.6233080630e-1,    1.6798698900e-1,    1.7703688140e-1,
        1.8212744590e-1,    1.9004610180e-1,    1.9400542970e-1,    1.9796475770e-1,
        1.9909599420e-1,    2.0000000300e-1,    3.0000001190e-1,    3.0056563020e-1,
        3.0339372160e-1,    3.0678743120e-1,    3.1187799570e-1,    3.1753417850e-1,
        3.2432159780e-1,    3.3167463540e-1,    3.3902767300e-1,    3.4694632890e-1,
        3.5543060300e-1,    3.6391487720e-1,    3.7239915130e-1,    3.8088342550e-1,
        3.8993331790e-1,    3.9898321030e-1,    4.0803310280e-1,    4.1708299520e-1,
        4.2613288760e-1,    4.3518278000e-1,    4.4423267250e-1,    4.5384818320e-1,
        4.6289807560e-1,    4.7194796800e-1,    4.8156347870e-1,    4.9061337110e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
104,1,2,3
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n104_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(104, FilterType::MultipleBand, &bands, 3);

    let expected_impulse_response = vec![
        2.2896514910e-7,   -2.0597082080e-7,   -1.1595847130e-6,    1.0854134870e-6,
        3.8644775490e-6,   -3.7852971670e-6,   -1.0310948710e-5,    1.0584360100e-5,
        2.3700144080e-5,   -2.5639397790e-5,   -4.8792608140e-5,    5.5916843850e-5,
        9.2061927720e-5,   -1.1224553600e-4,   -1.6148443680e-4,    2.1054749960e-4,
        2.6591686760e-4,   -3.7295633230e-4,   -4.1385125950e-4,    6.2872108540e-4,
        6.1145774090e-4,   -1.0147371800e-3,   -8.5991714150e-4,    1.5755303900e-3,
        1.1521347330e-3,   -2.3626128680e-3,   -1.4688915110e-3,    3.4334117080e-3,
        1.7747008240e-3,   -4.8499880360e-3,   -2.0132744680e-3,    6.6784285010e-3,
        2.1024968010e-3,   -8.9900642630e-3,   -1.9278107210e-3,    1.1867001650e-2,
        1.3317675330e-3,   -1.5416162090e-2,   -9.3815848230e-5,    1.9801303740e-2,
        -2.1135052670e-3,   -2.5315502660e-2,    5.8292746540e-3,    3.2558042560e-2,
        -1.2085906230e-2,   -4.2936477810e-2,    2.3319844160e-2,    6.0498893260e-2,
        -4.7425836320e-2,   -1.0253518820e-1,    1.3428519670e-1,    4.6444803480e-1,
    ];
    let expected_deviations = vec![
        5.2657551210e-8,    5.2657551210e-8,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.4557078550e2,
    ];
    let expected_extremal_frequencies = vec![
        3.2051282470e-3,    6.4102564940e-3,    1.2820512990e-2,    1.6025641930e-2,
        2.2435897960e-2,    3.5256411880e-2,    5.1282051950e-2,    5.4487179960e-2,
        6.0897435990e-2,    6.7307695750e-2,    7.6923079790e-2,    9.2948719860e-2,
        1.0576923190e-1,    1.1217948790e-1,    1.1538461600e-1,    1.2179487200e-1,
        1.3461540640e-1,    1.4102567730e-1,    1.5064108370e-1,    1.6025649010e-1,
        1.6666676100e-1,    1.7628216740e-1,    1.7948730290e-1,    1.8589757380e-1,
        1.9230784480e-1,    1.9551298020e-1,    2.0000000300e-1,    3.0000001190e-1,
        3.0320513250e-1,    3.0641025300e-1,    3.0961537360e-1,    3.1282049420e-1,
        3.2243585590e-1,    3.2884609700e-1,    3.3525633810e-1,    3.4166657920e-1,
        3.5128194090e-1,    3.6089730260e-1,    3.6730754380e-1,    3.7692290540e-1,
        3.8653826710e-1,    3.9294850830e-1,    4.0256387000e-1,    4.1217923160e-1,
        4.2179459330e-1,    4.3140995500e-1,    4.4102531670e-1,    4.5064067840e-1,
        4.5705091950e-1,    4.6666628120e-1,    4.7628164290e-1,    4.8589700460e-1,
        4.9551236630e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
104,1,2,4
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n104_grid4() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(104, FilterType::MultipleBand, &bands, 4);

    let expected_impulse_response = vec![
        -1.0287480020e-7,    1.6789925890e-7,    4.8228309880e-7,   -8.1139882010e-7,
        -1.6794635940e-6,    2.5429642390e-6,    4.7964595070e-6,   -6.3958186730e-6,
        -1.1978325350e-5,    1.3835310710e-5,    2.6955331120e-5,   -2.6728994270e-5,
        -5.5819557020e-5,    4.7075562180e-5,    1.0790002120e-4,   -7.6485681350e-5,
        -1.9673933280e-4,    1.1528180770e-4,    3.4102459900e-4,   -1.6124271500e-4,
        -5.6540558580e-4,    2.0794247390e-4,    9.0111425380e-4,   -2.4267358820e-4,
        -1.3862161430e-3,    2.4413829670e-4,    2.0656737030e-3,   -1.7979030960e-4,
        -2.9912274330e-3,    2.9501970860e-6,    4.2216754520e-3,    3.5068183210e-4,
        -5.8243023230e-3,   -9.6781249160e-4,    7.8790485860e-3,    1.9653968050e-3,
        -1.0488335970e-2,   -3.5053547470e-3,    1.3798735100e-2,    5.8251335290e-3,
        -1.8048500640e-2,   -9.3065164980e-3,    2.3678109050e-2,    1.4646789060e-2,
        -3.1619422140e-2,   -2.3353457450e-2,    4.4214684520e-2,    3.9568264040e-2,
        -6.9279134270e-2,   -8.0660298470e-2,    1.5663644670e-1,    4.4208961730e-1,
    ];
    let expected_deviations = vec![
        4.4436514910e-8,    4.4436514910e-8,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.4704519650e2,
    ];
    let expected_extremal_frequencies = vec![
        2.4038462430e-3,    7.2115389630e-3,    1.9230769950e-2,    2.4038461970e-2,
        3.3653847870e-2,    4.5673087240e-2,    6.0096174480e-2,    6.4903870220e-2,
        6.9711565970e-2,    7.9326957460e-2,    8.4134653210e-2,    9.6153892580e-2,
        1.0336543620e-1,    1.0817313190e-1,    1.1538467560e-1,    1.2500005960e-1,
        1.2740390000e-1,    1.3221158090e-1,    1.4182694260e-1,    1.5384614470e-1,
        1.6346150640e-1,    1.6586534680e-1,    1.7307686810e-1,    1.8749991060e-1,
        1.9230759140e-1,    1.9471143190e-1,    1.9711527230e-1,    2.0000000300e-1,
        3.0000001190e-1,    3.0240386720e-1,    3.0480772260e-1,    3.1201928850e-1,
        3.1682699920e-1,    3.2403856520e-1,    3.3125013110e-1,    3.4086555240e-1,
        3.4807711840e-1,    3.5769253970e-1,    3.6730796100e-1,    3.7451952700e-1,
        3.8413494830e-1,    3.9375036950e-1,    4.0336579080e-1,    4.1057735680e-1,
        4.2019277810e-1,    4.2980819940e-1,    4.3942362070e-1,    4.4903904200e-1,
        4.5865446330e-1,    4.6826988460e-1,    4.7548145060e-1,    4.8509687190e-1,
        4.9471229310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
104,1,2,16
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n104_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(104, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        2.1663130670e-7,   -3.5653295070e-7,   -1.3183645250e-6,    1.4412797780e-6,
        4.4942380550e-6,   -4.0726445150e-6,   -1.1467403970e-5,    1.0211058000e-5,
        2.5163753890e-5,   -2.3149767000e-5,   -4.9668778960e-5,    4.8037316450e-5,
        8.9690467580e-5,   -9.3188777100e-5,   -1.5040210560e-4,    1.7102592390e-4,
        2.3679574950e-4,   -2.9865963730e-4,   -3.5188772020e-4,    4.9900868910e-4,
        4.9468333600e-4,   -8.0197933130e-4,   -6.5839145100e-4,    1.2444294990e-3,
        8.2753208700e-4,   -1.8700866490e-3,   -9.7455218200e-4,    2.7299136850e-3,
        1.0571592720e-3,   -3.8818335160e-3,   -1.0152794420e-3,    5.3907171820e-3,
        7.6668732800e-4,   -7.3309051800e-3,   -2.0137289540e-4,    9.7926370800e-3,
        -8.2797417420e-4,   -1.2895923110e-2,    2.5263936260e-3,    1.6822144390e-2,
        -5.2014524120e-3,   -2.1886032070e-2,    9.3688471240e-3,    2.8710866350e-2,
        -1.6035865990e-2,   -3.8728259500e-2,    2.7614692230e-2,    5.6008655580e-2,
        -5.1978770640e-2,   -9.7860582170e-2,    1.3899274170e-1,    4.5969921350e-1,
    ];
    let expected_deviations = vec![
        6.7786906750e-8,    6.7786906750e-8,
    ];
    let expected_deviation_dbs = vec![
        1.0354386860e-6,    -1.4337709050e2,
    ];
    let expected_extremal_frequencies = vec![
        1.8028847410e-3,    6.6105769950e-3,    1.8028853460e-2,    1.9230777400e-2,
        2.0432701330e-2,    3.1850975010e-2,    3.7259615960e-2,    4.6874977650e-2,
        5.4086498920e-2,    6.1298020180e-2,    6.8509578700e-2,    7.7524036170e-2,
        8.5937529800e-2,    9.2548131940e-2,    1.1057704690e-1,    1.2319728730e-1,
        1.3281270860e-1,    1.4062523840e-1,    1.4783680440e-1,    1.5745222570e-1,
        1.6526475550e-1,    1.7427921300e-1,    1.8269270660e-1,    1.8990427260e-1,
        1.9531294700e-1,    1.9891873000e-1,    2.0000000300e-1,    3.0000001190e-1,
        3.0060097580e-1,    3.0360579490e-1,    3.0721157790e-1,    3.1201928850e-1,
        3.1862989070e-1,    3.2584145670e-1,    3.3365398650e-1,    3.4266844390e-1,
        3.5168290140e-1,    3.6009639500e-1,    3.6911085250e-1,    3.7812530990e-1,
        3.8713976740e-1,    3.9615422490e-1,    4.0516868230e-1,    4.1418313980e-1,
        4.2259663340e-1,    4.3221205470e-1,    4.4122651220e-1,    4.5024096970e-1,
        4.5925542710e-1,    4.6826988460e-1,    4.7728434210e-1,    4.8629879950e-1,
        4.9531325700e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
104,1,2,17
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n104_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(104, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        1.3903002250e-3,   -6.8450719120e-4,   -2.4218654730e-3,    2.0593712110e-3,
        1.4871268070e-3,   -3.5011963920e-3,   -5.5479630830e-4,    5.2549983370e-3,
        1.1073874780e-3,   -7.9925190660e-3,   -6.3133169900e-3,    1.3232162220e-2,
        2.2084806110e-2,   -2.3987092080e-2,   -5.8291181920e-2,    4.5645140110e-2,
        1.2989920380e-1,   -8.7010852990e-2,   -2.5771132110e-1,    1.6134405140e-1,
        4.6832114460e-1,   -2.8713268040e-1,   -7.9293358330e-1,    4.8827028270e-1,
        1.2648258210e0,   -7.9330015180e-1,    -1.9154431820e0,     1.2334649560e0,
        2.7694144250e0,    -1.8394664530e0,    -3.8390636440e0,     2.6370840070e0,
        5.1192469600e0,    -3.6420803070e0,    -6.5834593770e0,     4.8550834660e0,
        8.1821222310e0,    -6.2573213580e0,    -9.8437204360e0,     7.8081183430e0,
        1.1479064940e1,    -9.4449510570e0,    -1.2988411900e1,     1.1086581230e1,
        1.4270594600e1,    -1.2639521600e1,    -1.5232538220e1,     1.4008321760e1,
        1.5795203210e1,    -1.5116160390e1,    -1.5865180020e1,     1.6179929730e1,
    ];
    let expected_deviations = vec![
        5.4373927580e-8,    5.4373927580e-8,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.4529219060e2,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.7873303780e-3,    1.4705882410e-2,    1.5837104990e-2,
        2.2624434900e-2,    2.9977375640e-2,    3.9027150720e-2,    4.5814480630e-2,
        5.4298643020e-2,    5.9389140460e-2,    6.9004528220e-2,    7.6357468960e-2,
        8.3144798870e-2,    8.8235296310e-2,    1.0011312370e-1,    1.0690045360e-1,
        1.1595022680e-1,    1.2104072420e-1,    1.3291865590e-1,    1.3970607520e-1,
        1.4819034930e-1,    1.5610900520e-1,    1.6459327940e-1,    1.7307755350e-1,
        1.8099620940e-1,    1.9004610180e-1,    1.9683352110e-1,    2.0000000300e-1,
        3.0169686680e-1,    3.0791866780e-1,    3.1696856020e-1,    3.2545283440e-1,
        3.3393710850e-1,    3.4242138270e-1,    3.5034003850e-1,    3.5825869440e-1,
        3.6617735030e-1,    3.7409600620e-1,    3.8201466200e-1,    3.8993331790e-1,
        3.9785197380e-1,    4.0577062960e-1,    4.1368928550e-1,    4.2160794140e-1,
        4.2952659730e-1,    4.3744525310e-1,    4.4536390900e-1,    4.5384818320e-1,
        4.6176683900e-1,    4.6968549490e-1,    4.7816976900e-1,    4.8778527980e-1,
        4.9626955390e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,1,2,3
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n128_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 3);

    let expected_impulse_response = vec![
        1.2134856660e-5,   -9.4922415880e-7,   -2.4764700360e-5,    5.5226810220e-6,
        2.7363032130e-5,   -2.1574353010e-5,   -3.6218410970e-5,    6.9269728560e-5,
        6.1583821660e-5,   -1.9239273390e-4,   -1.2564715870e-4,    4.7662347790e-4,
        2.7269998100e-4,   -1.0760027220e-3,   -5.8501731840e-4,    2.2480832410e-3,
        1.2058158170e-3,   -4.3967561800e-3,   -2.3700294550e-3,    8.1193810330e-3,
        4.4422955250e-3,   -1.4252329250e-2,   -7.9597299920e-3,    2.3906450720e-2,
        1.3674821700e-2,   -3.8481783120e-2,   -2.2591151300e-2,    5.9650205080e-2,
        3.5982605070e-2,   -8.9295715090e-2,   -5.5385492740e-2,    1.2940533460e-1,
        8.2553319630e-2,   -1.8190917370e-1,   -1.1936606470e-1,    2.4847579000e-1,
        1.6769063470e-1,   -3.3027726410e-1,   -2.2919505830e-1,    4.2774692180e-1,
        3.0512726310e-1,   -5.4035890100e-1,   -3.9607670900e-1,    6.6646236180e-1,
        5.0174450870e-1,   -8.0320245030e-1,   -6.2075066570e-1,    9.4655483960e-1,
        7.5050759320e-1,    -1.0914915800e0,   -8.8718062640e-1,     1.2322880030e0,
        1.0257409810e0,    -1.3629758360e0,    -1.1600804330e0,     1.4779745340e0,
        1.2830615040e0,    -1.5730798240e0,    -1.3860292430e0,     1.6477599140e0,
        1.4552294020e0,    -1.7161500450e0,    -1.4357572790e0,     2.0721986290e0,
    ];
    let expected_deviations = vec![
        4.6792663970e-8,    4.6792663970e-8,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.4659643550e2,
    ];
    let expected_extremal_frequencies = vec![
        2.6041667440e-3,    7.8125000000e-3,    1.5625000000e-2,    2.0833332090e-2,
        3.6458332090e-2,    3.9062500000e-2,    4.4270835820e-2,    4.6875003730e-2,
        5.2083339540e-2,    5.4687507450e-2,    6.2500007450e-2,    7.0312500000e-2,
        7.5520828370e-2,    8.0729156730e-2,    8.5937485100e-2,    9.8958306010e-2,
        1.0416663440e-1,    1.1197912690e-1,    1.1979161950e-1,    1.2499994780e-1,
        1.3281245530e-1,    1.3802079860e-1,    1.4322914180e-1,    1.5364582840e-1,
        1.6145834330e-1,    1.6666668650e-1,    1.7187502980e-1,    1.7447920140e-1,
        1.7968754470e-1,    1.8489588800e-1,    1.9010423120e-1,    1.9531257450e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0260416870e-1,    3.0520832540e-1,
        3.0781248210e-1,    3.1041663890e-1,    3.1562495230e-1,    3.2083326580e-1,
        3.2604157920e-1,    3.3124989270e-1,    3.3645820620e-1,    3.4427067640e-1,
        3.4947898980e-1,    3.5729146000e-1,    3.6249977350e-1,    3.7031224370e-1,
        3.7812471390e-1,    3.8593718410e-1,    3.9114549760e-1,    3.9895796780e-1,
        4.0677043800e-1,    4.1458290820e-1,    4.2239537840e-1,    4.2760369180e-1,
        4.3541616200e-1,    4.4322863220e-1,    4.5104110240e-1,    4.5885357260e-1,
        4.6666604280e-1,    4.7447851300e-1,    4.8229098320e-1,    4.8749929670e-1,
        4.9531176690e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
128,1,2,4
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n128_grid4() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 4);

    let expected_impulse_response = vec![
        -4.9079631030e-6,   -3.7006091130e-6,    9.9486933320e-6,    1.3435306760e-5,
        -1.0668485630e-5,   -3.2348762030e-5,    1.3592200050e-5,    7.3866372990e-5,
        -2.3204051100e-5,   -1.6391406830e-4,    5.0092647140e-5,    3.4787313780e-4,
        -1.1592023660e-4,   -6.9879757940e-4,    2.6045468990e-4,    1.3256946110e-3,
        -5.5052136300e-4,   -2.3799987980e-3,    1.0901981730e-3,    4.0578329940e-3,
        -2.0306762310e-3,   -6.5952185540e-3,    3.5776302680e-3,    1.0253731160e-2,
        -5.9931995350e-3,   -1.5294963490e-2,    9.5894997940e-3,    2.1943699570e-2,
        -1.4710919930e-2,   -3.0341936280e-2,    2.1703565490e-2,    4.0498323740e-2,
        -3.0872121450e-2,   -5.2239988000e-2,    4.2426623400e-2,    6.5175235270e-2,
        -5.6424677370e-2,   -7.8676424920e-2,    7.2716735300e-2,    9.1890782120e-2,
        -9.0904206040e-2,   -1.0378492620e-1,    1.1032041160e-1,    1.1322417860e-1,
        -1.3004343210e-1,   -1.1908350140e-1,    1.4894768600e-1,    1.2038206310e-1,
        -1.6579747200e-1,   -1.1643116920e-1,    1.7938363550e-1,    1.0698707400e-1,
        -1.8870645760e-1,   -9.2413499950e-2,    1.9322699310e-1,    7.3902681470e-2,
        -1.9328407940e-1,   -5.3965307770e-2,    1.9111025330e-1,    3.8192130630e-2,
        -1.9494476910e-1,   -4.5731194320e-2,    2.5556340810e-1,    3.7399479750e-1,
    ];
    let expected_deviations = vec![
        5.8612602770e-8,    5.8612602770e-8,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.4464018250e2,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    9.7656250000e-3,    1.7578125000e-2,    1.9531250000e-2,
        2.5390625000e-2,    2.7343750000e-2,    3.7109375000e-2,    4.6875000000e-2,
        5.0781250000e-2,    5.4687500000e-2,    6.2500000000e-2,    6.4453125000e-2,
        6.8359375000e-2,    7.4218750000e-2,    8.5937500000e-2,    8.7890625000e-2,
        1.0156250000e-1,    1.0937500000e-1,    1.1523437500e-1,    1.2109375000e-1,
        1.2890625000e-1,    1.3671875000e-1,    1.4257812500e-1,    1.5039062500e-1,
        1.5820312500e-1,    1.6406250000e-1,    1.7187500000e-1,    1.7773437500e-1,
        1.8359375000e-1,    1.8750000000e-1,    1.9140625000e-1,    1.9531250000e-1,
        1.9726562500e-1,    2.0000000300e-1,    3.0195313690e-1,    3.0585938690e-1,
        3.0781251190e-1,    3.0976563690e-1,    3.1562501190e-1,    3.2148438690e-1,
        3.2734376190e-1,    3.3320313690e-1,    3.3906251190e-1,    3.4492188690e-1,
        3.5078126190e-1,    3.5664063690e-1,    3.6445313690e-1,    3.7031251190e-1,
        3.7812501190e-1,    3.8593751190e-1,    3.9375001190e-1,    3.9960938690e-1,
        4.0742188690e-1,    4.1523438690e-1,    4.2304688690e-1,    4.2890626190e-1,
        4.3671876190e-1,    4.4453126190e-1,    4.5234376190e-1,    4.5820313690e-1,
        4.6601563690e-1,    4.7382813690e-1,    4.8164063690e-1,    4.8945313690e-1,
        4.9531251190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
128,1,2,16
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n128_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -1.5079915270e-7,   -8.1042298920e-7,   -1.5023859990e-6,    3.8520993260e-6,
        9.4089100460e-6,   -1.2349444660e-5,   -3.1739571570e-5,    3.3966200140e-5,
        8.5824918640e-5,   -8.3969353000e-5,   -2.0350694830e-4,    1.8991598340e-4,
        4.3878558790e-4,   -3.9816653590e-4,   -8.7712879760e-4,    7.8268261860e-4,
        1.6465839000e-3,   -1.4552678910e-3,   -2.9294118290e-3,    2.5764643210e-3,
        4.9723722040e-3,   -4.3661389500e-3,   -8.0934744330e-3,    7.1126110850e-3,
        1.2683628130e-2,   -1.1177381500e-2,   -1.9200244920e-2,    1.6993086790e-2,
        2.8150238100e-2,   -2.5052927430e-2,   -4.0061220530e-2,    3.5890024160e-2,
        5.5441349740e-2,   -5.0045244400e-2,   -7.4728794400e-2,    6.8023853000e-2,
        9.8233424130e-2,   -9.0243846180e-2,   -1.2607598300e-1,    1.1697976290e-1,
        1.5813110770e-1,   -1.4830701050e-1,   -1.9398063420e-1,    1.8405362960e-1,
        2.3288273810e-1,   -2.2376818950e-1,   -2.7376261350e-1,    2.6671147350e-1,
        3.1522560120e-1,   -3.1188157200e-1,   -3.5558861490e-1,    3.5808292030e-1,
        3.9291274550e-1,   -4.0406361220e-1,   -4.2499640580e-1,    4.4878020880e-1,
        4.4920447470e-1,   -4.9200463290e-1,   -4.6167188880e-1,    5.3626894950e-1,
        4.5335510370e-1,   -5.9757769110e-1,   -3.7359198930e-1,    9.7037732600e-1,
    ];
    let expected_deviations = vec![
        6.2998914530e-8,    6.2998914530e-8,
    ];
    let expected_deviation_dbs = vec![
        1.0354386860e-6,    -1.4401333620e2,
    ];
    let expected_extremal_frequencies = vec![
        2.9296875000e-3,    5.8593750000e-3,    1.3671875000e-2,    1.6113281250e-2,
        2.6855468750e-2,    2.7343750000e-2,    3.6621093750e-2,    4.2480468750e-2,
        4.9316406250e-2,    5.4687500000e-2,    6.2500000000e-2,    6.7382812500e-2,
        7.5195312500e-2,    8.0566406250e-2,    8.5937500000e-2,    9.2773437500e-2,
        9.6191406250e-2,    1.0595703120e-1,    1.1132812500e-1,    1.2695312500e-1,
        1.3623046880e-1,    1.4501953120e-1,    1.5087890620e-1,    1.5771484380e-1,
        1.6601562500e-1,    1.7285156250e-1,    1.7822265620e-1,    1.8554687500e-1,
        1.8896484380e-1,    1.9580078120e-1,    1.9726562500e-1,    1.9873046880e-1,
        2.0000000300e-1,    3.0000001190e-1,    3.0048829320e-1,    3.0244141820e-1,
        3.0537110570e-1,    3.0927735570e-1,    3.1367188690e-1,    3.1855469940e-1,
        3.2392579320e-1,    3.3027344940e-1,    3.3662110570e-1,    3.4345704320e-1,
        3.4980469940e-1,    3.5664063690e-1,    3.6396485570e-1,    3.7080079320e-1,
        3.7763673070e-1,    3.8496094940e-1,    3.9228516820e-1,    3.9960938690e-1,
        4.0693360570e-1,    4.1425782440e-1,    4.2158204320e-1,    4.2890626190e-1,
        4.3623048070e-1,    4.4404298070e-1,    4.5136719940e-1,    4.5869141820e-1,
        4.6650391820e-1,    4.7382813690e-1,    4.8115235570e-1,    4.8896485570e-1,
        4.9628907440e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
128,1,2,17
0.000000000000,0.200000002980,0.300000011921,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_halfband_n128_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.2f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.3f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        -2.9506096240e-1,   -2.1212249990e-2,    5.8941400050e-1,    6.4105123280e-2,
        -5.8710730080e-1,   -1.0826389490e-1,    5.8260315660e-1,    1.5398690100e-1,
        -5.7475268840e-1,   -2.0046664770e-1,    5.6169056890e-1,    2.4479238690e-1,
        -5.4074454310e-1,   -2.8048068280e-1,    5.0855582950e-1,    2.9553884270e-1,
        -4.6155744790e-1,   -2.7016532420e-1,    3.9695847030e-1,    1.7431899910e-1,
        -3.1436061860e-1,    3.4455463290e-2,    2.1804933250e-1,   -4.1243261100e-1,
        -1.1986750360e-1,     1.0297607180e0,    4.2422235010e-2,    -1.9686007500e0,
        -2.2201061250e-2,     3.3189015390e0,    1.1200821400e-1,    -5.1716127400e0,
        -3.8202762600e-1,     7.6095027920e0,    9.1882276540e-1,    -1.0696158410e1,
        -1.8217124940e0,     1.4464180950e1,     3.1962008480e0,    -1.8903953550e1,
        -5.1445565220e0,     2.3954568860e1,     7.7540636060e0,    -2.9498527530e1,
        -1.1084004400e1,     3.5361541750e1,     1.5152784350e1,    -4.1318332670e1,
        -1.9926954270e1,     4.7104545590e1,     2.5313835140e1,    -5.2434234620e1,
        -3.1159355160e1,     5.7021335600e1,     3.7252147670e1,    -6.0603088380e1,
        -4.3334381100e1,     6.2962570190e1,     4.9119308470e1,    -6.3946720120e1,
        -5.4316658020e1,     6.3470027920e1,     5.8693199160e1,    -6.1256652830e1,
    ];
    let expected_deviations = vec![
        6.5447551380e-8,    6.5447551380e-8,
    ];
    let expected_deviation_dbs = vec![
        1.0354386860e-6,    -1.4368212890e2,
    ];
    let expected_extremal_frequencies = vec![
        3.2169118060e-3,    8.2720592620e-3,    1.1029414830e-2,    1.2408092620e-2,
        1.8382363020e-2,    2.4816192690e-2,    3.1709581610e-2,    3.6764733490e-2,
        4.5036800210e-2,    4.9632392820e-2,    5.6525781750e-2,    6.2040492890e-2,
        6.8474322560e-2,    7.4908152220e-2,    8.1341981890e-2,    8.8694930080e-2,
        9.1911844910e-2,    1.0202214870e-1,    1.1305157100e-1,    1.2224275620e-1,
        1.2959562240e-1,    1.3648889960e-1,    1.4338217680e-1,    1.5027545390e-1,
        1.5670917930e-1,    1.6360245650e-1,    1.7049573360e-1,    1.7738901080e-1,
        1.8382273610e-1,    1.9025646150e-1,    1.9577108320e-1,    1.9898794590e-1,
        2.0000000300e-1,    3.0045956370e-1,    3.0275732280e-1,    3.0781239270e-1,
        3.1424611810e-1,    3.2067984340e-1,    3.2757312060e-1,    3.3354729410e-1,
        3.3906191590e-1,    3.4365743400e-1,    3.4825295210e-1,    3.5376757380e-1,
        3.5928219560e-1,    3.6525636910e-1,    3.7123054270e-1,    3.7812381980e-1,
        3.8455754520e-1,    3.9099127050e-1,    3.9788454770e-1,    4.0477782490e-1,
        4.1167110200e-1,    4.1856437920e-1,    4.2545765640e-1,    4.3281048540e-1,
        4.3970376250e-1,    4.4659703970e-1,    4.5394986870e-1,    4.6084314580e-1,
        4.6773642300e-1,    4.7508925200e-1,    4.8198252920e-1,    4.8933535810e-1,
        4.9622863530e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
5,1,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n5_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(5, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        1.7225147780e-1,    2.7264657620e-1,    2.8129759430e-1,
    ];
    let expected_deviations = vec![
        1.7109370230e-1,    1.7109370230e-1,
    ];
    let expected_deviation_dbs = vec![
        1.3718329670e0,    -1.5335319520e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0000000150e-1,    2.0000000300e-1,    3.1458330150e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
5,1,2,17
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n5_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(5, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        1.7220878600e-1,    2.7275833490e-1,    2.8115129470e-1,
    ];
    let expected_deviations = vec![
        1.7108555140e-1,    1.7108555140e-1,
    ];
    let expected_deviation_dbs = vec![
        1.3717728850e0,    -1.5335733410e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0000000150e-1,    2.0000000300e-1,    3.1764706970e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
5,1,2,23
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n5_grid23() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(5, FilterType::MultipleBand, &bands, 23);

    let expected_impulse_response = vec![
        1.7224420610e-1,    2.7266561990e-1,    2.8127267960e-1,
    ];
    let expected_deviations = vec![
        1.7109231650e-1,    1.7109231650e-1,
    ];
    let expected_deviation_dbs = vec![
        1.3718223570e0,    -1.5335390090e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0000000150e-1,    2.0000000300e-1,    3.1594201920e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
5,1,2,24
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n5_grid24() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(5, FilterType::MultipleBand, &bands, 24);

    let expected_impulse_response = vec![
        1.7219616470e-1,    2.7279138570e-1,    2.8110802170e-1,
    ];
    let expected_deviations = vec![
        1.7108313740e-1,    1.7108313740e-1,
    ];
    let expected_deviation_dbs = vec![
        1.3717542890e0,    -1.5335855480e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0000000150e-1,    2.0000000300e-1,    3.1805559990e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
6,1,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n6_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(6, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        5.7075034830e-2,    2.5446560980e-1,    2.7620327470e-1,
    ];
    let expected_deviations = vec![
        1.7548783120e-1,    1.7548783120e-1,
    ];
    let expected_deviation_dbs = vec![
        1.4043630360e0,    -1.5115059850e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0000000150e-1,    2.0000000300e-1,    3.3541661500e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
6,1,2,17
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n6_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(6, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        5.7066582140e-2,    2.5447952750e-1,    2.7619627120e-1,
    ];
    let expected_deviations = vec![
        1.7548479140e-1,    1.7548479140e-1,
    ];
    let expected_deviation_dbs = vec![
        1.4043402670e0,    -1.5115211490e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0000000150e-1,    2.0000000300e-1,    3.2745099070e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
6,1,2,23
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n6_grid23() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(6, FilterType::MultipleBand, &bands, 23);

    let expected_impulse_response = vec![
        5.7109892370e-2,    2.5440812110e-1,    2.7623218300e-1,
    ];
    let expected_deviations = vec![
        1.7550042270e-1,    1.7550042270e-1,
    ];
    let expected_deviation_dbs = vec![
        1.4044556620e0,    -1.5114436150e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0000000150e-1,    2.0000000300e-1,    3.3043476940e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
6,1,2,24
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n6_grid24() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(6, FilterType::MultipleBand, &bands, 24);

    let expected_impulse_response = vec![
        5.7112943380e-2,    2.5440308450e-1,    2.7623471620e-1,
    ];
    let expected_deviations = vec![
        1.7550152540e-1,    1.7550152540e-1,
    ];
    let expected_deviation_dbs = vec![
        1.4044644830e0,    -1.5114381790e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0000000150e-1,    2.0000000300e-1,    3.3194449540e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n16_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        1.6626054420e-2,   -6.9641619920e-3,   -3.4036625180e-2,   -4.8550553620e-2,
        -1.4346834270e-2,    8.0486699940e-2,    2.0301045480e-1,    2.8957736490e-1,
    ];
    let expected_deviations = vec![
        2.8395211320e-2,    2.8395211320e-2,
    ];
    let expected_deviation_dbs = vec![
        2.4320062990e-1,    -3.0935096740e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    7.0312500000e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.2343750300e-1,    2.7812498810e-1,    3.4062498810e-1,    4.0312498810e-1,
        4.6562498810e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,17
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n16_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        1.6604714100e-2,   -6.9665480410e-3,   -3.4042973070e-2,   -4.8545256260e-2,
        -1.4339629560e-2,    8.0507189040e-2,    2.0299847420e-1,    2.8958114980e-1,
    ];
    let expected_deviations = vec![
        2.8405755760e-2,    2.8405755760e-2,
    ];
    let expected_deviation_dbs = vec![
        2.4329023060e-1,    -3.0931873320e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.9852948190e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.2573532160e-1,    2.7720594410e-1,    3.3970600370e-1,    4.0220606330e-1,
        4.6838259700e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,23
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n16_grid23() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 23);

    let expected_impulse_response = vec![
        1.6619490460e-2,   -6.9735571740e-3,   -3.4042555840e-2,   -4.8542119560e-2,
        -1.4334453270e-2,    8.0497525630e-2,    2.0297560100e-1,    2.8958702090e-1,
    ];
    let expected_deviations = vec![
        2.8426080940e-2,    2.8426080940e-2,
    ];
    let expected_deviation_dbs = vec![
        2.4346140030e-1,    -3.0925659180e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    7.0652164520e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.2445651890e-1,    2.7880448100e-1,    3.3858740330e-1,    4.0380513670e-1,
        4.6630546450e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,24
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n16_grid24() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 24);

    let expected_impulse_response = vec![
        1.6611199830e-2,   -6.9524645810e-3,   -3.4039940690e-2,   -4.8550907520e-2,
        -1.4349140230e-2,    8.0499984320e-2,    2.0299473400e-1,    2.8958344460e-1,
    ];
    let expected_deviations = vec![
        2.8406174850e-2,    2.8406174850e-2,
    ];
    let expected_deviation_dbs = vec![
        2.4329325560e-1,    -3.0931743620e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    7.0312500000e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.2343754770e-1,    2.7812498810e-1,    3.3802059290e-1,    4.0312451120e-1,
        4.6822842960e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
43,1,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n43_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(43, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        2.4129023950e-4,    4.4630825870e-5,   -6.1199563790e-4,   -1.2621002970e-3,
        -7.3168548990e-4,    1.5280211810e-3,    3.9739278150e-3,    3.3880707340e-3,
        -1.9927457910e-3,   -9.0329209340e-3,   -1.0291606190e-2,   -4.0036690190e-4,
        1.6221683470e-2,    2.4718172850e-2,    1.0618302040e-2,   -2.4088792500e-2,
        -5.3955029700e-2,   -4.2595136910e-2,    3.0296690760e-2,    1.4724741880e-1,
        2.5623118880e-1,    3.0067047480e-1,
    ];
    let expected_deviations = vec![
        2.3551993950e-4,    2.3551993950e-4,
    ];
    let expected_deviation_dbs = vec![
        2.0457860080e-3,    -7.2559448240e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.2727277130e-2,    4.5454539360e-2,    6.5340884030e-2,
        8.3806775510e-2,    9.5170401040e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.0426136260e-1,    2.1562498810e-1,    2.3267042640e-1,    2.5113633280e-1,
        2.7244335410e-1,    2.9375037550e-1,    3.1647786500e-1,    3.3920535450e-1,
        3.6193284390e-1,    3.8466033340e-1,    4.0738782290e-1,    4.3011531230e-1,
        4.5426326990e-1,    4.7699075940e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
43,1,2,17
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n43_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(43, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        2.4092210510e-4,    4.4559055820e-5,   -6.1146239750e-4,   -1.2610098350e-3,
        -7.3084910400e-4,    1.5275449260e-3,    3.9719203490e-3,    3.3859610560e-3,
        -1.9931907300e-3,   -9.0309465300e-3,   -1.0288252500e-2,   -3.9826903960e-4,
        1.6220826660e-2,    2.4714753030e-2,    1.0614889670e-2,   -2.4089926850e-2,
        -5.3953111170e-2,   -4.2591419070e-2,    3.0299793930e-2,    1.4724789560e-1,
        2.5622850660e-1,    3.0066657070e-1,
    ];
    let expected_deviations = vec![
        2.3516784130e-4,    2.3516784130e-4,
    ];
    let expected_deviation_dbs = vec![
        2.0426805130e-3,    -7.2572441100e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.2727275270e-2,    4.5454550530e-2,    6.5508022900e-2,
        8.2887656990e-2,    9.4919711350e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.0401071010e-1,    2.1604283150e-1,    2.3208566010e-1,    2.5213918090e-1,
        2.7219271660e-1,    2.9358315470e-1,    3.1631049510e-1,    3.3903783560e-1,
        3.6176517610e-1,    3.8449251650e-1,    4.0721985700e-1,    4.2994719740e-1,
        4.5401144030e-1,    4.7673878070e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
43,1,2,23
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n43_grid23() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(43, FilterType::MultipleBand, &bands, 23);

    let expected_impulse_response = vec![
        2.4126203790e-4,    4.5265009250e-5,   -6.1105226630e-4,   -1.2618330070e-3,
        -7.3328893630e-4,    1.5251546870e-3,    3.9724740200e-3,    3.3904986920e-3,
        -1.9872609990e-3,   -9.0292450040e-3,   -1.0294138450e-2,   -4.0859534060e-4,
        1.6214596110e-2,    2.4719521400e-2,    1.0628590360e-2,   -2.4077795450e-2,
        -5.3953576830e-2,   -4.2605768890e-2,    3.0282547700e-2,    1.4724227790e-1,
        2.5623986120e-1,    3.0068576340e-1,
    ];
    let expected_deviations = vec![
        2.3527225130e-4,    2.3527225130e-4,
    ];
    let expected_deviation_dbs = vec![
        2.0437156780e-3,    -7.2568588260e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.2727277130e-2,    4.5454528180e-2,    6.6205486660e-2,
        8.3003878590e-2,    9.5849707720e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.0395256580e-1,    2.1581025420e-1,    2.3260864620e-1,    2.5138333440e-1,
        2.7213460210e-1,    2.9387402530e-1,    3.1660160420e-1,    3.3834102750e-1,
        3.6106860640e-1,    3.8478434090e-1,    4.0751191970e-1,    4.3023949860e-1,
        4.5395523310e-1,    4.7668281200e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
43,1,2,24
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n43_grid24() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(43, FilterType::MultipleBand, &bands, 24);

    let expected_impulse_response = vec![
        2.4161032340e-4,    4.5395419870e-5,   -6.1145058130e-4,   -1.2627495450e-3,
        -7.3404295840e-4,    1.5253627210e-3,    3.9738481860e-3,    3.3921648280e-3,
        -1.9866686780e-3,   -9.0304585170e-3,   -1.0296533820e-2,   -4.1045487160e-4,
        1.6214815900e-2,    2.4721838530e-2,    1.0631569660e-2,   -2.4076398460e-2,
        -5.3954940290e-2,   -4.2608961460e-2,    3.0279736970e-2,    1.4724186060e-1,
        2.5624206660e-1,    3.0068916080e-1,
    ];
    let expected_deviations = vec![
        2.3564674480e-4,    2.3564674480e-4,
    ];
    let expected_deviation_dbs = vec![
        2.0468214060e-3,    -7.2554771420e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.3674234750e-2,    4.5454517010e-2,    6.6287830470e-2,
        8.3333268760e-2,    9.5643863080e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.0473484690e-1,    2.1609847250e-1,    2.3219694200e-1,    2.5208330150e-1,
        2.7291661500e-1,    2.9374992850e-1,    3.1647717950e-1,    3.3920443060e-1,
        3.6193168160e-1,    3.8465893270e-1,    4.0738618370e-1,    4.3011343480e-1,
        4.5378765460e-1,    4.7651490570e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
44,1,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n44_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(44, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        2.5466995430e-4,    2.3565371520e-4,   -3.0157432780e-4,   -1.1912961490e-3,
        -1.3641365800e-3,    2.7868960750e-4,    3.1936042940e-3,    4.5401183890e-3,
        1.2596715240e-3,   -6.0329320840e-3,   -1.1279707770e-2,   -6.9000981750e-3,
        8.0819651480e-3,    2.3027777670e-2,    2.1152526140e-2,   -5.3891157730e-3,
        -4.2278151960e-2,   -5.5633954700e-2,   -1.3684924690e-2,    8.6157023910e-2,
        2.0672370490e-1,    2.8906264900e-1,
    ];
    let expected_deviations = vec![
        1.7569570630e-4,    1.7569570630e-4,
    ];
    let expected_deviation_dbs = vec![
        1.5261025400e-3,    -7.5104774480e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.5568187240e-2,    4.8295445740e-2,    6.8181790410e-2,
        8.3806775510e-2,    9.5170401040e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.0426136260e-1,    2.1562498810e-1,    2.3124997320e-1,    2.4971586470e-1,
        2.6818194990e-1,    2.8948897120e-1,    3.1079599260e-1,    3.3210301400e-1,
        3.5483050350e-1,    3.7613752480e-1,    3.9886501430e-1,    4.2159250380e-1,
        4.4431999330e-1,    4.6562701460e-1,    4.8835450410e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
44,1,2,17
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n44_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(44, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        2.5493226710e-4,    2.3592231450e-4,   -3.0151745890e-4,   -1.1916277000e-3,
        -1.3650400800e-3,    2.7779288940e-4,    3.1934953290e-3,    4.5415489000e-3,
        1.2618808540e-3,   -6.0317083260e-3,   -1.1281008830e-2,   -6.9036218340e-3,
        8.0788834020e-3,    2.3027811200e-2,    2.1156564350e-2,   -5.3839418110e-3,
        -4.2276132850e-2,   -5.5637210610e-2,   -1.3691391800e-2,    8.6152493950e-2,
        2.0672501620e-1,    2.8906881810e-1,
    ];
    let expected_deviations = vec![
        1.7605912580e-4,    1.7605912580e-4,
    ];
    let expected_deviation_dbs = vec![
        1.5292083840e-3,    -7.5086830140e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.5401072580e-2,    4.8128347840e-2,    6.8181812760e-2,
        8.4224551920e-2,    9.6256606280e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.0401071010e-1,    2.1470592920e-1,    2.3074875770e-1,    2.4946539100e-1,
        2.6818200950e-1,    2.8957244750e-1,    3.1096288560e-1,    3.3235332370e-1,
        3.5508066420e-1,    3.7647110220e-1,    3.9919844270e-1,    4.2192578320e-1,
        4.4331622120e-1,    4.6604356170e-1,    4.8877090220e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
44,1,2,23
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n44_grid23() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(44, FilterType::MultipleBand, &bands, 23);

    let expected_impulse_response = vec![
        2.5533165900e-4,    2.3662275630e-4,   -3.0104361940e-4,   -1.1923700800e-3,
        -1.3673594220e-3,    2.7509441130e-4,    3.1932266430e-3,    4.5453449710e-3,
        1.2680948710e-3,   -6.0281790790e-3,   -1.1284573000e-2,   -6.9135362280e-3,
        8.0697182570e-3,    2.3028098050e-2,    2.1168159320e-2,   -5.3688990880e-3,
        -4.2270205910e-2,   -5.5646620690e-2,   -1.3709884140e-2,    8.6139515040e-2,
        2.0672865210e-1,    2.8908669950e-1,
    ];
    let expected_deviations = vec![
        1.7632094390e-4,    1.7632094390e-4,
    ];
    let expected_deviation_dbs = vec![
        1.5312789470e-3,    -7.5073921200e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.4703562260e-2,    4.7430809590e-2,    6.8181768060e-2,
        8.3992019300e-2,    9.5849707720e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.0395256580e-1,    2.1482211350e-1,    2.3063236470e-1,    2.4940703810e-1,
        2.6917013530e-1,    2.8992140290e-1,    3.1067267060e-1,    3.3241209390e-1,
        3.5415151720e-1,    3.7687909600e-1,    3.9861851930e-1,    4.2134609820e-1,
        4.4407367710e-1,    4.6581310030e-1,    4.8854067920e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
44,1,2,24
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n44_grid24() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(44, FilterType::MultipleBand, &bands, 24);

    let expected_impulse_response = vec![
        2.5526719400e-4,    2.3682278700e-4,   -3.0059504210e-4,   -1.1920725230e-3,
        -1.3678068060e-3,    2.7390313330e-4,    3.1921758780e-3,    4.5458069070e-3,
        1.2702895330e-3,   -6.0258228330e-3,   -1.1284487320e-2,   -6.9166123870e-3,
        8.0655748020e-3,    2.3026695470e-2,    2.1171566100e-2,   -5.3630736660e-3,
        -4.2266920210e-2,   -5.5649347600e-2,   -1.3716723770e-2,    8.6134292190e-2,
        2.0672965050e-1,    2.8909331560e-1,
    ];
    let expected_deviations = vec![
        1.7621229930e-4,    1.7621229930e-4,
    ];
    let expected_deviation_dbs = vec![
        1.5302435490e-3,    -7.5079277040e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.5568172340e-2,    4.8295423390e-2,    6.8181768060e-2,
        8.4280237560e-2,    9.5643863080e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.0378787820e-1,    2.1515150370e-1,    2.3030300440e-1,    2.4924238030e-1,
        2.6912873980e-1,    2.8996205330e-1,    3.1079536680e-1,    3.3257564900e-1,
        3.5435593130e-1,    3.7708318230e-1,    3.9886346460e-1,    4.2159071560e-1,
        4.4337099790e-1,    4.6609824900e-1,    4.8882550000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
107,1,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n107_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(107, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -3.2031977070e-8,   -1.1279968960e-7,   -1.3722747380e-7,    1.5442859080e-7,
        8.6048396500e-7,    1.3387629000e-6,    1.8289131280e-7,   -3.2945740710e-6,
        -6.6891975620e-6,   -4.3344521150e-6,    7.5733582890e-6,    2.2738096960e-5,
        2.3006497940e-5,   -7.4966592360e-6,   -5.7651617680e-5,   -7.9653385910e-5,
        -2.0970919650e-5,    1.1090256880e-4,    2.1092007230e-4,    1.3523641970e-4,
        -1.4966847080e-4,   -4.5186682840e-4,   -4.3602063670e-4,    7.3462106230e-5,
        7.9582404580e-4,    1.0542634410e-3,    3.1941043560e-4,   -1.1275778520e-3,
        -2.0963796410e-3,   -1.3376544230e-3,    1.1457334040e-3,    3.5407647960e-3,
        3.3469593620e-3,   -3.0626426450e-4,   -5.1006805150e-3,   -6.6527808090e-3,
        -2.1825232540e-3,    6.0860579830e-3,    1.1339232330e-2,    7.3240338820e-3,
        -5.2762962880e-3,   -1.7125137150e-2,   -1.6382556410e-2,    6.8027159430e-4,
        2.3315064610e-2,    3.1506434080e-2,    1.1545025740e-2,   -2.8900763020e-2,
        -5.9447083620e-2,   -4.4086005540e-2,    3.2808169720e-2,    1.4998741450e-1,
        2.5629872080e-1,    2.9911974070e-1,
    ];
    let expected_deviations = vec![
        3.0463320750e-9,    3.0463320750e-9,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.7032446290e2,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    4.0509263050e-3,    1.5624996270e-2,    3.0671283600e-2,
        3.8194447760e-2,    4.9768552180e-2,    5.7870425280e-2,    6.5972276030e-2,
        7.2337992490e-2,    8.1018514930e-2,    8.5648126900e-2,    9.0856440370e-2,
        9.4328649340e-2,    9.7800858320e-2,    9.8958261310e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0057870450e-1,    2.0289351050e-1,    2.0636571940e-1,
        2.1099533140e-1,    2.1678234640e-1,    2.2256936130e-1,    2.2951377930e-1,
        2.3645819720e-1,    2.4398131670e-1,    2.5208315250e-1,    2.6018497350e-1,
        2.6828679440e-1,    2.7638861540e-1,    2.8506913780e-1,    2.9374966030e-1,
        3.0243018270e-1,    3.1111070510e-1,    3.1979122760e-1,    3.2847175000e-1,
        3.3715227250e-1,    3.4641149640e-1,    3.5509201880e-1,    3.6435124280e-1,
        3.7303176520e-1,    3.8229098920e-1,    3.9097151160e-1,    4.0023073550e-1,
        4.0891125800e-1,    4.1817048190e-1,    4.2742970590e-1,    4.3611022830e-1,
        4.4536945220e-1,    4.5462867620e-1,    4.6388790010e-1,    4.7256842260e-1,
        4.8182764650e-1,    4.9108687040e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
107,1,2,17
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n107_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(107, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        -2.7173491900e-8,   -9.0176612840e-8,   -8.8603336220e-8,    1.9940178220e-7,
        8.1341181610e-7,    1.1107907770e-6,   -1.6014563190e-7,   -3.4135082390e-6,
        -6.1347054730e-6,   -3.0745545700e-6,    8.6508180170e-6,    2.2098487530e-5,
        1.9874334610e-5,   -1.1488582460e-5,   -5.8525882200e-5,   -7.4114970630e-5,
        -1.0798226870e-5,    1.1750611880e-4,    2.0459800730e-4,    1.1526982420e-4,
        -1.6968365530e-4,   -4.5113055970e-4,   -4.0515675210e-4,    1.1702883780e-4,
        8.1389938710e-4,    1.0183502450e-3,    2.4370430040e-4,   -1.1840726950e-3,
        -2.0724283530e-3,   -1.2295541820e-3,    1.2617823670e-3,    3.5580068360e-3,
        3.2218983400e-3,   -4.9505662170e-4,   -5.1954542290e-3,   -6.5455227160e-3,
        -1.9268698520e-3,    6.2900749040e-3,    1.1299810370e-2,    7.0339958180e-3,
        -5.6014116850e-3,   -1.7205763610e-2,   -1.6115676610e-2,    1.1061740810e-3,
        2.3551192130e-2,    3.1332273040e-2,    1.1073134840e-2,   -2.9293831440e-2,
        -5.9425845740e-2,   -4.3646294620e-2,    3.3318281170e-2,    1.5014818310e-1,
        2.5597062710e-1,    2.9856628180e-1,
    ];
    let expected_deviations = vec![
        2.7447721870e-9,    2.7447721870e-9,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.7122987370e2,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0348583570e-2,    1.7973855140e-2,    2.5054456670e-2,
        4.1939005260e-2,    5.3376939150e-2,    5.7734247300e-2,    6.6448837520e-2,
        7.1350775660e-2,    8.0065332350e-2,    8.6601249870e-2,    9.2047847810e-2,
        9.5315806570e-2,    9.8039105530e-2,    9.9128425120e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0054467020e-1,    2.0272333920e-1,    2.0653600990e-1,
        2.1089334790e-1,    2.1634002030e-1,    2.2287602720e-1,    2.2941203420e-1,
        2.3649270830e-1,    2.4411804970e-1,    2.5174337630e-1,    2.5991338490e-1,
        2.6808339360e-1,    2.7625340220e-1,    2.8496807810e-1,    2.9368275400e-1,
        3.0239742990e-1,    3.1111210580e-1,    3.1982678170e-1,    3.2854145770e-1,
        3.3725613360e-1,    3.4597080950e-1,    3.5523015260e-1,    3.6394482850e-1,
        3.7320417170e-1,    3.8191884760e-1,    3.9117819070e-1,    3.9989286660e-1,
        4.0915220980e-1,    4.1841155290e-1,    4.2712622880e-1,    4.3638557200e-1,
        4.4564491510e-1,    4.5435959100e-1,    4.6361893420e-1,    4.7287827730e-1,
        4.8159295320e-1,    4.9085229640e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
107,1,2,23
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n107_grid23() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(107, FilterType::MultipleBand, &bands, 23);

    let expected_impulse_response = vec![
        -3.0573296780e-8,   -1.0587321240e-7,   -1.2250437460e-7,    1.6669288530e-7,
        8.4084479110e-7,    1.2590573990e-6,    6.9628704580e-8,   -3.3209094000e-6,
        -6.4751952780e-6,   -3.8826751730e-6,    7.9299607020e-6,    2.2447611630e-5,
        2.1819347240e-5,   -8.9326285890e-6,   -5.7847490100e-5,   -7.7437434810e-5,
        -1.7124099030e-5,    1.1320393970e-4,    2.0818441410e-4,    1.2740060630e-4,
        -1.5717267520e-4,   -4.5102127480e-4,   -4.2350424340e-4,    9.0436835310e-5,
        8.0216192870e-4,    1.0391341060e-3,    2.8910514080e-4,   -1.1492491470e-3,
        -2.0853737370e-3,   -1.2934322000e-3,    1.1918306120e-3,    3.5459289790e-3,
        3.2947412690e-3,   -3.8286013300e-4,   -5.1374132740e-3,   -6.6067865120e-3,
        -2.0772307180e-3,    6.1682006340e-3,    1.1320500630e-2,    7.2031854650e-3,
        -5.4095271040e-3,   -1.7156042160e-2,   -1.6270227730e-2,    8.5662526540e-4,
        2.3411171510e-2,    3.1432207670e-2,    1.1348421690e-2,   -2.9063213620e-2,
        -5.9436839070e-2,   -4.3902289120e-2,    3.3020328730e-2,    1.5005381410e-1,
        2.5616177920e-1,    2.9888913040e-1,
    ];
    let expected_deviations = vec![
        2.9660449650e-9,    2.9660449650e-9,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.7055645750e2,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    4.4283410530e-3,    1.8115941440e-2,    3.0998412520e-2,
        4.0660265830e-2,    5.0322119150e-2,    5.5958200250e-2,    6.5620049830e-2,
        7.2061285380e-2,    8.0917984250e-2,    8.5346333680e-2,    9.1384992000e-2,
        9.4605609770e-2,    9.8228804770e-2,    9.9436536430e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0080514250e-1,    2.0281799140e-1,    2.0644111930e-1,
        2.1086938680e-1,    2.1650536360e-1,    2.2254391010e-1,    2.2938759620e-1,
        2.3663385210e-1,    2.4428267780e-1,    2.5193151830e-1,    2.5998291370e-1,
        2.6803430910e-1,    2.7648827430e-1,    2.8494223950e-1,    2.9339620470e-1,
        3.0225273970e-1,    3.1070670490e-1,    3.1956323980e-1,    3.2841977480e-1,
        3.3727630970e-1,    3.4613284470e-1,    3.5498937960e-1,    3.6424848440e-1,
        3.7310501930e-1,    3.8196155430e-1,    3.9122065900e-1,    4.0007719400e-1,
        4.0933629870e-1,    4.1819283370e-1,    4.2745193840e-1,    4.3630847330e-1,
        4.4556757810e-1,    4.5442411300e-1,    4.6368321780e-1,    4.7253975270e-1,
        4.8179885750e-1,    4.9105796220e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
107,1,2,24
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n107_grid24() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(107, FilterType::MultipleBand, &bands, 24);

    let expected_impulse_response = vec![
        -2.4267876600e-8,   -7.7271096190e-8,   -6.5018646470e-8,    2.0615856040e-7,
        7.4214256070e-7,    9.2328082250e-7,   -3.6418700230e-7,   -3.3449805410e-6,
        -5.5077148320e-6,   -2.0464169670e-6,    9.1797528510e-6,    2.0924253480e-5,
        1.6795334890e-5,   -1.4473167540e-5,   -5.7815730540e-5,   -6.7618246250e-5,
        -1.5822636210e-6,    1.2111664550e-4,    1.9503336810e-4,    9.4712173450e-5,
        -1.8589297540e-4,   -4.4377322770e-4,   -3.6967598130e-4,    1.5833329230e-4,
        8.2234409640e-4,    9.7132480000e-4,    1.6452290580e-4,   -1.2311202010e-3,
        -2.0304333880e-3,   -1.1076866650e-3,    1.3747323540e-3,    3.5538068040e-3,
        3.0705546960e-3,   -6.9453241300e-4,   -5.2735456270e-3,   -6.4025586470e-3,
        -1.6417590670e-3,    6.4928587530e-3,    1.1224819350e-2,    6.6972142090e-3,
        -5.9493514710e-3,   -1.7264721920e-2,   -1.5794560310e-2,    1.5803889840e-3,
        2.3792231460e-2,    3.1112454830e-2,    1.0536068120e-2,   -2.9723659160e-2,
        -5.9383913870e-2,   -4.3141439560e-2,    3.3891070630e-2,    1.5032270550e-1,
        2.5559639930e-1,    2.9794016480e-1,
    ];
    let expected_deviations = vec![
        2.5517237210e-9,    2.5517237210e-9,
    ];
    let expected_deviation_dbs = vec![
        0.0000000000e0,    -1.7186332700e2,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    3.0864197760e-3,    1.5817895530e-2,    2.7006160470e-2,
        3.6265414210e-2,    4.3981458990e-2,    5.5941328410e-2,    6.5200604500e-2,
        7.0987693970e-2,    7.9861231150e-2,    8.6034126580e-2,    9.1435410080e-2,
        9.5293469730e-2,    9.8765723410e-2,    9.9537335340e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0077161490e-1,    2.0308645070e-1,    2.0655870440e-1,
        2.1080257000e-1,    2.1658965950e-1,    2.2276255490e-1,    2.2932125630e-1,
        2.3665156960e-1,    2.4398188290e-1,    2.5169792770e-1,    2.5979954000e-1,
        2.6790115240e-1,    2.7638855580e-1,    2.8487595920e-1,    2.9336336260e-1,
        3.0223655700e-1,    3.1072396040e-1,    3.1959715490e-1,    3.2847034930e-1,
        3.3734354380e-1,    3.4621673820e-1,    3.5508993270e-1,    3.6396312710e-1,
        3.7283632160e-1,    3.8209530710e-1,    3.9096850160e-1,    4.0022748710e-1,
        4.0910068150e-1,    4.1797387600e-1,    4.2723286150e-1,    4.3610605600e-1,
        4.4536504150e-1,    4.5462402700e-1,    4.6349722150e-1,    4.7275620700e-1,
        4.8162940140e-1,    4.9088838700e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
108,1,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n108_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(108, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        1.2131450830e-6,    3.4488512030e-6,    2.6248910670e-6,   -8.7666912800e-6,
        -3.0772767790e-5,   -3.6031244240e-5,    2.2332194930e-5,    1.4514214130e-4,
        2.1108017240e-4,    2.3236978450e-5,   -4.4801196780e-4,   -8.1318808950e-4,
        -4.1048805000e-4,    9.7215839200e-4,    2.3533469070e-3,    1.8746622370e-3,
        -1.3825503410e-3,   -5.4257260640e-3,   -5.7904524730e-3,    4.9298838710e-4,
        1.0206883770e-2,    1.4036823060e-2,    4.1821529160e-3,   -1.5591505910e-2,
        -2.8252664950e-2,   -1.6398750250e-2,    1.8293639640e-2,    4.8444394020e-2,
        4.0211372080e-2,   -1.2580775660e-2,   -7.1403175590e-2,   -7.7952235940e-2,
        -8.6461585020e-3,    8.9886516330e-2,    1.2775817510e-1,    5.1200002430e-2,
        -9.3607194720e-2,   -1.8194213510e-1,   -1.1610157040e-1,    7.2461962700e-2,
        2.2767610850e-1,    1.9687679410e-1,   -2.1264929320e-2,   -2.5080230830e-1,
        -2.7961912750e-1,   -5.5776122960e-2,    2.4259865280e-1,    3.4733879570e-1,
        1.4211297040e-1,   -2.1124321220e-1,   -3.9215397830e-1,   -1.9808754320e-1,
        2.6282659170e-1,    6.4355528350e-1,
    ];
    let expected_deviations = vec![
        6.0889064460e-8,    6.0889064460e-8,
    ];
    let expected_deviation_dbs = vec![
        1.0354386860e-6,    -1.4430921940e2,
    ];
    let expected_extremal_frequencies = vec![
        1.1574074160e-3,    4.6296301300e-3,    1.1574072760e-2,    1.2152776120e-2,
        1.7939809710e-2,    2.6041656730e-2,    3.2407395540e-2,    3.5879626870e-2,
        4.3402794750e-2,    5.0347257410e-2,    6.1342656610e-2,    7.1759290990e-2,
        7.8703708950e-2,    8.5069425400e-2,    9.1435141860e-2,    9.6643455330e-2,
        9.8958261310e-2,    1.0000000150e-1,    2.0115740600e-1,    2.0578701790e-1,
        2.1273143590e-1,    2.2141195830e-1,    2.2835637630e-1,    2.3356468980e-1,
        2.3935170470e-1,    2.4687482420e-1,    2.5497666000e-1,    2.6307848100e-1,
        2.7118030190e-1,    2.7928212290e-1,    2.8738394380e-1,    2.9606446620e-1,
        3.0358758570e-1,    3.1226810810e-1,    3.2094863060e-1,    3.2905045150e-1,
        3.3773097400e-1,    3.4641149640e-1,    3.5509201880e-1,    3.6377254130e-1,
        3.7245306370e-1,    3.8113358620e-1,    3.8981410860e-1,    3.9849463110e-1,
        4.0717515350e-1,    4.1585567590e-1,    4.2453619840e-1,    4.3379542230e-1,
        4.4247594480e-1,    4.5115646720e-1,    4.6041569110e-1,    4.6909621360e-1,
        4.7777673600e-1,    4.8645725850e-1,    4.9571648240e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
108,1,2,17
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n108_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(108, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        -2.1217538230e-6,   -3.5501507230e-6,   -2.0200943710e-6,    2.2454366900e-6,
        1.4981753340e-5,    3.0794446500e-5,    1.9309192790e-5,   -4.9980586480e-5,
        -1.4214047410e-4,   -1.2890508510e-4,    1.0299443970e-4,    4.4491622250e-4,
        5.0356582510e-4,   -7.9466219180e-5,   -1.0725072350e-3,   -1.4846712580e-3,
        -3.2672588710e-4,    2.0692055110e-3,    3.5674790850e-3,    1.7681899480e-3,
        -3.1699361280e-3,   -7.2417412880e-3,   -5.2880356090e-3,    3.5280513110e-3,
        1.2638092040e-2,    1.2135788800e-2,   -1.5442078000e-3,   -1.9044533370e-2,
        -2.3266706620e-2,   -4.9931039100e-3,    2.4502560500e-2,    3.8587041200e-2,
        1.8301256000e-2,   -2.5784347210e-2,   -5.6211024520e-2,   -3.9575692270e-2,
        1.8988307570e-2,    7.2126351300e-2,    6.7855551840e-2,   -7.3070824150e-4,
        -8.0581389370e-2,   -9.9276423450e-2,   -3.0476380140e-2,    7.5110085310e-2,
        1.2708580490e-1,    7.3482915760e-2,   -4.9198240040e-2,   -1.4203807710e-1,
        -1.2603485580e-1,   -8.1041306260e-3,    1.2852610650e-1,    1.9954767820e-1,
        1.8992951510e-1,    1.5498283510e-1,
    ];
    let expected_deviations = vec![
        6.3787389590e-8,    6.3787389590e-8,
    ];
    let expected_deviation_dbs = vec![
        1.0354386860e-6,    -1.4390530400e2,
    ];
    let expected_extremal_frequencies = vec![
        5.4466229630e-4,    5.9912847360e-3,    9.2592583970e-3,    9.8039209840e-3,
        1.6339870170e-2,    1.8518516790e-2,    2.7233103290e-2,    3.4858379510e-2,
        4.7930303960e-2,    5.6644920260e-2,    6.6448837520e-2,    7.5163394210e-2,
        8.2788631320e-2,    8.8779889050e-2,    9.3681827190e-2,    9.7494445740e-2,
        9.9128425120e-2,    1.0000000150e-1,    2.0054467020e-1,    2.0326800640e-1,
        2.0708067720e-1,    2.1252734960e-1,    2.1906335650e-1,    2.2668869790e-1,
        2.3431403930e-1,    2.4248404800e-1,    2.5065404180e-1,    2.5936871770e-1,
        2.6753872630e-1,    2.7570873500e-1,    2.8442341090e-1,    2.9259341960e-1,
        3.0130809550e-1,    3.0947810410e-1,    3.1819278000e-1,    3.2690745590e-1,
        3.3507746460e-1,    3.4379214050e-1,    3.5250681640e-1,    3.6122149230e-1,
        3.6993616820e-1,    3.7865084410e-1,    3.8791018720e-1,    3.9662486310e-1,
        4.0533953910e-1,    4.1405421500e-1,    4.2331355810e-1,    4.3202823400e-1,
        4.4128757720e-1,    4.5054692030e-1,    4.5926159620e-1,    4.6852093940e-1,
        4.7778028250e-1,    4.8758429290e-1,    4.9629896880e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
108,1,2,23
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n108_grid23() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(108, FilterType::MultipleBand, &bands, 23);

    let expected_impulse_response = vec![
        -3.8333490690e-7,   -1.1039578570e-6,   -9.5967266130e-7,    2.0147017490e-6,
        7.3114383670e-6,    8.6512336570e-6,   -2.8402430420e-6,   -2.6225507100e-5,
        -3.9740887590e-5,   -1.1341014210e-5,    6.3258012230e-5,    1.2647725820e-4,
        8.2926671890e-5,   -1.0260672570e-4,   -3.0934874670e-4,   -2.9790040570e-4,
        7.2514347270e-5,    6.0351123100e-4,    7.8061630480e-4,    1.9190416790e-4,
        -9.3434297010e-4,   -1.6439694450e-3,   -9.5713086190e-4,    1.0597066720e-3,
        2.8829537330e-3,    2.5310793430e-3,   -5.2433443490e-4,   -4.2308443230e-3,
        -5.1076039670e-3,   -1.2901725710e-3,    5.0410879780e-3,    8.5346568380e-3,
        4.9746441650e-3,   -4.2824032720e-3,   -1.2082571160e-2,   -1.0775052010e-2,
        7.0927618070e-4,    1.4316802840e-2,    1.8235299740e-2,    6.8134129980e-3,
        -1.3132718390e-2,   -2.5922430680e-2,   -1.8990049140e-2,    5.8185588570e-3,
        3.1260937450e-2,    3.6154784260e-2,    1.1594252660e-2,   -2.9972158370e-2,
        -5.9903915970e-2,   -4.9879781900e-2,    1.0714983570e-2,    1.0705752670e-1,
        2.0482607190e-1,    2.6595675950e-1,
    ];
    let expected_deviations = vec![
        6.2474207140e-8,    6.2474207140e-8,
    ];
    let expected_deviation_dbs = vec![
        1.0354386860e-6,    -1.4408598330e2,
    ];
    let expected_extremal_frequencies = vec![
        8.0515276640e-3,    9.2592565340e-3,    1.1674714270e-2,    1.4492748310e-2,
        1.8518518660e-2,    2.6570063080e-2,    3.3816453070e-2,    3.6634493620e-2,
        4.5491192490e-2,    6.1996858570e-2,    7.4476748700e-2,    8.1723138690e-2,
        8.8164374230e-2,    9.2592723670e-2,    9.7021073100e-2,    9.8631381990e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0080514250e-1,    2.0362313090e-1,
        2.0724625890e-1,    2.1247966590e-1,    2.1851821240e-1,    2.2536189850e-1,
        2.3220558460e-1,    2.3985441030e-1,    2.4790580570e-1,    2.5595721600e-1,
        2.6400861140e-1,    2.7246257660e-1,    2.8091654180e-1,    2.8937050700e-1,
        2.9822704200e-1,    3.0668100710e-1,    3.1553754210e-1,    3.2439407710e-1,
        3.3325061200e-1,    3.4210714700e-1,    3.5096368190e-1,    3.5982021690e-1,
        3.6907932160e-1,    3.7793585660e-1,    3.8679239150e-1,    3.9605149630e-1,
        4.0490803120e-1,    4.1376456620e-1,    4.2302367090e-1,    4.3188020590e-1,
        4.4113931060e-1,    4.4999584560e-1,    4.5925495030e-1,    4.6811148520e-1,
        4.7737059000e-1,    4.8622712490e-1,    4.9548622970e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
108,1,2,24
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n108_grid24() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(108, FilterType::MultipleBand, &bands, 24);

    let expected_impulse_response = vec![
        -4.6491754800e-7,   -1.4354424140e-6,   -1.3834375070e-6,    2.4978853620e-6,
        9.9464750750e-6,    1.2591852280e-5,   -2.5850908970e-6,   -3.6087734770e-5,
        -5.8556586740e-5,   -2.2560721850e-5,    8.4132079790e-5,    1.8480274590e-4,
        1.3908457190e-4,   -1.2027187770e-4,   -4.3852912490e-4,   -4.6959135220e-4,
        2.5300876590e-5,    8.1229908390e-4,    1.1678732700e-3,    4.4007605170e-4,
        -1.1576483960e-3,   -2.3311208930e-3,   -1.6116744370e-3,    1.1099434920e-3,
        3.8576759400e-3,    3.8027903070e-3,   -9.0271933000e-5,   -5.3054573950e-3,
        -7.0955930280e-3,   -2.5610830630e-3,    5.8382125570e-3,    1.1109545830e-2,
        7.3201344350e-3,   -4.3317186650e-3,   -1.4852008780e-2,   -1.4159057290e-2,
        -3.5873800520e-4,    1.6732161860e-2,    2.2290747610e-2,    9.0683652090e-3,
        -1.4715532770e-2,   -3.0064444990e-2,   -2.2122796620e-2,    6.3946554440e-3,
        3.4940220420e-2,    3.9594508710e-2,    1.1796547100e-2,   -3.2938148830e-2,
        -6.3092261550e-2,   -5.0348330290e-2,    1.3129522090e-2,    1.0973242670e-1,
        2.0505243540e-1,    2.6363885400e-1,
    ];
    let expected_deviations = vec![
        6.5921959450e-8,    6.5921959450e-8,
    ];
    let expected_deviation_dbs = vec![
        1.0354386860e-6,    -1.4361940000e2,
    ];
    let expected_extremal_frequencies = vec![
        3.4722222480e-3,    9.2592565340e-3,    9.6450587730e-3,    1.3503081170e-2,
        1.9290115680e-2,    2.7777764950e-2,    3.6265414210e-2,    4.3595656750e-2,
        5.2083306010e-2,    6.4428992570e-2,    6.7901246250e-2,    7.1759305890e-2,
        7.2916723790e-2,    8.8734768330e-2,    9.5293469730e-2,    9.7994111480e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0077161490e-1,    2.0424386860e-1,
        2.0964515210e-1,    2.1620385350e-1,    2.2199094300e-1,    2.2816383840e-1,
        2.3472253980e-1,    2.4166704710e-1,    2.4938316640e-1,    2.5709900260e-1,
        2.6520061490e-1,    2.7368801830e-1,    2.8217542170e-1,    2.9066282510e-1,
        2.9953601960e-1,    3.0802342300e-1,    3.1689661740e-1,    3.2538402080e-1,
        3.3425721530e-1,    3.4313040970e-1,    3.5200360420e-1,    3.6087679860e-1,
        3.6974999310e-1,    3.7862318750e-1,    3.8749638200e-1,    3.9636957650e-1,
        4.0562856200e-1,    4.1450175640e-1,    4.2337495090e-1,    4.3224814530e-1,
        4.4150713090e-1,    4.5038032530e-1,    4.5925351980e-1,    4.6851250530e-1,
        4.7738569970e-1,    4.8664468530e-1,    4.9551787970e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
128,1,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n128_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -3.9431947130e-7,   -1.5167981930e-6,   -1.4715758430e-6,    3.8804014370e-6,
        1.3791733180e-5,    1.5020612410e-5,   -1.2044572940e-5,   -6.4237407060e-5,
        -8.6882479080e-5,    6.6485517890e-7,    2.0175839020e-4,    3.4413492540e-4,
        1.5132541010e-4,   -4.5178312580e-4,   -1.0322452290e-3,   -7.8930612650e-4,
        6.6925317510e-4,    2.4675354360e-3,    2.6174932720e-3,   -2.6298721790e-4,
        -4.7965133560e-3,   -6.7037525590e-3,   -2.1349946040e-3,    7.4986713010e-3,
        1.4160348100e-2,    8.8345967230e-3,   -8.7264403700e-3,   -2.5366447870e-2,
        -2.2764455530e-2,    4.8079411500e-3,    3.8829311730e-2,    4.6409003440e-2,
        9.5912562680e-3,   -5.0125356760e-2,   -8.0070748930e-2,   -4.0219023820e-2,
        5.1633726810e-2,    1.2000815570e-1,    9.0878054500e-2,   -3.3734809610e-2,
        -1.5740796920e-1,   -1.6058363020e-1,   -1.2339051810e-2,    1.7917639020e-1,
        2.4144816400e-1,    9.0633802120e-2,   -1.7099983990e-1,   -3.1864511970e-1,
        -1.9700932500e-1,    1.2215401980e-1,    3.7328621750e-1,    3.1781637670e-1,
        -3.0565485360e-2,   -3.8803181050e-1,   -4.3203431370e-1,   -9.3724809590e-2,
        3.5450008510e-1,    5.1741307970e-1,    2.2734135390e-1,   -2.8415262700e-1,
        -5.6348234420e-1,   -3.2434508200e-1,    2.8810524940e-1,    7.9965615270e-1,
    ];
    let expected_deviations = vec![
        6.3256088370e-8,    6.3256088370e-8,
    ];
    let expected_deviation_dbs = vec![
        1.0354386860e-6,    -1.4397795100e2,
    ];
    let expected_extremal_frequencies = vec![
        5.3710937500e-3,    1.1230468750e-2,    1.6113281250e-2,    1.7578125000e-2,
        2.0019531250e-2,    2.6855468750e-2,    3.0761718750e-2,    3.6621093750e-2,
        4.2480468750e-2,    4.9316406250e-2,    5.4687500000e-2,    6.0058593750e-2,
        6.6406250000e-2,    7.3730468750e-2,    7.7636718750e-2,    8.1054687500e-2,
        9.2285156250e-2,    9.6679687500e-2,    9.9121093750e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0048828420e-1,    2.0244140920e-1,    2.0488281550e-1,
        2.0878906550e-1,    2.1367187800e-1,    2.1904297170e-1,    2.2441406550e-1,
        2.3076172170e-1,    2.3662109670e-1,    2.4345703420e-1,    2.4980469050e-1,
        2.5664061310e-1,    2.6396483180e-1,    2.7080076930e-1,    2.7812498810e-1,
        2.8496092560e-1,    2.9228514430e-1,    2.9960936310e-1,    3.0693358180e-1,
        3.1425780060e-1,    3.2158201930e-1,    3.2890623810e-1,    3.3671873810e-1,
        3.4404295680e-1,    3.5185545680e-1,    3.5917967560e-1,    3.6650389430e-1,
        3.7431639430e-1,    3.8164061310e-1,    3.8945311310e-1,    3.9677733180e-1,
        4.0458983180e-1,    4.1240233180e-1,    4.1972655060e-1,    4.2753905060e-1,
        4.3486326930e-1,    4.4267576930e-1,    4.5048826930e-1,    4.5781248810e-1,
        4.6562498810e-1,    4.7343748810e-1,    4.8076170680e-1,    4.8857420680e-1,
        4.9638670680e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
128,1,2,17
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn multipleband_lowpass_n128_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        5.6462539530e-7,    2.1475527770e-6,    2.7946657610e-6,   -2.5387512320e-6,
        -1.5399080440e-5,   -2.3636843250e-5,   -2.5791778170e-6,    5.7177723650e-5,
        1.1053364870e-4,    6.2760496800e-5,   -1.3656046940e-4,   -3.6417599770e-4,
        -3.2490363810e-4,    1.8876412650e-4,    9.2226377460e-4,    1.1049972380e-3,
        5.6835066060e-5,   -1.8526244680e-3,   -2.9049883600e-3,   -1.2477622370e-3,
        2.9146266170e-3,    6.2768668870e-3,    4.5156525450e-3,   -3.2208755150e-3,
        -1.1452095580e-2,   -1.1352542790e-2,    9.5394207160e-4,    1.7766844480e-2,
        2.3113876580e-2,    6.6106552260e-3,   -2.3088734600e-2,   -4.0155701340e-2,
        -2.2542573510e-2,    2.3598847910e-2,    6.0856733470e-2,    4.9113824960e-2,
        -1.4280043540e-2,   -8.0990403890e-2,   -8.6353845890e-2,   -9.7688529640e-3,
        9.3930371110e-2,    1.3086265330e-1,    5.1460355520e-2,   -9.1922447090e-2,
        -1.7550410330e-1,   -1.0991543530e-1,    6.8204298620e-2,    2.1034273510e-1,
        1.7947687210e-1,   -1.9303772600e-2,   -2.2466954590e-1,   -2.5003618000e-1,
        -5.3373012690e-2,    2.0936809480e-1,    3.0860435960e-1,    1.4352782070e-1,
        -1.5812967720e-1,   -3.4110319610e-1,   -2.4276858570e-1,    6.2761932610e-2,
        3.2927161460e-1,    3.5661143060e-1,    1.6536368430e-1,   -3.1210182230e-2,
    ];
    let expected_deviations = vec![
        6.4412283510e-8,    6.4412283510e-8,
    ];
    let expected_deviation_dbs = vec![
        1.0354386860e-6,    -1.4382063290e2,
    ];
    let expected_extremal_frequencies = vec![
        3.2169118060e-3,    8.2720592620e-3,    1.1029414830e-2,    1.2408092620e-2,
        1.8382363020e-2,    2.2058837120e-2,    2.5735311210e-2,    3.2628700140e-2,
        5.0091952090e-2,    5.5606663230e-2,    6.1580933630e-2,    6.6636085510e-2,
        7.4448592960e-2,    8.0422863360e-2,    8.4558896720e-2,    9.2371404170e-2,
        9.5128759740e-2,    9.8345674570e-2,    9.9264793100e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0045955480e-1,    2.0229776200e-1,    2.0551462470e-1,
        2.0919103920e-1,    2.1378655730e-1,    2.1884162720e-1,    2.2435624900e-1,
        2.3033042250e-1,    2.3676414790e-1,    2.4319787320e-1,    2.5009116530e-1,
        2.5652489070e-1,    2.6387771960e-1,    2.7077099680e-1,    2.7766427400e-1,
        2.8501710300e-1,    2.9236993190e-1,    2.9972276090e-1,    3.0707558990e-1,
        3.1442841890e-1,    3.2178124790e-1,    3.2913407680e-1,    3.3648690580e-1,
        3.4429928660e-1,    3.5165211560e-1,    3.5900494460e-1,    3.6681732540e-1,
        3.7417015430e-1,    3.8198253510e-1,    3.8933536410e-1,    3.9714774490e-1,
        4.0450057390e-1,    4.1231295470e-1,    4.1966578360e-1,    4.2747816440e-1,
        4.3529054520e-1,    4.4264337420e-1,    4.5045575500e-1,    4.5780858400e-1,
        4.6562096480e-1,    4.7343334560e-1,    4.8078617450e-1,    4.8859855530e-1,
        4.9595138430e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
7,1,3,16
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n7_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(7, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        4.6523381020e-2,   -2.0000003280e-1,   -4.6523325150e-2,    2.0978866520e-1,
    ];
    let expected_deviations = vec![
        1.9021129610e-1,    5.7063388820e-1,    1.9021129610e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.4415273670e1,     3.9214992520e0,    -1.4415273670e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5000000600e-1,    2.0000000300e-1,    4.0000000600e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
7,1,3,17
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n7_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(7, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        4.6523381020e-2,   -2.0000003280e-1,   -4.6523325150e-2,    2.0978866520e-1,
    ];
    let expected_deviations = vec![
        1.9021129610e-1,    5.7063388820e-1,    1.9021129610e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.4415273670e1,     3.9214992520e0,    -1.4415273670e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5000000600e-1,    2.0000000300e-1,    4.0000000600e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
7,1,3,31
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n7_grid31() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(7, FilterType::MultipleBand, &bands, 31);

    let expected_impulse_response = vec![
        4.6523381020e-2,   -2.0000003280e-1,   -4.6523325150e-2,    2.0978866520e-1,
    ];
    let expected_deviations = vec![
        1.9021129610e-1,    5.7063388820e-1,    1.9021129610e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.4415273670e1,     3.9214992520e0,    -1.4415273670e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5000000600e-1,    2.0000000300e-1,    4.0000000600e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
7,1,3,32
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n7_grid32() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(7, FilterType::MultipleBand, &bands, 32);

    let expected_impulse_response = vec![
        4.6523381020e-2,   -2.0000003280e-1,   -4.6523325150e-2,    2.0978866520e-1,
    ];
    let expected_deviations = vec![
        1.9021129610e-1,    5.7063388820e-1,    1.9021129610e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.4415273670e1,     3.9214992520e0,    -1.4415273670e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5000000600e-1,    2.0000000300e-1,    4.0000000600e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
8,1,3,16
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n8_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(8, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        1.6926960650e-1,   -1.0437540710e-1,   -2.8724658490e-1,    2.2231715920e-1,
    ];
    let expected_deviations = vec![
        1.1953954400e-1,    3.5861861710e-1,    1.1953954400e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.8449768070e1,     2.6619510650e0,    -1.8449768070e1,
    ];
    let expected_extremal_frequencies = vec![
        8.5937500000e-2,    1.5000000600e-1,    2.0000000300e-1,    4.0000000600e-1,
        4.5468750600e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
8,1,3,17
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n8_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(8, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        1.6940073670e-1,   -1.0426492990e-1,   -2.8715848920e-1,    2.2253540160e-1,
    ];
    let expected_deviations = vec![
        1.1954069140e-1,    3.5862207410e-1,    1.1954069140e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.8449684140e1,     2.6619734760e0,    -1.8449684140e1,
    ];
    let expected_extremal_frequencies = vec![
        8.8235296310e-2,    1.5000000600e-1,    2.0000000300e-1,    4.0000000600e-1,
        4.5147064330e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
8,1,3,31
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n8_grid31() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(8, FilterType::MultipleBand, &bands, 31);

    let expected_impulse_response = vec![
        1.6928377750e-1,   -1.0424983500e-1,   -2.8709349040e-1,    2.2243610020e-1,
    ];
    let expected_deviations = vec![
        1.1959364270e-1,    3.5878092050e-1,    1.1959364270e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.8445838930e1,     2.6629886630e0,    -1.8445838930e1,
    ];
    let expected_extremal_frequencies = vec![
        8.8709652420e-2,    1.5000000600e-1,    2.0000000300e-1,    4.0000000600e-1,
        4.5241931080e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
8,1,3,32
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n8_grid32() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(8, FilterType::MultipleBand, &bands, 32);

    let expected_impulse_response = vec![
        1.6935145850e-1,   -1.0425144430e-1,   -2.8712207080e-1,    2.2249950470e-1,
    ];
    let expected_deviations = vec![
        1.1956638100e-1,    3.5869914290e-1,    1.1956638100e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.8447818760e1,     2.6624655720e0,    -1.8447818760e1,
    ];
    let expected_extremal_frequencies = vec![
        8.9843750000e-2,    1.5000000600e-1,    2.0000000300e-1,    4.0000000600e-1,
        4.5468750600e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,3,16
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n16_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -4.0730793030e-2,   -1.2830013410e-2,   -1.2071261180e-2,   -1.6163289550e-3,
        1.7818668480e-1,   -1.0700613260e-1,   -3.0037075280e-1,    2.5910782810e-1,
    ];
    let expected_deviations = vec![
        7.4661538000e-2,    2.2398461400e-1,    7.4661538000e-2,
    ];
    let expected_deviation_dbs = vec![
        -2.2538061140e1,     1.7555190320e0,    -2.2538061140e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.2500000000e-2,    1.2109375000e-1,    1.5000000600e-1,
        2.0000000300e-1,    2.8593748810e-1,    3.4999999400e-1,    4.0000000600e-1,
        4.3515625600e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,3,17
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n16_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        -4.0728058670e-2,   -1.2864863500e-2,   -1.2010248380e-2,   -1.6363300380e-3,
        1.7819005250e-1,   -1.0701917110e-1,   -3.0033701660e-1,    2.5908222790e-1,
    ];
    let expected_deviations = vec![
        7.4646823110e-2,    2.2394046190e-1,    7.4646823110e-2,
    ];
    let expected_deviation_dbs = vec![
        -2.2539772030e1,     1.7552061080e0,    -2.2539772030e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.2500000000e-2,    1.2132358550e-1,    1.5000000600e-1,
        2.0000000300e-1,    2.8455889230e-1,    3.4999999400e-1,    4.0000000600e-1,
        4.3308827280e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,3,31
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n16_grid31() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 31);

    let expected_impulse_response = vec![
        -4.0729813280e-2,   -1.2838121500e-2,   -1.2059140950e-2,   -1.6224272550e-3,
        1.7816682160e-1,   -1.0699984430e-1,   -3.0034434800e-1,    2.5909090040e-1,
    ];
    let expected_deviations = vec![
        7.4671939020e-2,    2.2401580210e-1,    7.4671939020e-2,
    ];
    let expected_deviation_dbs = vec![
        -2.2536851880e1,     1.7557406430e0,    -2.2536851880e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.2499970200e-2,    1.2096765640e-1,    1.5000000600e-1,
        2.0000000300e-1,    2.8467735650e-1,    3.4999999400e-1,    4.0000000600e-1,
        4.3427416680e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,3,32
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n16_grid32() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 32);

    let expected_impulse_response = vec![
        -4.0737364440e-2,   -1.2833038340e-2,   -1.2054493650e-2,   -1.6187988220e-3,
        1.7816303670e-1,   -1.0699553790e-1,   -3.0034494400e-1,    2.5908750300e-1,
    ];
    let expected_deviations = vec![
        7.4667297300e-2,    2.2400189940e-1,    7.4667297300e-2,
    ];
    let expected_deviation_dbs = vec![
        -2.2537391660e1,     1.7556418180e0,    -2.2537391660e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.2500000000e-2,    1.2109375000e-1,    1.5000000600e-1,
        2.0000000300e-1,    2.8398436310e-1,    3.4999999400e-1,    4.0000000600e-1,
        4.3515625600e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
37,1,3,16
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n37_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(37, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -7.6172472910e-3,    3.7963695360e-3,    2.6189079510e-3,    6.8600727250e-3,
        1.2505672870e-2,   -2.6415402070e-2,   -8.9944992210e-3,    1.9582133740e-2,
        -5.0852075220e-5,    3.5998895760e-2,   -1.8624186520e-2,   -6.3116729260e-2,
        2.8667811300e-2,   -5.4394663310e-3,    7.7281787990e-2,    9.0501762930e-2,
        -2.8178721670e-1,   -5.7383574550e-2,    3.9199963210e-1,
    ];
    let expected_deviations = vec![
        8.7680816650e-3,    2.6304245000e-2,    8.7680816650e-3,
    ];
    let expected_deviation_dbs = vec![
        -4.1141910550e1,    2.2552251820e-1,    -4.1141910550e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    4.2763158680e-2,    7.2368443010e-2,    9.8684251310e-2,
        1.2335532160e-1,    1.4144736530e-1,    1.5000000600e-1,    2.0000000300e-1,
        2.1151311700e-1,    2.3618407550e-1,    2.7236816290e-1,    3.1184169650e-1,
        3.3815738560e-1,    3.4999999400e-1,    4.0000000600e-1,    4.0822365880e-1,
        4.2631569500e-1,    4.4934192300e-1,    4.7401288150e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
37,1,3,17
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n37_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(37, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        -7.6094171960e-3,    3.7917175800e-3,    2.6155025700e-3,    6.8725580350e-3,
        1.2501673770e-2,   -2.6417639110e-2,   -8.9805098250e-3,    1.9570332020e-2,
        -5.5847922340e-5,    3.6010079090e-2,   -1.8633922560e-2,   -6.3116565350e-2,
        2.8668329120e-2,   -5.4463143460e-3,    7.7288173140e-2,    9.0503647920e-2,
        -2.8179275990e-1,   -5.7383097710e-2,    3.9199757580e-1,
    ];
    let expected_deviations = vec![
        8.7693938990e-3,    2.6308180760e-2,    8.7693938990e-3,
    ];
    let expected_deviation_dbs = vec![
        -4.1140609740e1,    2.2555580740e-1,    -4.1140609740e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    4.1795678440e-2,    7.1207441390e-2,    9.9071167410e-2,
        1.2229093910e-1,    1.4241482320e-1,    1.5000000600e-1,    2.0000000300e-1,
        2.1083594860e-1,    2.3715181650e-1,    2.7120763060e-1,    3.1145542860e-1,
        3.3777129650e-1,    3.4999999400e-1,    4.0000000600e-1,    4.0773996710e-1,
        4.2631587390e-1,    4.4953575730e-1,    4.7430363300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
37,1,3,31
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n37_grid31() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(37, FilterType::MultipleBand, &bands, 31);

    let expected_impulse_response = vec![
        -7.6259085910e-3,    3.7832176310e-3,    2.6211775840e-3,    6.8616382780e-3,
        1.2505715710e-2,   -2.6411056520e-2,   -8.9943744240e-3,    1.9577264790e-2,
        -5.3560052040e-5,    3.6004431550e-2,   -1.8613712860e-2,   -6.3115820290e-2,
        2.8661366550e-2,   -5.4390868170e-3,    7.7285215260e-2,    9.0512409810e-2,
        -2.8178489210e-1,   -5.7386219500e-2,    3.9199790360e-1,
    ];
    let expected_deviations = vec![
        8.7735177950e-3,    2.6320552450e-2,    8.7735177950e-3,
    ];
    let expected_deviation_dbs = vec![
        -4.1136528020e1,    2.2566072640e-1,    -4.1136528020e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    4.2444836350e-2,    7.2156220670e-2,    9.8472021520e-2,
        1.2309002880e-1,    1.4176560940e-1,    1.5000000600e-1,    2.0000000300e-1,
        2.1103556450e-1,    2.3650224510e-1,    2.7215561270e-1,    3.1120452280e-1,
        3.3836898210e-1,    3.4999999400e-1,    4.0000000600e-1,    4.0764001010e-1,
        4.2631557580e-1,    4.5008447770e-1,    4.7470226880e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
37,1,3,32
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n37_grid32() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(37, FilterType::MultipleBand, &bands, 32);

    let expected_impulse_response = vec![
        -7.6231346470e-3,    3.7855841220e-3,    2.6194953830e-3,    6.8615111520e-3,
        1.2505315240e-2,   -2.6411429050e-2,   -8.9920135210e-3,    1.9577020780e-2,
        -5.2807481550e-5,    3.6005057390e-2,   -1.8617995080e-2,   -6.3117623330e-2,
        2.8663374480e-2,   -5.4392977620e-3,    7.7285490930e-2,    9.0510442850e-2,
        -2.8178697820e-1,   -5.7385645810e-2,    3.9199846980e-1,
    ];
    let expected_deviations = vec![
        8.7711857630e-3,    2.6313558220e-2,    8.7711857630e-3,
    ];
    let expected_deviation_dbs = vec![
        -4.1138835910e1,    2.2560119630e-1,    -4.1138835910e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    4.2763173580e-2,    7.1546047930e-2,    9.8684102300e-2,
        1.2253269550e-1,    1.4226946230e-1,    1.5000000600e-1,    2.0000000300e-1,
        2.1151311700e-1,    2.3700644080e-1,    2.7154579760e-1,    3.1101933120e-1,
        3.3815738560e-1,    3.4999999400e-1,    4.0000000600e-1,    4.0822365880e-1,
        4.2631569500e-1,    4.5016428830e-1,    4.7483524680e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
38,1,3,16
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n38_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(38, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -6.8362783640e-3,   -3.1098246110e-3,    6.6664060580e-3,   -1.9826511850e-3,
        1.4424593190e-2,   -5.5908225480e-3,   -2.9471138490e-2,    1.3240320610e-2,
        6.9278236480e-3,    1.6568847000e-2,    2.5513689970e-2,   -6.5501570700e-2,
        -1.3812240210e-2,    2.4961978200e-2,    3.6959573630e-3,    1.4203712340e-1,
        -9.6060618760e-2,   -2.8610336780e-1,    2.5024831300e-1,
    ];
    let expected_deviations = vec![
        8.3669349550e-3,    2.5100806730e-2,    8.3669349550e-3,
    ];
    let expected_deviation_dbs = vec![
        -4.1548671720e1,    2.1533168850e-1,    -4.1548671720e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.6315785940e-2,    5.0986848770e-2,    7.7302657070e-2,
        1.0197372730e-1,    1.2335532160e-1,    1.4144736530e-1,    1.5000000600e-1,
        2.0000000300e-1,    2.1151311700e-1,    2.4111826720e-1,    2.8388127680e-1,
        3.1513115760e-1,    3.3980211620e-1,    3.4999999400e-1,    4.0000000600e-1,
        4.0822365880e-1,    4.2796042560e-1,    4.5263138410e-1,    4.8388126490e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
38,1,3,17
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n38_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(38, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        -6.8096276370e-3,   -3.1368408820e-3,    6.6216918640e-3,   -1.9723719450e-3,
        1.4458631170e-2,   -5.5700596420e-3,   -2.9486458750e-2,    1.3217032890e-2,
        6.9154920060e-3,    1.6566269100e-2,    2.5549381970e-2,   -6.5472230320e-2,
        -1.3849459590e-2,    2.4931468070e-2,    3.6959405990e-3,    1.4205859600e-1,
        -9.6033245330e-2,   -2.8610584140e-1,    2.5022789840e-1,
    ];
    let expected_deviations = vec![
        8.3874911070e-3,    2.5162473320e-2,    8.3874911070e-3,
    ];
    let expected_deviation_dbs = vec![
        -4.1527359010e1,    2.1585388480e-1,    -4.1527359010e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.6315793400e-2,    5.1083609460e-2,    7.7399380510e-2,
        1.0061915220e-1,    1.2383892390e-1,    1.4241482320e-1,    1.5000000600e-1,
        2.0000000300e-1,    2.1083594860e-1,    2.4024780090e-1,    2.8513956070e-1,
        3.1609940530e-1,    3.3931928870e-1,    3.4999999400e-1,    4.0000000600e-1,
        4.0773996710e-1,    4.2786386610e-1,    4.5263174180e-1,    4.8359158640e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
38,1,3,31
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n38_grid31() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(38, FilterType::MultipleBand, &bands, 31);

    let expected_impulse_response = vec![
        -6.8081570790e-3,   -3.1372685920e-3,    6.6267056390e-3,   -1.9656831860e-3,
        1.4446566810e-2,   -5.5819302800e-3,   -2.9482107610e-2,    1.3228946370e-2,
        6.9191576910e-3,    1.6558811070e-2,    2.5543684140e-2,   -6.5471574660e-2,
        -1.3841167090e-2,    2.4932477620e-2,    3.6943852900e-3,    1.4205491540e-1,
        -9.6035934980e-2,   -2.8611171250e-1,    2.5023475290e-1,
    ];
    let expected_deviations = vec![
        8.3902021870e-3,    2.5170607490e-2,    8.3902021870e-3,
    ];
    let expected_deviation_dbs = vec![
        -4.1524551390e1,    2.1592257920e-1,    -4.1524551390e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.5466900320e-2,    5.1782701160e-2,    7.6400704680e-2,
        1.0101871190e-1,    1.2393892560e-1,    1.4261449870e-1,    1.5000000600e-1,
        2.0000000300e-1,    2.1188445390e-1,    2.4074669180e-1,    2.8488895300e-1,
        3.1544896960e-1,    3.3921787140e-1,    3.4999999400e-1,    4.0000000600e-1,
        4.0848889950e-1,    4.2801335450e-1,    4.5348003510e-1,    4.8319116230e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
38,1,3,32
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n38_grid32() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(38, FilterType::MultipleBand, &bands, 32);

    let expected_impulse_response = vec![
        -6.8058418110e-3,   -3.1393023670e-3,    6.6152093930e-3,   -1.9587883730e-3,
        1.4459590430e-2,   -5.5966842920e-3,   -2.9481895270e-2,    1.3237884270e-2,
        6.9118775430e-3,    1.6556246210e-2,    2.5546018030e-2,   -6.5464228390e-2,
        -1.3838995250e-2,    2.4919902910e-2,    3.6960132420e-3,    1.4206035440e-1,
        -9.6038036050e-2,   -2.8610739110e-1,    2.5023055080e-1,
    ];
    let expected_deviations = vec![
        8.3950795230e-3,    2.5185238570e-2,    8.3950795230e-3,
    ];
    let expected_deviation_dbs = vec![
        -4.1519504550e1,    2.1604679520e-1,    -4.1519504550e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.5493424390e-2,    5.1809232680e-2,    7.6480239630e-2,
        1.0115119810e-1,    1.2417742610e-1,    1.4226946230e-1,    1.5000000600e-1,
        2.0000000300e-1,    2.1151311700e-1,    2.4111826720e-1,    2.8470364210e-1,
        3.1513115760e-1,    3.3897975090e-1,    3.4999999400e-1,    4.0000000600e-1,
        4.0822365880e-1,    4.2796042560e-1,    4.5345374940e-1,    4.8388126490e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
95,1,3,16
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n95_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(95, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -2.9520822860e-5,    2.3667080030e-7,   -7.5317635490e-5,    2.1481959270e-4,
        1.8214067680e-4,   -3.4650936140e-4,   -2.5621537130e-5,   -2.8325253520e-4,
        1.2241059450e-4,    1.1638670690e-3,   -5.3156953070e-4,   -4.7182169510e-4,
        -4.2616098650e-4,   -1.1816986370e-3,    2.6595303790e-3,    9.4058021200e-4,
        -1.5394178920e-3,    8.9186374910e-7,   -4.0531028060e-3,    1.9056017040e-3,
        5.6179882960e-3,   -2.0424798130e-3,    1.1766594830e-3,   -6.1203781520e-3,
        -5.4070539770e-3,    1.1434163900e-2,    1.1308406250e-3,    2.3385980170e-3,
        -2.0226598720e-3,   -2.0006397740e-2,    9.8452921960e-3,    1.0823925030e-2,
        2.1308881700e-3,    1.3120249850e-2,   -3.2824974510e-2,   -1.2230303140e-2,
        2.5098077950e-2,   -4.5509460730e-7,    3.9229977880e-2,   -2.1495651450e-2,
        -7.0486180480e-2,    3.2794147730e-2,   -2.8186314740e-3,    7.6474621890e-2,
        9.1023191810e-2,   -2.8493866320e-1,   -5.7976778600e-2,    3.9587539430e-1,
    ];
    let expected_deviations = vec![
        6.3613268140e-5,    1.9083981170e-4,    6.3613268140e-5,
    ];
    let expected_deviation_dbs = vec![
        -8.3929046630e1,    1.6575793270e-3,    -8.3929046630e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.2369792910e-2,    2.4088535460e-2,    3.5807285460e-2,
        4.6875014900e-2,    5.7942744340e-2,    6.9010436530e-2,    7.9427063470e-2,
        8.9843690400e-2,    1.0026031730e-1,    1.1002590510e-1,    1.1914045360e-1,
        1.2825503950e-1,    1.3671864570e-1,    1.4322911200e-1,    1.4843748510e-1,
        1.5000000600e-1,    2.0000000300e-1,    2.0195314290e-1,    2.0716151600e-1,
        2.1497407560e-1,    2.2408872840e-1,    2.3320338130e-1,    2.4362012740e-1,
        2.5403678420e-1,    2.6445329190e-1,    2.7486979960e-1,    2.8528630730e-1,
        2.9570281510e-1,    3.0611932280e-1,    3.1653583050e-1,    3.2630130650e-1,
        3.3476471900e-1,    3.4257709980e-1,    3.4778535370e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.0195310120e-1,    4.0651032330e-1,    4.1302064060e-1,
        4.2148405310e-1,    4.3059849740e-1,    4.3971294160e-1,    4.4947841760e-1,
        4.5924389360e-1,    4.6966040130e-1,    4.7942587730e-1,    4.8984238510e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
95,1,3,17
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n95_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(95, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        -2.9664293830e-5,   -3.8374687960e-7,   -7.5307456430e-5,    2.1541309250e-4,
        1.8207459650e-4,   -3.4576066540e-4,   -2.6088089730e-5,   -2.8497583120e-4,
        1.2347727900e-4,    1.1640433220e-3,   -5.3049443520e-4,   -4.6980188930e-4,
        -4.3020580780e-4,   -1.1824137760e-3,    2.6607133910e-3,    9.3988084700e-4,
        -1.5331108590e-3,   -1.2152954700e-6,   -4.0589608250e-3,    1.9076979950e-3,
        5.6139370430e-3,   -2.0356148020e-3,    1.1836092450e-3,   -6.1302492400e-3,
        -5.4062693380e-3,    1.1428554540e-2,    1.1295349100e-3,    2.3563425060e-3,
        -2.0284319760e-3,   -2.0009884610e-2,    9.8433177920e-3,    1.0808963330e-2,
        2.1490459330e-3,    1.3128103690e-2,   -3.2832641150e-2,   -1.2224876320e-2,
        2.5075970220e-2,    1.5765909890e-6,    3.9251621810e-2,   -2.1502446380e-2,
        -7.0475645360e-2,    3.2780669630e-2,   -2.8396565470e-3,    7.6496854420e-2,
        9.1023363170e-2,   -2.8493040800e-1,   -5.7970177380e-2,    3.9584341650e-1,
    ];
    let expected_deviations = vec![
        6.3558072720e-5,    1.9067421090e-4,    6.3558072720e-5,
    ];
    let expected_deviation_dbs = vec![
        -8.3936584470e1,    1.6555087640e-3,    -8.3936584470e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.2254902160e-2,    2.3897059260e-2,    3.5539228470e-2,
        4.7181420030e-2,    5.8210864660e-2,    6.8627521400e-2,    7.9656898980e-2,
        9.0073533360e-2,    9.9877424540e-2,    1.0968131570e-1,    1.1948520690e-1,
        1.2867639960e-1,    1.3664215800e-1,    1.4338241520e-1,    1.4828442040e-1,
        1.5000000600e-1,    2.0000000300e-1,    2.0183825490e-1,    2.0735301080e-1,
        2.1470601860e-1,    2.2389727830e-1,    2.3370128870e-1,    2.4350529910e-1,
        2.5392195580e-1,    2.6433846350e-1,    2.7475497130e-1,    2.8517147900e-1,
        2.9558798670e-1,    3.0600449440e-1,    3.1642100210e-1,    3.2622477410e-1,
        3.3480307460e-1,    3.4276863930e-1,    3.4767052530e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.0183821320e-1,    4.0612736340e-1,    4.1348019240e-1,
        4.2144575720e-1,    4.3002405760e-1,    4.3982782960e-1,    4.4963160160e-1,
        4.5943537350e-1,    4.6923914550e-1,    4.7965565320e-1,    4.9007216100e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
98,1,3,16
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n98_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(98, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        3.7872963730e-5,   -7.9897159590e-5,   -2.8738460970e-5,   -3.6309800630e-5,
        7.6334545160e-6,    3.9500708230e-4,   -1.7148815100e-4,   -3.3938128040e-4,
        4.3183463280e-6,   -3.9681443010e-4,    9.7908789760e-4,    5.9965142280e-4,
        -1.1005694980e-3,   -1.6796693670e-5,   -1.3342448510e-3,    5.7931907940e-4,
        3.1737214890e-3,   -1.3653531200e-3,   -4.0886178610e-4,   -1.9577792850e-3,
        -2.9951045290e-3,    6.3653560360e-3,    1.3942229560e-3,   -1.3816426510e-3,
        -3.0199973840e-4,   -1.0021563620e-2,    4.7997208310e-3,    9.1188140210e-3,
        -1.9774925900e-3,    5.0668409090e-3,   -1.5168691050e-2,   -8.9402934540e-3,
        1.8635857850e-2,    5.4600415750e-4,    1.2495992710e-2,   -6.9810301070e-3,
        -3.7359789010e-2,    1.8107928340e-2,    9.8406467590e-3,    1.6073992480e-2,
        2.9666289690e-2,   -7.2061873970e-2,   -1.7113812270e-2,    3.0740382150e-2,
        3.5031717270e-3,    1.4265882970e-1,   -9.6291027960e-2,   -2.9007321600e-1,
        2.5308877230e-1,
    ];
    let expected_deviations = vec![
        4.8679103200e-5,    1.4603731690e-4,    4.8679103200e-5,
    ];
    let expected_deviation_dbs = vec![
        -8.6253143310e1,    1.2683197860e-3,    -8.6253143310e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0204081420e-2,    2.1045913920e-2,    3.1249990690e-2,
        4.1454069320e-2,    5.2295900880e-2,    6.2499977650e-2,    7.2704054420e-2,
        8.2270376380e-2,    9.2474453150e-2,    1.0204077510e-1,    1.1160709710e-1,
        1.2117341910e-1,    1.2946423890e-1,    1.3711729650e-1,    1.4413259920e-1,
        1.4859688280e-1,    1.5000000600e-1,    2.0000000300e-1,    2.0191326740e-1,
        2.0701530580e-1,    2.1403060850e-1,    2.2295917570e-1,    2.3188774290e-1,
        2.4209181960e-1,    2.5229594110e-1,    2.6313802600e-1,    2.7398011090e-1,
        2.8482219580e-1,    2.9566428070e-1,    3.0650636550e-1,    3.1671068070e-1,
        3.2627722620e-1,    3.3520600200e-1,    3.4285923840e-1,    3.4796139600e-1,
        3.4999999400e-1,    4.0000000600e-1,    4.0191331510e-1,    4.0573993330e-1,
        4.1275539990e-1,    4.2040863630e-1,    4.2869964240e-1,    4.3762841820e-1,
        4.4655719400e-1,    4.5612373950e-1,    4.6569028500e-1,    4.7589460020e-1,
        4.8546114560e-1,    4.9502769110e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
98,1,3,17
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n98_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(98, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        3.7625188270e-5,   -7.9702556830e-5,   -2.8217567890e-5,   -3.7250603780e-5,
        7.6863216240e-6,    3.9527600170e-4,   -1.7182310690e-4,   -3.3740358780e-4,
        2.8910362740e-6,   -3.9841630500e-4,    9.8060746680e-4,    5.9783027970e-4,
        -1.0969731960e-3,   -1.5501776940e-5,   -1.3398814480e-3,    5.8161257770e-4,
        3.1709910840e-3,   -1.3637496160e-3,   -3.9999524600e-4,   -1.9660070540e-3,
        -2.9956551730e-3,    6.3646789640e-3,    1.3880212790e-3,   -1.3655347280e-3,
        -3.0435528610e-4,   -1.0029329920e-2,    4.8039485700e-3,    9.1032078490e-3,
        -1.9644265990e-3,    5.0802165640e-3,   -1.5182683240e-2,   -8.9322505520e-3,
        1.8618740140e-2,    5.4074358200e-4,    1.2524879540e-2,   -6.9918222730e-3,
        -3.7354007360e-2,    1.8102917820e-2,    9.8116956650e-3,    1.6102306540e-2,
        2.9670331630e-2,   -7.2064355020e-2,   -1.7099507150e-2,    3.0701844020e-2,
        3.5093966870e-3,    1.4267985520e-1,   -9.6301086250e-2,   -2.9004710910e-1,
        2.5306531790e-1,
    ];
    let expected_deviations = vec![
        4.8876485380e-5,    1.4662944890e-4,    4.8876485380e-5,
    ];
    let expected_deviation_dbs = vec![
        -8.6218002320e1,    1.2734963090e-3,    -8.6218002320e1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0204082350e-2,    2.1008398380e-2,    3.1212465840e-2,
        4.1416563090e-2,    5.2220903340e-2,    6.2425002460e-2,    7.2629101570e-2,
        8.2833200690e-2,    9.2437058690e-2,    1.0204091670e-1,    1.1164477470e-1,
        1.2064839150e-1,    1.2965194880e-1,    1.3745498660e-1,    1.4405755700e-1,
        1.4825919270e-1,    1.5000000600e-1,    2.0000000300e-1,    2.0180070400e-1,
        2.0720280710e-1,    2.1440561120e-1,    2.2280888260e-1,    2.3241262140e-1,
        2.4201636020e-1,    2.5222039220e-1,    2.6302486660e-1,    2.7382934090e-1,
        2.8463381530e-1,    2.9543828960e-1,    3.0624276400e-1,    3.1704723830e-1,
        3.2665121560e-1,    3.3565494420e-1,    3.4285792710e-1,    3.4826016430e-1,
        3.4999999400e-1,    4.0000000600e-1,    4.0180075170e-1,    4.0600249170e-1,
        4.1260522600e-1,    4.2040845750e-1,    4.2881193760e-1,    4.3781566620e-1,
        4.4681939480e-1,    4.5642337200e-1,    4.6602734920e-1,    4.7563132640e-1,
        4.8523530360e-1,    4.9483928080e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
128,1,3,16
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n128_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -1.1280264520e-5,    5.3813237170e-6,    1.0350037880e-5,    3.3040532800e-6,
        2.3351432900e-5,   -5.7080935220e-5,   -3.3578817240e-5,    6.2335369880e-5,
        -8.7662738220e-7,    1.1213921970e-4,   -4.8340938520e-5,   -2.5134225140e-4,
        1.0659681720e-4,    6.5233907660e-6,    2.1438227850e-4,    3.0148442600e-4,
        -6.3823367240e-4,   -1.2930514640e-4,    6.5141241070e-5,    4.2339379430e-5,
        1.2176329040e-3,   -5.7497870880e-4,   -1.0252343490e-3,    1.3824434430e-4,
        -7.7252014310e-4,    2.0940138490e-3,    1.1741175550e-3,   -2.3315637370e-3,
        -4.2663887140e-5,   -2.0767236130e-3,    1.0236850940e-3,    5.1212431860e-3,
        -2.3055702910e-3,   -8.7373727000e-4,   -2.6098145170e-3,   -4.0185875260e-3,
        8.7886350230e-3,    1.8330109310e-3,   -2.2400161250e-3,   -4.0783686560e-4,
        -1.2158143330e-2,    5.9895291920e-3,    1.1113859710e-2,   -2.6911429600e-3,
        5.4483828130e-3,   -1.7028149220e-2,   -9.8587907850e-3,    2.0993664860e-2,
        6.0055404900e-4,    1.2860749850e-2,   -7.5209178030e-3,   -3.9439935240e-2,
        1.9329000260e-2,    1.0712904860e-2,    1.6119496900e-2,    3.0058022590e-2,
        -7.3577865960e-2,   -1.7307654020e-2,    3.1788956370e-2,    3.5417750480e-3,
        1.4278630910e-1,   -9.6701867880e-2,   -2.9059514400e-1,    2.5364086030e-1,
    ];
    let expected_deviations = vec![
        3.1094371020e-6,    9.3283115350e-6,    3.1094371020e-6,
    ];
    let expected_deviation_dbs = vec![
        -1.1014636230e2,    8.0763849840e-5,    -1.1014636230e2,
    ];
    let expected_extremal_frequencies = vec![
        7.3242187500e-3,    1.5136718750e-2,    2.2460937500e-2,    3.0273437500e-2,
        3.7597656250e-2,    4.5410156250e-2,    5.2734375000e-2,    6.0546875000e-2,
        6.7871093750e-2,    7.5195312500e-2,    8.2519531250e-2,    8.9843750000e-2,
        9.7167968750e-2,    1.0449218750e-1,    1.1132812500e-1,    1.1816406250e-1,
        1.2500000000e-1,    1.3134765620e-1,    1.3720703120e-1,    1.4208984380e-1,
        1.4648437500e-1,    1.4892578120e-1,    1.5000000600e-1,    2.0000000300e-1,
        2.0097656550e-1,    2.0390625300e-1,    2.0878906550e-1,    2.1416015920e-1,
        2.2050781550e-1,    2.2734375300e-1,    2.3466797170e-1,    2.4199219050e-1,
        2.4931640920e-1,    2.5761717560e-1,    2.6542967560e-1,    2.7421873810e-1,
        2.8300780060e-1,    2.9130858180e-1,    2.9960936310e-1,    3.0791014430e-1,
        3.1523436310e-1,    3.2255858180e-1,    3.2939451930e-1,    3.3574217560e-1,
        3.4111326930e-1,    3.4599608180e-1,    3.4892576930e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.0097656850e-1,    4.0390625600e-1,    4.0781250600e-1,
        4.1269531850e-1,    4.1855469350e-1,    4.2490234970e-1,    4.3173828720e-1,
        4.3857422470e-1,    4.4541016220e-1,    4.5224609970e-1,    4.5957031850e-1,
        4.6689453720e-1,    4.7421875600e-1,    4.8154297470e-1,    4.8886719350e-1,
        4.9619141220e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,1,3,17
0.000000000000,0.150000005960,0.200000002980,0.349999994040
0.400000005960,0.500000000000
0.000000000000,1.000000000000,0.000000000000
3.000000000000,1.000000000000,3.000000000000
*/
#[test]
fn multipleband_bandpass_n128_grid17() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.35f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 3f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 17);

    let expected_impulse_response = vec![
        -1.1234141310e-5,    5.4508009270e-6,    1.0288878910e-5,    3.2620250750e-6,
        2.3342918210e-5,   -5.7072120400e-5,   -3.3329761210e-5,    6.2336563130e-5,
        -1.0683161240e-6,    1.1209334480e-4,   -4.8575588150e-5,   -2.5094015290e-4,
        1.0705220480e-4,    6.1607133830e-6,    2.1432938230e-4,    3.0087959020e-4,
        -6.3829350980e-4,   -1.2802076530e-4,    6.5001309850e-5,    4.2310653950e-5,
        1.2170417470e-3,   -5.7645112970e-4,   -1.0235883530e-3,    1.3904941440e-4,
        -7.7255361250e-4,    2.0944685680e-3,    1.1711033290e-3,   -2.3313201960e-3,
        -4.0507758970e-5,   -2.0768912510e-3,    1.0259652040e-3,    5.1185726190e-3,
        -2.3086131550e-3,   -8.7102199900e-4,   -2.6102664410e-3,   -4.0153302250e-3,
        8.7895635520e-3,    1.8269368450e-3,   -2.2386692000e-3,   -4.0839123540e-4,
        -1.2156703510e-2,    5.9957355260e-3,    1.1108291340e-2,   -2.6928684670e-3,
        5.4484922440e-3,   -1.7031211410e-2,   -9.8500624300e-3,    2.0993333310e-2,
        5.9599988160e-4,    1.2862394560e-2,   -7.5280088930e-3,   -3.9434932170e-2,
        1.9335599620e-2,    1.0707909240e-2,    1.6122702510e-2,    3.0051408340e-2,
        -7.3581427340e-2,   -1.7297763380e-2,    3.1786378470e-2,    3.5450067370e-3,
        1.4278556410e-1,   -9.6712566910e-2,   -2.9058849810e-1,    2.5364184380e-1,
    ];
    let expected_deviations = vec![
        3.1133174620e-6,    9.3399521570e-6,    3.1133174620e-6,
    ];
    let expected_deviation_dbs = vec![
        -1.1013553620e2,    8.0763849840e-5,    -1.1013553620e2,
    ];
    let expected_extremal_frequencies = vec![
        6.8933824080e-3,    1.4705888930e-2,    2.2518396380e-2,    3.0330903830e-2,
        3.7683852020e-2,    4.5496359470e-2,    5.2849307660e-2,    6.0202255850e-2,
        6.8014763300e-2,    7.5367711480e-2,    8.2720659670e-2,    9.0073607860e-2,
        9.7426556050e-2,    1.0431994500e-1,    1.1167289320e-1,    1.1810672280e-1,
        1.2500010430e-1,    1.3143382970e-1,    1.3694845140e-1,    1.4200352130e-1,
        1.4613948760e-1,    1.4889679850e-1,    1.5000000600e-1,    2.0000000300e-1,
        2.0091910660e-1,    2.0413596930e-1,    2.0873148740e-1,    2.1424610910e-1,
        2.2067983450e-1,    2.2757311170e-1,    2.3446638880e-1,    2.4181921780e-1,
        2.4963159860e-1,    2.5744399430e-1,    2.6571592690e-1,    2.7398785950e-1,
        2.8271934390e-1,    2.9145082830e-1,    2.9972276090e-1,    3.0753514170e-1,
        3.1534752250e-1,    3.2270035150e-1,    3.2913407680e-1,    3.3556780220e-1,
        3.4154197570e-1,    3.4567794200e-1,    3.4889480470e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.0091910960e-1,    4.0367642050e-1,    4.0781238680e-1,
        4.1286745670e-1,    4.1884163020e-1,    4.2481580380e-1,    4.3170908090e-1,
        4.3814280630e-1,    4.4549563530e-1,    4.5238891240e-1,    4.5974174140e-1,
        4.6709457040e-1,    4.7398784760e-1,    4.8134067650e-1,    4.8915305730e-1,
        4.9650588630e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
9,1,2,15
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n9_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(9, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        1.1682424690e-2,   -1.4297962190e-1,    1.7489455640e-1,   -2.3147127030e-1,
        2.5129693750e-1,
    ];
    let expected_deviations = vec![
        1.2445089970e-1,    3.7335270640e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.8100038530e1,     2.7556421760e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5333332120e-1,    2.8666654230e-1,    3.4999999400e-1,
        4.0000000600e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
9,1,2,16
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n9_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(9, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        1.1621468700e-2,   -1.4301973580e-1,    1.7492446300e-1,   -2.3141521220e-1,
        2.5134310130e-1,
    ];
    let expected_deviations = vec![
        1.2443497030e-1,    3.7330490350e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.8101150510e1,     2.7553391460e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5624998510e-1,    2.8749987480e-1,    3.4999999400e-1,
        4.0000000600e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
9,1,2,34
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n9_grid34() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(9, FilterType::MultipleBand, &bands, 34);

    let expected_impulse_response = vec![
        1.1708342470e-2,   -1.4298257230e-1,    1.7485940460e-1,   -2.3147755860e-1,
        2.5132462380e-1,
    ];
    let expected_deviations = vec![
        1.2446011600e-1,    3.7338036300e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.8099395750e1,     2.7558171750e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5294118230e-1,    2.8529429440e-1,    3.4999999400e-1,
        4.0000000600e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
9,1,2,35
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n9_grid35() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(9, FilterType::MultipleBand, &bands, 35);

    let expected_impulse_response = vec![
        1.1701596900e-2,   -1.4298683400e-1,    1.7486292120e-1,   -2.3147150870e-1,
        2.5132930280e-1,
    ];
    let expected_deviations = vec![
        1.2445833530e-1,    3.7337499860e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.8099519730e1,     2.7557826040e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5428572890e-1,    2.8571456670e-1,    3.4999999400e-1,
        4.0000000600e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
10,1,2,15
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n10_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(10, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        2.2558239100e-1,   -2.1365028620e-1,    2.2796037790e-1,   -1.7560601230e-1,
        6.5808534620e-2,
    ];
    let expected_deviations = vec![
        2.6019006970e-1,    7.8057020900e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.1694185260e1,     5.0111823080e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0666667670e-1,    2.0666660370e-1,    2.9999986290e-1,
        3.4999999400e-1,    4.9333325030e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
10,1,2,16
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n10_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(10, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        2.2863647340e-1,   -2.1660661700e-1,    2.3109629750e-1,   -1.7797660830e-1,
        6.6746354100e-2,
    ];
    let expected_deviations = vec![
        2.6379191880e-1,    7.9137581590e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.1574771880e1,     5.0637345310e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0625001790e-1,    2.0624993740e-1,    2.9999986290e-1,
        3.4999999400e-1,    4.9374991660e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
10,1,2,34
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n10_grid34() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(10, FilterType::MultipleBand, &bands, 34);

    let expected_impulse_response = vec![
        2.3136226830e-1,   -2.1931596100e-1,    2.3391616340e-1,   -1.8008559940e-1,
        6.7657053470e-2,
    ];
    let expected_deviations = vec![
        2.6706761120e-1,    8.0120283370e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.1467576030e1,     5.1112527850e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0294117780e-1,    2.0588235560e-1,    3.0000025030e-1,
        3.4999999400e-1,    4.9411812420e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
10,1,2,35
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n10_grid35() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(10, FilterType::MultipleBand, &bands, 35);

    let expected_impulse_response = vec![
        2.3267714680e-1,   -2.2055725750e-1,    2.3524519800e-1,   -1.8110948800e-1,
        6.8035542960e-2,
    ];
    let expected_deviations = vec![
        2.6858186720e-1,    8.0574566130e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.1418466570e1,     5.1331310270e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.0285710540e-1,    2.0571440460e-1,    3.0000030990e-1,
        3.4999999400e-1,    4.9428591130e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,15
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n16_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        -1.7613439260e-1,    1.1915747820e-1,   -1.4936113360e-1,    1.7100140450e-1,
        -1.7483925820e-1,    1.5374338630e-1,   -1.0594111680e-1,    3.7968516350e-2,
    ];
    let expected_deviations = vec![
        2.4881014230e-1,    7.4643045660e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.2082638740e1,     4.8430256840e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.6666670140e-2,    1.2916670740e-1,    1.9583331050e-1,
        2.5833326580e-1,    3.2083320620e-1,    4.0000000600e-1,    4.5416662100e-1,
        4.9583324790e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,16
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n16_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -1.5233190360e-1,    1.0389186440e-1,   -1.3624292610e-1,    1.6400679950e-1,
        -1.7501091960e-1,    1.5884864330e-1,   -1.1163556580e-1,    4.0271520610e-2,
    ];
    let expected_deviations = vec![
        2.1640506390e-1,    6.4921522140e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.3294651030e1,     4.3455467220e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.6406250000e-2,    1.2890625000e-1,    1.9531250000e-1,
        2.5781250000e-1,    3.2421875000e-1,    4.0000000600e-1,    4.5468750600e-1,
        4.9375000600e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,34
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n16_grid34() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 34);

    let expected_impulse_response = vec![
        -1.9716601070e-1,    1.3285900650e-1,   -1.6108924150e-1,    1.7732942100e-1,
        -1.7483615880e-1,    1.4924806360e-1,   -1.0080891850e-1,    3.5677909850e-2,
    ];
    let expected_deviations = vec![
        2.7757203580e-1,    8.3271610740e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.1132486340e1,     5.2619037630e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.4338266850e-2,    1.3051480050e-1,    1.9485309720e-1,
        2.5919139390e-1,    3.2169145350e-1,    4.0000000600e-1,    4.5514711740e-1,
        4.9742656950e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,2,35
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n16_grid35() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 35);

    let expected_impulse_response = vec![
        -1.8361677230e-1,    1.2417693440e-1,   -1.5367421510e-1,    1.7335110900e-1,
        -1.7492020130e-1,    1.5214437250e-1,   -1.0411292310e-1,    3.7103295330e-2,
    ];
    let expected_deviations = vec![
        2.5909674170e-1,    7.7729022500e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.1730761530e1,     4.9951672550e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    6.4285717900e-2,    1.3035725060e-1,    1.9464282690e-1,
        2.5892847780e-1,    3.2142886520e-1,    4.0000000600e-1,    4.5357176660e-1,
        4.9642917510e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
41,1,2,15
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n41_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(41, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        1.7766407690e-3,   -7.0793586780e-3,    5.8723590340e-3,   -4.0236408820e-3,
        -2.1520084700e-3,    8.8611673560e-3,   -1.2099668380e-2,    8.1340121110e-3,
        2.9912954200e-3,   -1.6282893720e-2,    2.3380929600e-2,   -1.7039805650e-2,
        -3.6575859410e-3,    3.0411293730e-2,   -4.7613680360e-2,    3.8941107690e-2,
        4.0911613030e-3,   -7.5984261930e-2,    1.5727502110e-1,   -2.2141438720e-1,
        2.4574780460e-1,
    ];
    let expected_deviations = vec![
        5.4767699910e-3,    1.6430309040e-2,
    ];
    let expected_deviation_dbs = vec![
        -4.5229511260e1,    1.4155170320e-1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.6984127240e-2,    5.3968254480e-2,    8.0952383580e-2,
        1.0634920750e-1,    1.3333334030e-1,    1.5873016420e-1,    1.8412698810e-1,
        2.0952381190e-1,    2.3492063580e-1,    2.6031747460e-1,    2.8412699700e-1,
        3.0634921790e-1,    3.2698413730e-1,    3.4285715220e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.0793651340e-1,    4.2857143280e-1,    4.5079365370e-1,
        4.7460317610e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
41,1,2,16
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n41_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(41, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        1.7687934450e-3,   -7.0819593970e-3,    5.8801481500e-3,   -4.0319291870e-3,
        -2.1434004880e-3,    8.8539095600e-3,   -1.2096347290e-2,    8.1376768650e-3,
        2.9817002360e-3,   -1.6275763510e-2,    2.3377316070e-2,   -1.7039319500e-2,
        -3.6498878620e-3,    3.0407220130e-2,   -4.7612585130e-2,    3.8942590360e-2,
        4.0835081600e-3,   -7.5980894270e-2,    1.5727521480e-1,   -2.2141908110e-1,
        2.4575860800e-1,
    ];
    let expected_deviations = vec![
        5.4875644860e-3,    1.6462692990e-2,
    ];
    let expected_deviation_dbs = vec![
        -4.5212406160e1,    1.4182880520e-1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.6785714550e-2,    5.3571455180e-2,    8.0357201400e-2,
        1.0714294760e-1,    1.3392864170e-1,    1.5922616420e-1,    1.8452368680e-1,
        2.0982120930e-1,    2.3511873190e-1,    2.5892817970e-1,    2.8422570230e-1,
        3.0654704570e-1,    3.2738029960e-1,    3.4374928470e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.0892854330e-1,    4.2827370760e-1,    4.5059505110e-1,
        4.7589257360e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
41,1,2,34
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n41_grid34() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(41, FilterType::MultipleBand, &bands, 34);

    let expected_impulse_response = vec![
        1.7712987030e-3,   -7.0835691880e-3,    5.8771837500e-3,   -4.0313885550e-3,
        -2.1452782680e-3,    8.8549712670e-3,   -1.2097217140e-2,    8.1342663620e-3,
        2.9864397370e-3,   -1.6276886690e-2,    2.3374643180e-2,   -1.7035700380e-2,
        -3.6525325850e-3,    3.0406367030e-2,   -4.7611732040e-2,    3.8939584050e-2,
        4.0900367310e-3,   -7.5984559950e-2,    1.5727692840e-1,   -2.2141800820e-1,
        2.4575537440e-1,
    ];
    let expected_deviations = vec![
        5.4949484770e-3,    1.6484845430e-2,
    ];
    let expected_deviation_dbs = vec![
        -4.5200729370e1,    1.4201827350e-1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.7310924600e-2,    5.3921569140e-2,    8.0532215540e-2,
        1.0714285820e-1,    1.3305322830e-1,    1.5896359090e-1,    1.8487395350e-1,
        2.1008403600e-1,    2.3529411850e-1,    2.5980412960e-1,    2.8361415860e-1,
        3.0602359770e-1,    3.2703244690e-1,    3.4313923120e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.0840354560e-1,    4.2801180480e-1,    4.5112153890e-1,
        4.7493156790e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
41,1,2,35
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n41_grid35() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(41, FilterType::MultipleBand, &bands, 35);

    let expected_impulse_response = vec![
        1.7712883420e-3,   -7.0832683700e-3,    5.8761890980e-3,   -4.0302653800e-3,
        -2.1462561560e-3,    8.8554127140e-3,   -1.2096938680e-2,    8.1342235210e-3,
        2.9878471980e-3,   -1.6277534890e-2,    2.3375069720e-2,   -1.7035653810e-2,
        -3.6536476110e-3,    3.0406724660e-2,   -4.7611672430e-2,    3.8939110930e-2,
        4.0908893570e-3,   -7.5985297560e-2,    1.5727691350e-1,   -2.2141808270e-1,
        2.4575527010e-1,
    ];
    let expected_deviations = vec![
        5.4946490560e-3,    1.6483947630e-2,
    ];
    let expected_deviation_dbs = vec![
        -4.5201202390e1,    1.4201012250e-1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.7210878210e-2,    5.3741469980e-2,    8.0952435730e-2,
        1.0680289570e-1,    1.3333353400e-1,    1.5918371080e-1,    1.8503388760e-1,
        2.1020379660e-1,    2.3537370560e-1,    2.5986334680e-1,    2.8367272020e-1,
        3.0612155800e-1,    3.2720986010e-1,    3.4353628750e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.0884348750e-1,    4.2789098620e-1,    4.5102009180e-1,
        4.7550973300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
42,1,2,15
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n42_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(42, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        1.8928360940e-1,   -1.4494249220e-1,    1.3349333410e-1,   -7.6683938500e-2,
        -2.1567940710e-3,    6.4673364160e-2,   -7.7267110350e-2,    3.0070096250e-2,
        5.3290486340e-2,   -1.2796947360e-1,    1.5022617580e-1,   -1.0438263420e-1,
        1.4008045200e-2,    6.8011462690e-2,   -8.7534666060e-2,    2.1032035350e-2,
        1.0818630460e-1,   -2.3787617680e-1,    2.9783159490e-1,   -2.4671649930e-1,
        9.5542192460e-2,
    ];
    let expected_deviations = vec![
        2.4023792150e-1,    7.2071379420e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.2387168880e1,     4.7141733170e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.3809524250e-2,    4.7619048510e-2,    7.3015876110e-2,
        9.6825398500e-2,    1.2063492090e-1,    1.4444445070e-1,    1.6825397310e-1,
        1.9206349550e-1,    2.1587301790e-1,    2.3968254030e-1,    2.6349207760e-1,
        2.8730160000e-1,    3.0952382090e-1,    3.3174604180e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.1428571940e-1,    4.3809524180e-1,    4.6190476420e-1,
        4.8412698510e-1,    4.9841269850e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
42,1,2,16
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n42_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(42, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        1.8219369650e-1,   -1.4004945760e-1,    1.2913438680e-1,   -7.4437290430e-2,
        -1.6081333160e-3,    6.1956048010e-2,   -7.3886394500e-2,    2.8007090090e-2,
        5.2841663360e-2,   -1.2511688470e-1,    1.4646703000e-1,   -1.0181277990e-1,
        1.4116048810e-2,    6.5088689330e-2,   -8.3260476590e-2,    1.7597496510e-2,
        1.0891509060e-1,   -2.3544329400e-1,    2.9344397780e-1,   -2.4261987210e-1,
        9.3898177150e-2,
    ];
    let expected_deviations = vec![
        2.3084990680e-1,    6.9254970550e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.2733406070e1,     4.5708289150e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.3809524250e-2,    4.7619067130e-2,    7.2916716340e-2,
        9.6726268530e-2,    1.2053582070e-1,    1.4434526860e-1,    1.6815470160e-1,
        1.9196413460e-1,    2.1577356760e-1,    2.3958300050e-1,    2.6339244840e-1,
        2.8720188140e-1,    3.0952322480e-1,    3.3184456830e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.1488090160e-1,    4.3720224500e-1,    4.6101167800e-1,
        4.8482111100e-1,    4.9821391700e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
42,1,2,34
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n42_grid34() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(42, FilterType::MultipleBand, &bands, 34);

    let expected_impulse_response = vec![
        2.0156918470e-1,   -1.5355195110e-1,    1.4093506340e-1,   -8.0526947980e-2,
        -3.1912326810e-3,    6.9627761840e-2,   -8.3039283750e-2,    3.3564001320e-2,
        5.4517626760e-2,   -1.3340798020e-1,    1.5731322770e-1,   -1.0940128560e-1,
        1.4287233350e-2,    7.2737157340e-2,   -9.4931006430e-2,    2.7045905590e-2,
        1.0675120350e-1,   -2.4198836090e-1,    3.0536520480e-1,   -2.5373768810e-1,
        9.8386168480e-2,
    ];
    let expected_deviations = vec![
        2.5664734840e-1,    7.6994204520e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.1813264850e1,     4.9591813090e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.3809524250e-2,    4.8319328580e-2,    7.2128854690e-2,
        9.6638657150e-2,    1.2044817950e-1,    1.4425770940e-1,    1.6876751180e-1,
        1.9257703420e-1,    2.1638655660e-1,    2.4019607900e-1,    2.6330560450e-1,
        2.8711563350e-1,    3.0952507260e-1,    3.3193451170e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.1470620040e-1,    4.3781593440e-1,    4.6162596340e-1,
        4.8473569750e-1,    4.9874159690e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
42,1,2,35
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n42_grid35() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(42, FilterType::MultipleBand, &bands, 35);

    let expected_impulse_response = vec![
        2.2585810720e-1,   -1.7043720190e-1,    1.5581834320e-1,   -8.8215589520e-2,
        -5.0718784330e-3,    7.9134821890e-2,   -9.4548136000e-2,    4.0413737300e-2,
        5.6601613760e-2,   -1.4388871190e-1,    1.7097777130e-1,   -1.1905103920e-1,
        1.4550745490e-2,    8.2338213920e-2,   -1.0954290630e-1,    3.8982570170e-2,
        1.0399746890e-1,   -2.5012588500e-1,    3.2027208810e-1,   -2.6766562460e-1,
        1.0401380060e-1,
    ];
    let expected_deviations = vec![
        2.8882554170e-1,    8.6647665500e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.0787287710e1,     5.4204516410e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.3809520530e-2,    4.8299297690e-2,    7.2108857330e-2,
        9.6598766740e-2,    1.2040840090e-1,    1.4421781900e-1,    1.6870746020e-1,
        1.9251683350e-1,    2.1632620690e-1,    2.4013558030e-1,    2.6394495370e-1,
        2.8707405920e-1,    3.0952289700e-1,    3.3197173480e-1,    3.4999999400e-1,
        4.0000000600e-1,    4.1428563000e-1,    4.3741473560e-1,    4.6190437670e-1,
        4.8503348230e-1,    4.9931910630e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
77,1,2,15
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n77_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(77, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        -1.8132856350e-4,    7.7392796810e-6,    2.3316335860e-4,   -5.3586676950e-4,
        6.6217791750e-4,   -3.7932884880e-4,   -3.3392160550e-4,    1.1562615400e-3,
        -1.5153929130e-3,    9.2738907550e-4,    5.8624963280e-4,   -2.3030368610e-3,
        3.0640277550e-3,   -1.9658578090e-3,   -8.8664318900e-4,    4.1040279900e-3,
        -5.5735940110e-3,    3.7189738360e-3,    1.2282244860e-3,   -6.8393750120e-3,
        9.5097478480e-3,   -6.5823188050e-3,   -1.5810369510e-3,    1.1001836510e-2,
        -1.5746673570e-2,    1.1335324500e-2,    1.9099896310e-3,   -1.7756735910e-2,
        2.6465246450e-2,   -2.0046371970e-2,   -2.1785087880e-3,    3.0972121280e-2,
        -4.9672406170e-2,    4.1367132220e-2,    2.3541406260e-3,   -7.5490236280e-2,
        1.5800082680e-1,   -2.2296504680e-1,    2.4758480490e-1,
    ];
    let expected_deviations = vec![
        2.7334870540e-4,    8.2004605790e-4,
    ];
    let expected_deviation_dbs = vec![
        -7.1265655520e1,    7.1198642250e-3,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.2820512990e-2,    2.6495726780e-2,    3.9316240700e-2,
        5.2136752750e-2,    6.5811969340e-2,    7.8632481400e-2,    9.1452993450e-2,
        1.0427350550e-1,    1.1794871840e-1,    1.3076923790e-1,    1.4358974990e-1,
        1.5641026200e-1,    1.7008547480e-1,    1.8290598690e-1,    1.9572649900e-1,
        2.0854701100e-1,    2.2136752310e-1,    2.3418803510e-1,    2.4786324800e-1,
        2.5982907410e-1,    2.7264958620e-1,    2.8547009830e-1,    2.9743590950e-1,
        3.0940172080e-1,    3.2136753200e-1,    3.3162394170e-1,    3.4102565050e-1,
        3.4786325690e-1,    3.4999999400e-1,    4.0000000600e-1,    4.0256410840e-1,
        4.1025641560e-1,    4.2136752610e-1,    4.3333333730e-1,    4.4529914860e-1,
        4.5897436140e-1,    4.7264957430e-1,    4.8632478710e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
77,1,2,16
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n77_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(77, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -1.8205103700e-4,    8.2603346530e-6,    2.3302176850e-4,   -5.3663458680e-4,
        6.6396786130e-4,   -3.8147351010e-4,   -3.3247572720e-4,    1.1562599100e-3,
        -1.5170135300e-3,    9.3002163340e-4,    5.8399129190e-4,   -2.3028869180e-3,
        3.0665404630e-3,   -1.9697821700e-3,   -8.8333862370e-4,    4.1032512670e-3,
        -5.5763162670e-3,    3.7241957620e-3,    1.2233325980e-3,   -6.8374918770e-3,
        9.5121636990e-3,   -6.5881623890e-3,   -1.5750892930e-3,    1.0999009940e-2,
        -1.5748882670e-2,    1.1341707780e-2,    1.9028905080e-3,   -1.7752900720e-2,
        2.6467168700e-2,   -2.0053377380e-2,   -2.1701401570e-3,    3.0967129390e-2,
        -4.9673628060e-2,    4.1374202820e-2,    2.3450034200e-3,   -7.5484082100e-2,
        1.5800103550e-1,   -2.2297158840e-1,    2.4759396910e-1,
    ];
    let expected_deviations = vec![
        2.7430983030e-4,    8.2292954900e-4,
    ];
    let expected_deviation_dbs = vec![
        -7.1235168460e1,    7.1446942170e-3,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.2820512990e-2,    2.6442307980e-2,    3.9262838660e-2,
        5.2083380520e-2,    6.5705187620e-2,    7.8525669870e-2,    9.1346152130e-2,
        1.0496791450e-1,    1.1778839680e-1,    1.3060888650e-1,    1.4342936870e-1,
        1.5705113110e-1,    1.6987161340e-1,    1.8269209560e-1,    1.9551257790e-1,
        2.0913434030e-1,    2.2195482250e-1,    2.3477530480e-1,    2.4759578700e-1,
        2.6041644810e-1,    2.7323716880e-1,    2.8525659440e-1,    2.9807731510e-1,
        3.1009674070e-1,    3.2131487130e-1,    3.3173170690e-1,    3.4134724740e-1,
        3.4775760770e-1,    3.4999999400e-1,    4.0000000600e-1,    4.0320518610e-1,
        4.1041684150e-1,    4.2083367710e-1,    4.3285310270e-1,    4.4567382340e-1,
        4.5849454400e-1,    4.7211655970e-1,    4.8653987050e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
78,1,2,15
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n78_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(78, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        2.2772736850e-1,   -2.3223893340e-1,    1.9673305750e-1,   -4.2791008950e-2,
        -1.6546028850e-1,    3.0281209950e-1,   -2.5848993660e-1,    1.9136914980e-2,
        2.9634165760e-1,   -4.8916351800e-1,    4.0270328520e-1,   -4.0457606320e-2,
        -4.0539455410e-1,    6.5167069440e-1,   -5.0191593170e-1,   -6.0852169990e-3,
        5.9002256390e-1,   -8.7780737880e-1,    6.4454543590e-1,    1.9121170040e-2,
        -7.3000055550e-1,     1.0331826210e0,   -6.9379580020e-1,   -1.2807357310e-1,
        9.4675409790e-1,    -1.2402226920e0,    7.8603160380e-1,    1.7390799520e-1,
        -1.0597387550e0,     1.3067837950e0,   -7.2782224420e-1,   -3.5056990390e-1,
        1.2761294840e0,    -1.4745154380e0,    8.1514525410e-1,    2.9765343670e-1,
        -1.1733530760e0,     1.2702544930e0,   -5.4005038740e-1,
    ];
    let expected_deviations = vec![
        2.3741826420e-1,    7.1225482230e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.2489717480e1,     4.6713681220e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.2820512990e-2,    2.5641025980e-2,    3.8461539890e-2,
        5.1282051950e-2,    6.4102567730e-2,    7.6923079790e-2,    8.9743591850e-2,
        1.0256410390e-1,    1.1538461600e-1,    1.2820513550e-1,    1.4102564750e-1,
        1.5384615960e-1,    1.6666667160e-1,    1.7948718370e-1,    1.9230769570e-1,
        2.0512820780e-1,    2.1794871990e-1,    2.3076923190e-1,    2.4358974400e-1,
        2.5641027090e-1,    2.6837608220e-1,    2.8119659420e-1,    2.9316240550e-1,
        3.0512821670e-1,    3.1709402800e-1,    3.2820513840e-1,    3.3846154810e-1,
        3.4700855610e-1,    3.4999999400e-1,    4.0000000600e-1,    4.0683761240e-1,
        4.1794872280e-1,    4.2991453410e-1,    4.4188034530e-1,    4.5470085740e-1,
        4.6752136950e-1,    4.7948718070e-1,    4.9145299200e-1,    4.9914529920e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
78,1,2,16
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n78_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(78, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        1.8322560190e-1,   -1.8751889470e-1,    1.5999081730e-1,   -3.7318646910e-2,
        -1.2959152460e-1,    2.4059791860e-1,   -2.0655815300e-1,    1.5912067150e-2,
        2.3693759740e-1,   -3.9301788810e-1,    3.2622903590e-1,   -3.8048982620e-2,
        -3.1888633970e-1,    5.1817947630e-1,   -4.0157699580e-1,   -2.5684237480e-3,
        4.7003537420e-1,   -7.0346480610e-1,    5.2151632310e-1,    5.8442354200e-3,
        -5.7430380580e-1,    8.2041484120e-1,   -5.5446088310e-1,   -9.8772674800e-2,
        7.5365495680e-1,   -9.9326843020e-1,    6.3744759560e-1,    1.2440121170e-1,
        -8.3211731910e-1,     1.0348742010e0,   -5.7998752590e-1,   -2.7679979800e-1,
        1.0188038350e0,    -1.1881078480e0,    6.7673587800e-1,    1.9922614100e-1,
        -8.9514970780e-1,    9.8189258580e-1,   -4.1901624200e-1,
    ];
    let expected_deviations = vec![
        1.9076956810e-1,    5.7230871920e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.4389818190e1,     3.9307570460e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.2820512990e-2,    2.5641025980e-2,    3.8461554800e-2,
        5.1282096650e-2,    6.4102627340e-2,    7.6923109590e-2,    8.9743591850e-2,
        1.0256407410e-1,    1.1538455640e-1,    1.2820504610e-1,    1.4102552830e-1,
        1.5384601060e-1,    1.6666649280e-1,    1.7948697510e-1,    1.9230745730e-1,
        2.0512793960e-1,    2.1794842180e-1,    2.3076890410e-1,    2.4358938630e-1,
        2.5640997290e-1,    2.6842939850e-1,    2.8125011920e-1,    2.9326954480e-1,
        3.0528897050e-1,    3.1730839610e-1,    3.2852652670e-1,    3.3894336220e-1,
        3.4695631270e-1,    3.4999999400e-1,    4.0000000600e-1,    4.0641036630e-1,
        4.1762849690e-1,    4.2964792250e-1,    4.4246864320e-1,    4.5448806880e-1,
        4.6730878950e-1,    4.8012951020e-1,    4.9134764080e-1,    4.9855929610e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,1,2,15
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n128_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        -3.3191770320e-1,    5.1927131410e-1,   -5.2128803730e-1,    1.0198867320e-1,
        6.8047648670e-1,    -1.4049494270e0,     1.4464393850e0,   -3.9692497250e-1,
        -1.4823305610e0,     3.1493253710e0,    -3.2545623780e0,     1.0592837330e0,
        2.7498390670e0,    -6.0296497340e0,     6.2439303400e0,    -2.1914014820e0,
        -4.6760978700e0,     1.0455475810e1,    -1.0832147600e1,     4.0110836030e0,
        7.3229608540e0,    -1.6678562160e1,     1.7269218440e1,    -6.5990791320e0,
        -1.0823745730e1,     2.4959840770e1,    -2.5805498120e1,     1.0107135770e1,
        1.5117004390e1,    -3.5260356900e1,     3.6377487180e1,    -1.4472938540e1,
        -2.0190822600e1,     4.7460372920e1,    -4.8834278110e1,     1.9681095120e1,
        2.5783851620e1,    -6.1038837430e1,     6.2613239290e1,    -2.5441535950e1,
        -3.1709892270e1,     7.5411743160e1,    -7.7095207210e1,     3.1549278260e1,
        3.7511108400e1,    -8.9600494380e1,     9.1265769960e1,    -3.7502601620e1,
        -4.2909126280e1,     1.0274861910e2,    -1.0425902560e2,     4.3021720890e1,
        4.7342048650e1,    -1.1369538120e2,     1.1490550230e2,    -4.7500854490e1,
        -5.0643707280e1,     1.2176974490e2,    -1.2258646390e2,     5.0869697570e1,
        5.2185508730e1,    -1.2588601680e2,     1.2619613650e2,    -5.2348228450e1,
    ];
    let expected_deviations = vec![
        2.3546537760e-1,    7.0639616250e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.2561458590e1,     4.6415967940e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    7.8125000000e-3,    1.5625005590e-2,    2.3437498140e-2,
        3.1249990690e-2,    3.9062485100e-2,    4.6874977650e-2,    5.4687470200e-2,
        6.2499962750e-2,    7.0312455300e-2,    7.8124947850e-2,    8.5937440400e-2,
        9.3749932940e-2,    1.0156242550e-1,    1.0937491800e-1,    1.1718741060e-1,
        1.2447907030e-1,    1.3229165970e-1,    1.4010426400e-1,    1.4791686830e-1,
        1.5572947260e-1,    1.6354207690e-1,    1.7135468130e-1,    1.7916728560e-1,
        1.8697988990e-1,    1.9479249420e-1,    2.0208425820e-1,    2.0989686250e-1,
        2.1770946680e-1,    2.2552207110e-1,    2.3333467540e-1,    2.4062643950e-1,
        2.4843904380e-1,    2.5625145440e-1,    2.6354300980e-1,    2.7135539050e-1,
        2.7864694600e-1,    2.8645932670e-1,    2.9375088210e-1,    3.0104243760e-1,
        3.0833399300e-1,    3.1562554840e-1,    3.2239627840e-1,    3.2916700840e-1,
        3.3541691300e-1,    3.4114599230e-1,    3.4583342080e-1,    3.4895837310e-1,
        3.4999999400e-1,    4.0000000600e-1,    4.0364578370e-1,    4.0885403750e-1,
        4.1510394220e-1,    4.2135384680e-1,    4.2864540220e-1,    4.3541613220e-1,
        4.4322851300e-1,    4.5052006840e-1,    4.5781162380e-1,    4.6562400460e-1,
        4.7291556000e-1,    4.8072794080e-1,    4.8801949620e-1,    4.9531105160e-1,
        4.9947765470e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,1,2,16
0.000000000000,0.349999994040,0.400000005960,0.500000000000
0.000000000000,1.000000000000
3.000000000000,1.000000000000
*/
#[test]
fn multipleband_highpass_n128_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 3f32,
    });
    bands.push(Band {
        lower_edge: 0.4f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -2.6572597030e-1,    4.1662013530e-1,   -4.1946554180e-1,    8.4764003750e-2,
        5.4175263640e-1,    -1.1234672070e0,     1.1594886780e0,   -3.2150995730e-1,
        -1.1832180020e0,     2.5214650630e0,    -2.6115064620e0,    8.5845756530e-1,
        2.1915423870e0,    -4.8240852360e0,     5.0059180260e0,    -1.7691044810e0,
        -3.7299537660e0,     8.3680591580e0,    -8.6860504150e0,     3.2370595930e0,
        5.8384456630e0,    -1.3345220570e1,     1.3841941830e1,    -5.3161640170e0,
        -8.6350545880e0,     1.9975784300e1,    -2.0685215000e1,     8.1390228270e0,
        1.2059251790e1,    -2.8216377260e1,     2.9151599880e1,    -1.1641567230e1,
        -1.6115837100e1,     3.7985321040e1,    -3.9133956910e1,     1.5824512480e1,
        2.0582244870e1,    -4.8850593570e1,     5.0165187840e1,    -2.0437843320e1,
        -2.5327636720e1,     6.0362911220e1,    -6.1767688750e1,     2.5335693360e1,
        2.9966247560e1,    -7.1717880250e1,     7.3106674190e1,    -3.0092491150e1,
        -3.4300086980e1,     8.2255645750e1,    -8.3515769960e1,     3.4512504580e1,
        3.7848945620e1,    -9.1014320370e1,     9.2024566650e1,    -3.8075119020e1,
        -4.0518249510e1,     9.7500083920e1,    -9.8187728880e1,     4.0782722470e1,
        4.1735420230e1,    -1.0076264190e2,     1.0102805330e2,    -4.1910453800e1,
    ];
    let expected_deviations = vec![
        1.8818597500e-1,    5.6455796960e-1,
    ];
    let expected_deviation_dbs = vec![
        -1.4508255000e1,     3.8878335950e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    7.8125000000e-3,    1.5625000000e-2,    2.3437500000e-2,
        3.1250000000e-2,    3.9062500000e-2,    4.6875000000e-2,    5.4687500000e-2,
        6.2500000000e-2,    7.0312500000e-2,    7.8125000000e-2,    8.5937500000e-2,
        9.3750000000e-2,    1.0156250000e-1,    1.0937500000e-1,    1.1669921880e-1,
        1.2451171880e-1,    1.3232421880e-1,    1.4013671880e-1,    1.4794921880e-1,
        1.5576171880e-1,    1.6357421880e-1,    1.7138671880e-1,    1.7919921880e-1,
        1.8701171880e-1,    1.9482421880e-1,    2.0214843750e-1,    2.0996093750e-1,
        2.1777343750e-1,    2.2558593750e-1,    2.3339843750e-1,    2.4072265620e-1,
        2.4853515620e-1,    2.5634765620e-1,    2.6367187500e-1,    2.7148437500e-1,
        2.7880859380e-1,    2.8613281250e-1,    2.9394531250e-1,    3.0126953120e-1,
        3.0859375000e-1,    3.1542968750e-1,    3.2226562500e-1,    3.2910156250e-1,
        3.3544921880e-1,    3.4082031250e-1,    3.4570312500e-1,    3.4863281250e-1,
        3.4999999400e-1,    4.0000000600e-1,    4.0390625600e-1,    4.0878906850e-1,
        4.1464844350e-1,    4.2148438100e-1,    4.2832031850e-1,    4.3564453720e-1,
        4.4296875600e-1,    4.5029297470e-1,    4.5810547470e-1,    4.6542969350e-1,
        4.7275391220e-1,    4.8056641220e-1,    4.8789063100e-1,    4.9472656850e-1,
        4.9912109970e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,1,3,15
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n11_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        -2.6046151300e-2,    3.7154424940e-2,   -4.9197159710e-2,    4.0683308240e-1,
        7.5243271890e-2,    6.2077325580e-1,
    ];
    let expected_deviations = vec![
        5.0874817370e-1,    1.2718704340e-1,    5.0874817370e-1,
    ];
    let expected_deviation_dbs = vec![
        3.5723352430e0,    -1.7911144260e1,     3.5723352430e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5000000600e-1,    1.8999999760e-1,    2.6777783040e-1,
        3.4999999400e-1,    3.8999998570e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,1,3,16
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n11_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -2.6042437180e-2,    3.7142094220e-2,   -4.9202147870e-2,    4.0684083100e-1,
        7.5244553390e-2,    6.2076473240e-1,
    ];
    let expected_deviations = vec![
        5.0873059030e-1,    1.2718264760e-1,    5.0873059030e-1,
    ];
    let expected_deviation_dbs = vec![
        3.5722346310e0,    -1.7911441800e1,     3.5722346310e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5000000600e-1,    1.8999999760e-1,    2.6812496780e-1,
        3.4999999400e-1,    3.8999998570e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,1,3,36
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n11_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::MultipleBand, &bands, 36);

    let expected_impulse_response = vec![
        -2.6055775580e-2,    3.7186365570e-2,   -4.9184218050e-2,    4.0681287650e-1,
        7.5239963830e-2,    6.2079536910e-1,
    ];
    let expected_deviations = vec![
        5.0879377130e-1,    1.2719844280e-1,    5.0879377130e-1,
    ];
    let expected_deviation_dbs = vec![
        3.5725979800e0,    -1.7910364150e1,     3.5725979800e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5000000600e-1,    1.8999999760e-1,    2.6638898250e-1,
        3.4999999400e-1,    3.8999998570e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,1,3,37
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n11_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::MultipleBand, &bands, 37);

    let expected_impulse_response = vec![
        -2.6054967190e-2,    3.7183687090e-2,   -4.9185309560e-2,    4.0681457520e-1,
        7.5240239500e-2,    6.2079352140e-1,
    ];
    let expected_deviations = vec![
        5.0878995660e-1,    1.2719748910e-1,    5.0878995660e-1,
    ];
    let expected_deviation_dbs = vec![
        3.5725760460e0,    -1.7910430910e1,     3.5725760460e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.5000000600e-1,    1.8999999760e-1,    2.6657652850e-1,
        3.4999999400e-1,    3.8999998570e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,1,3,15
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n12_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        -2.6326033470e-1,    4.3608725070e-2,   -2.8398296240e-1,    2.6139622930e-1,
        7.4695885180e-2,    3.7651741500e-1,
    ];
    let expected_deviations = vec![
        5.8204966780e-1,    1.4551241700e-1,    5.8204966780e-1,
    ];
    let expected_deviation_dbs = vec![
        3.9844021800e0,    -1.6741998670e1,     3.9844021800e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    8.8888891040e-2,    2.0666666330e-1,    2.6777783040e-1,
        3.3444467190e-1,    4.4000011680e-1,    4.9000024800e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,1,3,16
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n12_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -3.0497267840e-1,    3.7094950680e-2,   -3.2000121470e-1,    2.6505649090e-1,
        6.8599343300e-2,    3.8993895050e-1,
    ];
    let expected_deviations = vec![
        7.2856849430e-1,    1.8214212360e-1,    7.2856849430e-1,
    ];
    let expected_deviation_dbs = vec![
        4.7537322040e0,    -1.4791791920e1,     4.7537322040e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    9.3750007450e-2,    2.0562498270e-1,    2.6812496780e-1,
        3.3583343030e-1,    4.4208341840e-1,    4.9416685100e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,1,3,36
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n12_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::MultipleBand, &bands, 36);

    let expected_impulse_response = vec![
        -3.3218932150e-1,    3.2383322720e-2,   -3.4386533500e-1,    2.6723337170e-1,
        6.4659297470e-2,    3.9890527730e-1,
    ];
    let expected_deviations = vec![
        8.2574665550e-1,    2.0643666390e-1,    8.2574665550e-1,
    ];
    let expected_deviation_dbs = vec![
        5.2288103100e0,    -1.3704262730e1,     5.2288103100e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    9.4907373190e-2,    2.0388892290e-1,    2.7101859450e-1,
        3.3583316210e-1,    4.4092571740e-1,    4.9648106100e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,1,3,37
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n12_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::MultipleBand, &bands, 37);

    let expected_impulse_response = vec![
        -3.2448676230e-1,    3.3726304770e-2,   -3.3716833590e-1,    2.6664501430e-1,
        6.5785169600e-2,    3.9637690780e-1,
    ];
    let expected_deviations = vec![
        7.9824364190e-1,    1.9956091050e-1,    7.9824364190e-1,
    ];
    let expected_deviation_dbs = vec![
        5.0969705580e0,    -1.3998490330e1,     5.0969705580e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    9.4594560560e-2,    2.0351350310e-1,    2.7108103040e-1,
        3.3639630680e-1,    4.4180175660e-1,    4.9585577850e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,3,15
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n16_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        -1.5341199930e-1,    1.9037455320e-2,   -2.4213476480e-1,    6.0381025080e-2,
        -2.3996725680e-1,    2.6327878240e-1,    9.1471850870e-2,    3.6227935550e-1,
    ];
    let expected_deviations = vec![
        6.7813104390e-1,    1.6953276100e-1,    6.7813104390e-1,
    ];
    let expected_deviation_dbs = vec![
        4.4965181350e0,    -1.5414926530e1,     4.4965181350e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    7.5000010430e-2,    1.5000000600e-1,    1.8999999760e-1,
        2.6083326340e-1,    3.2749986650e-1,    3.8999998570e-1,    4.5249992610e-1,
        4.9416655300e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,3,16
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n16_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -1.7485190930e-1,    1.6738787290e-2,   -2.5495487450e-1,    6.5705090760e-2,
        -2.3779523370e-1,    2.6294434070e-1,    9.3916893010e-2,    3.6005747320e-1,
    ];
    let expected_deviations = vec![
        7.3647892480e-1,    1.8411973120e-1,    7.3647892480e-1,
    ];
    let expected_deviation_dbs = vec![
        4.7933902740e0,    -1.4697992320e1,     4.7933902740e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    7.4218750000e-2,    1.5000000600e-1,    1.9390624760e-1,
        2.6421874760e-1,    3.2671874760e-1,    3.8999998570e-1,    4.5249998570e-1,
        4.9546873570e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,3,36
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n16_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 36);

    let expected_impulse_response = vec![
        -2.1414285900e-1,    1.3535514470e-2,   -2.8171765800e-1,    7.4112415310e-2,
        -2.3335075380e-1,    2.6382428410e-1,    9.8952949050e-2,    3.5410338640e-1,
    ];
    let expected_deviations = vec![
        8.4936517480e-1,    2.1234129370e-1,    8.4936517480e-1,
    ];
    let expected_deviation_dbs = vec![
        5.3404541020e0,    -1.3459310530e1,     5.3404541020e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    7.4652805920e-2,    1.5000000600e-1,    2.0041662450e-1,
        2.6465249060e-1,    3.2541614770e-1,    3.8999998570e-1,    4.5249974730e-1,
        4.9763846400e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,1,3,37
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n16_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::MultipleBand, &bands, 37);

    let expected_impulse_response = vec![
        -2.2351484000e-1,    1.2732028960e-2,   -2.8832548860e-1,    7.5882017610e-2,
        -2.3218214510e-1,    2.6411110160e-1,    1.0031092170e-1,    3.5256278510e-1,
    ];
    let expected_deviations = vec![
        8.7684720750e-1,    2.1921180190e-1,    8.7684720750e-1,
    ];
    let expected_deviation_dbs = vec![
        5.4685788150e0,    -1.3182721140e1,     5.4685788150e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    7.4324309830e-2,    1.5000000600e-1,    2.0013517140e-1,
        2.6432460550e-1,    3.2513564830e-1,    3.8999998570e-1,    4.5418941970e-1,
        4.9810850620e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,1,3,15
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n47_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        2.4451230270e-5,   -5.3538070060e-3,    7.0190364490e-3,    1.0136167520e-2,
        3.0474893280e-3,    1.3121956030e-2,   -1.8595153470e-2,   -5.8298986410e-3,
        4.5848037120e-3,   -1.3994422040e-3,    2.9408397150e-2,   -1.9900503100e-3,
        -3.1889941540e-2,    2.6211633810e-3,   -2.6121933010e-2,    3.4241504970e-2,
        6.1914958060e-2,   -4.7044470910e-2,   -5.1204115150e-3,   -7.3058202860e-2,
        -7.1925193070e-2,    2.9349106550e-1,    4.7653507440e-2,    5.9110128880e-1,
    ];
    let expected_deviations = vec![
        2.8973283250e-2,    7.2433208120e-3,    2.8973283250e-2,
    ];
    let expected_deviation_dbs = vec![
        2.4808245900e-1,    -4.2801246640e1,    2.4808245900e-1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.5000000370e-2,    5.0000000750e-2,    7.5000032780e-2,
        9.8611205820e-2,    1.2222237880e-1,    1.4166687430e-1,    1.5000000600e-1,
        1.8999999760e-1,    1.9555556770e-1,    2.1083338560e-1,    2.2888898850e-1,
        2.4972237650e-1,    2.7055555580e-1,    2.8999984260e-1,    3.1083300710e-1,
        3.2888841630e-1,    3.4416607020e-1,    3.4999999400e-1,    3.8999998570e-1,
        3.9833325150e-1,    4.1916641590e-1,    4.4277733560e-1,    4.6916601060e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,1,3,16
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n47_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        2.2198970330e-5,   -5.3572631440e-3,    7.0450939240e-3,    1.0133954700e-2,
        3.0647201930e-3,    1.3114202770e-2,   -1.8601108340e-2,   -5.8125462380e-3,
        4.5743468220e-3,   -1.3784305190e-3,    2.9399219900e-2,   -1.9891588020e-3,
        -3.1887955960e-2,    2.6258963630e-3,   -2.6117892940e-2,    3.4252330660e-2,
        6.1894204470e-2,   -4.7044459730e-2,   -5.1184133630e-3,   -7.3052190240e-2,
        -7.1912772950e-2,    2.9348620770e-1,    4.7638375310e-2,    5.9107840060e-1,
    ];
    let expected_deviations = vec![
        2.9035540300e-2,    7.2588850740e-3,    2.9035540300e-2,
    ];
    let expected_deviation_dbs = vec![
        2.4860772490e-1,    -4.2782600400e1,    2.4860772490e-1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.6041669770e-2,    5.0781235100e-2,    7.5520828370e-2,
        9.8958373070e-2,    1.2239591780e-1,    1.4062502980e-1,    1.5000000600e-1,
        1.8999999760e-1,    1.9651038940e-1,    2.1083325150e-1,    2.2906234860e-1,
        2.4989560250e-1,    2.6942700150e-1,    2.9026049380e-1,    3.1109398600e-1,
        3.2932329180e-1,    3.4364631770e-1,    3.4999999400e-1,    3.8999998570e-1,
        3.9911463860e-1,    4.1864603760e-1,    4.4208371640e-1,    4.6812558170e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,1,3,36
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n47_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::MultipleBand, &bands, 36);

    let expected_impulse_response = vec![
        8.4171979320e-6,   -5.3504640240e-3,    7.0263254460e-3,    1.0145148260e-2,
        3.0602144540e-3,    1.3125980270e-2,   -1.8592745070e-2,   -5.8164400980e-3,
        4.5817322100e-3,   -1.3905955710e-3,    2.9396625240e-2,   -1.9951693250e-3,
        -3.1884144990e-2,    2.6255664420e-3,   -2.6108806950e-2,    3.4249208870e-2,
        6.1898499730e-2,   -4.7045722600e-2,   -5.1203793850e-3,   -7.3046915230e-2,
        -7.1911975740e-2,    2.9349085690e-1,    4.7646246850e-2,    5.9108662610e-1,
    ];
    let expected_deviations = vec![
        2.9069544750e-2,    7.2673861870e-3,    2.9069544750e-2,
    ];
    let expected_deviation_dbs = vec![
        2.4889451270e-1,    -4.2772438050e1,    2.4889451270e-1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.5462953370e-2,    5.0347257410e-2,    7.5231499970e-2,
        9.8958261310e-2,    1.2152761970e-1,    1.4120347800e-1,    1.5000000600e-1,
        1.8999999760e-1,    1.9636571410e-1,    2.1083325150e-1,    2.2935169940e-1,
        2.4960625170e-1,    2.6986080410e-1,    2.9069405790e-1,    3.1036990880e-1,
        3.2888835670e-1,    3.4393459560e-1,    3.4999999400e-1,    3.8999998570e-1,
        3.9925920960e-1,    4.1835635900e-1,    4.4208312030e-1,    4.6870338920e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,1,3,37
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n47_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::MultipleBand, &bands, 37);

    let expected_impulse_response = vec![
        1.0419848880e-5,   -5.3486265240e-3,    7.0315673950e-3,    1.0144844650e-2,
        3.0595662540e-3,    1.3121888040e-2,   -1.8597554420e-2,   -5.8117434380e-3,
        4.5805787670e-3,   -1.3858567690e-3,    2.9394453390e-2,   -2.0004548130e-3,
        -3.1884912400e-2,    2.6246057820e-3,   -2.6104880500e-2,    3.4252874550e-2,
        6.1898764220e-2,   -4.7046929600e-2,   -5.1215523850e-3,   -7.3046691720e-2,
        -7.1909978990e-2,    2.9349115490e-1,    4.7643534840e-2,    5.9108513590e-1,
    ];
    let expected_deviations = vec![
        2.9075294730e-2,    7.2688236830e-3,    2.9075294730e-2,
    ];
    let expected_deviation_dbs = vec![
        2.4894279240e-1,    -4.2770721440e1,    2.4894279240e-1,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.5337828320e-2,    5.0675652920e-2,    7.5450412930e-2,
        9.9099047480e-2,    1.2162155660e-1,    1.4132896070e-1,    1.5000000600e-1,
        1.8999999760e-1,    1.9619376960e-1,    2.1083359420e-1,    2.2941491010e-1,
        2.4968543650e-1,    2.6995542650e-1,    2.9078847170e-1,    3.1049540640e-1,
        3.2907623050e-1,    3.4371566770e-1,    3.4999999400e-1,    3.8999998570e-1,
        3.9900887010e-1,    4.1871580480e-1,    4.4180107120e-1,    4.6882772450e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,1,3,15
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n48_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        -2.3900672790e-1,   -7.9028904440e-3,   -1.8139892820e-1,    6.6845268010e-2,
        3.2522261140e-2,    7.8402459620e-3,    2.6386559010e-2,   -5.6955009700e-2,
        -5.3780496120e-2,    9.0290695430e-2,   -4.3645858760e-2,    8.6673498150e-2,
        -5.1025986670e-2,   -6.3306868080e-2,    1.8000006680e-2,    1.6656696800e-2,
        2.0853519440e-2,    1.1339002850e-1,   -1.8475610020e-1,    8.1543624400e-2,
        -1.8337672950e-1,    1.9054102900e-1,    1.4667963980e-1,    3.2881581780e-1,
    ];
    let expected_deviations = vec![
        6.7623454330e-1,    1.6905863580e-1,    6.7623454330e-1,
    ];
    let expected_deviation_dbs = vec![
        4.4866952900e0,    -1.5439252850e1,     4.4866952900e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.2222222760e-2,    4.3055556710e-2,    6.5277785060e-2,
        8.6111173030e-2,    1.0833345350e-1,    1.2916684150e-1,    1.5000000600e-1,
        1.8999999760e-1,    1.9972224530e-1,    2.1777784820e-1,    2.3722234370e-1,
        2.5666677950e-1,    2.7611106630e-1,    2.9555535320e-1,    3.1499964000e-1,
        3.3444392680e-1,    3.4999999400e-1,    3.8999998570e-1,    4.0249988440e-1,
        4.2333304880e-1,    4.4555509090e-1,    4.6638825540e-1,    4.8583254220e-1,
        4.9833244090e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,1,3,16
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n48_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -2.2569988670e-1,   -7.5398385520e-3,   -1.7167864740e-1,    6.4475119110e-2,
        2.9725700620e-2,    9.2669725420e-3,    2.3153930900e-2,   -5.3637772800e-2,
        -5.1523000000e-2,    8.6167722940e-2,   -4.1289240120e-2,    8.4336221220e-2,
        -5.1408112050e-2,   -5.8705091480e-2,    1.4385402200e-2,    1.7235934730e-2,
        2.1759092810e-2,    1.0973000530e-1,   -1.7961567640e-1,    8.0238044260e-2,
        -1.8550741670e-1,    1.9348216060e-1,    1.4241409300e-1,    3.3112907410e-1,
    ];
    let expected_deviations = vec![
        6.3820970060e-1,    1.5955242510e-1,    6.3820970060e-1,
    ];
    let expected_deviation_dbs = vec![
        4.2873902320e0,    -1.5941932680e1,     4.2873902320e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.2135417910e-2,    4.2968742550e-2,    6.5104141830e-2,
        8.5937514900e-2,    1.0807297380e-1,    1.3020840290e-1,    1.5000000600e-1,
        1.8999999760e-1,    2.0041662450e-1,    2.1734364330e-1,    2.3687481880e-1,
        2.5640606880e-1,    2.7593746780e-1,    2.9546886680e-1,    3.1500026580e-1,
        3.3453166480e-1,    3.4999999400e-1,    3.8999998570e-1,    4.0302091840e-1,
        4.2385441060e-1,    4.4468790290e-1,    4.6552139520e-1,    4.8635488750e-1,
        4.9807372690e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,1,3,36
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n48_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::MultipleBand, &bands, 36);

    let expected_impulse_response = vec![
        -3.0207303170e-1,   -1.0563433170e-2,   -2.2694429760e-1,    7.7853918080e-2,
        4.5299530030e-2,    2.7682781220e-3,    4.0004312990e-2,   -7.1312040090e-2,
        -6.4608067270e-2,    1.0918414590e-1,   -5.3387403490e-2,    9.6068143840e-2,
        -4.8592805860e-2,   -8.5278093810e-2,    3.4366130830e-2,    1.4585375790e-2,
        1.7287552360e-2,    1.3041424750e-1,   -2.0772731300e-1,    8.6364746090e-2,
        -1.7327558990e-1,    1.7623353000e-1,    1.6593229770e-1,    3.1728243830e-1,
    ];
    let expected_deviations = vec![
        8.6023569110e-1,    2.1505892280e-1,    8.6023569110e-1,
    ];
    let expected_deviation_dbs = vec![
        5.3913598060e0,    -1.3348851200e1,     5.3913598060e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.1412029860e-2,    4.3402794750e-2,    6.4814873040e-2,
        8.6226828400e-2,    1.0821748520e-1,    1.2962944810e-1,    1.5000000600e-1,
        1.8999999760e-1,    1.9983792300e-1,    2.1777766940e-1,    2.3687481880e-1,
        2.5655066970e-1,    2.7622652050e-1,    2.9590237140e-1,    3.1557822230e-1,
        3.3409667020e-1,    3.4999999400e-1,    3.8999998570e-1,    4.0331012010e-1,
        4.2414337400e-1,    4.4497662780e-1,    4.6638858320e-1,    4.8664313550e-1,
        4.9937456850e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,1,3,37
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n48_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::MultipleBand, &bands, 37);

    let expected_impulse_response = vec![
        -2.9231458900e-1,   -1.0113865140e-2,   -2.1993938090e-1,    7.6180100440e-2,
        4.3294608590e-2,    3.6133527760e-3,    3.7885725500e-2,   -6.9016993050e-2,
        -6.3012421130e-2,    1.0633099080e-1,   -5.1973402500e-2,    9.4664871690e-2,
        -4.8993766310e-2,   -8.1908166410e-2,    3.1834721570e-2,    1.4858424660e-2,
        1.7880678180e-2,    1.2780356410e-1,   -2.0414555070e-1,    8.5605859760e-2,
        -1.7484581470e-1,    1.7838203910e-1,    1.6294920440e-1,    3.1903839110e-1,
    ];
    let expected_deviations = vec![
        8.3188307290e-1,    2.0797076820e-1,    8.3188307290e-1,
    ];
    let expected_deviation_dbs = vec![
        5.2579550740e0,    -1.3639954570e1,     5.2579550740e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    2.1396389230e-2,    4.3355837460e-2,    6.4752221110e-2,
        8.6711667480e-2,    1.0810805110e-1,    1.3006755710e-1,    1.5000000600e-1,
        1.8999999760e-1,    2.0013526080e-1,    2.1759043630e-1,    2.3673482240e-1,
        2.5644209980e-1,    2.7614903450e-1,    2.9585596920e-1,    3.1556290390e-1,
        3.3414372800e-1,    3.4999999400e-1,    3.8999998570e-1,    4.0351331230e-1,
        4.2378330230e-1,    4.4517940280e-1,    4.6657550330e-1,    4.8684549330e-1,
        4.9923270940e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,1,3,15
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n81_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        -9.3938491770e-5,    6.5128685670e-4,    6.5719475970e-5,    2.0225462500e-3,
        -3.3656068260e-4,   -1.1134655210e-3,   -1.9240927940e-4,   -2.1061846060e-3,
        2.5171344170e-3,    2.7211816050e-3,   -1.5879099960e-3,    7.9955667020e-4,
        -5.1837111820e-3,   -2.0040881350e-3,    7.2130030020e-3,    4.4213179540e-6,
        3.9888764730e-3,   -2.6905515700e-3,   -1.3463689950e-2,    5.8983620260e-3,
        3.1480498150e-3,    5.6942310180e-3,    1.2577843850e-2,   -2.0972052590e-2,
        -8.5144136100e-3,    4.5881569390e-3,   -3.0451968780e-3,    3.4992080180e-2,
        -2.6588069740e-3,   -3.2530568540e-2,    1.7276351570e-3,   -2.9406137760e-2,
        3.6095414310e-2,    6.3683755700e-2,   -4.6786412600e-2,   -3.5930201410e-3,
        -7.5792655350e-2,   -7.2696149350e-2,    2.9346635940e-1,    4.7490984200e-2,
        5.9371137620e-1,
    ];
    let expected_deviations = vec![
        2.8686993760e-3,    7.1717484390e-4,    2.8686993760e-3,
    ];
    let expected_deviation_dbs = vec![
        2.4881128220e-2,    -6.2887496950e1,    2.4881128220e-2,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.3008131650e-2,    2.6016252120e-2,    3.9024371650e-2,
        5.2032489330e-2,    6.5040610730e-2,    7.8048728410e-2,    9.0243838730e-2,
        1.0325195640e-1,    1.1544706670e-1,    1.2764218450e-1,    1.3821128010e-1,
        1.4634135370e-1,    1.5000000600e-1,    1.8999999760e-1,    1.9243901970e-1,
        1.9975608590e-1,    2.0951217410e-1,    2.2008126970e-1,    2.3227638010e-1,
        2.4447149040e-1,    2.5666660070e-1,    2.6967471840e-1,    2.8268283610e-1,
        2.9569095370e-1,    3.0788606410e-1,    3.1926816700e-1,    3.3065027000e-1,
        3.4040635820e-1,    3.4772342440e-1,    3.4999999400e-1,    3.8999998570e-1,
        3.9325201510e-1,    4.0138208870e-1,    4.1276419160e-1,    4.2414629460e-1,
        4.3634140490e-1,    4.4934952260e-1,    4.6154463290e-1,    4.7455275060e-1,
        4.8756086830e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,1,3,16
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n81_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -9.2299014800e-5,    6.5133016320e-4,    6.2498795160e-5,    2.0201888400e-3,
        -3.3470313060e-4,   -1.1140673890e-3,   -1.8929934590e-4,   -2.1069939250e-3,
        2.5120044590e-3,    2.7242098000e-3,   -1.5886195470e-3,    8.0318195980e-4,
        -5.1788659770e-3,   -2.0138474650e-3,    7.2121876290e-3,    3.9013266360e-6,
        3.9872596970e-3,   -2.6769647370e-3,   -1.3468120250e-2,    5.8910874650e-3,
        3.1501927880e-3,    5.6852526030e-3,    1.2590469790e-2,   -2.0962595940e-2,
        -8.5255671290e-3,    4.5902924610e-3,   -3.0574898700e-3,    3.4988988190e-2,
        -2.6380305640e-3,   -3.2534282650e-2,    1.7294136340e-3,   -2.9409393670e-2,
        3.6075919870e-2,    6.3697949050e-2,   -4.6778976920e-2,   -3.5943919790e-3,
        -7.5783386830e-2,   -7.2718098760e-2,    2.9346123340e-1,    4.7506347300e-2,
        5.9370839600e-1,
    ];
    let expected_deviations = vec![
        2.8641964310e-3,    7.1604910770e-4,    2.8641964310e-3,
    ];
    let expected_deviation_dbs = vec![
        2.4842927230e-2,    -6.2901145940e1,    2.4842927230e-2,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.2957318690e-2,    2.5914626200e-2,    3.8871932770e-2,
        5.1829237490e-2,    6.4786545930e-2,    7.7743850650e-2,    9.0701155360e-2,
        1.0289626570e-1,    1.1585357040e-1,    1.2728649380e-1,    1.3795721530e-1,
        1.4634135370e-1,    1.5000000600e-1,    1.8999999760e-1,    1.9228658080e-1,
        1.9990852480e-1,    2.0905485750e-1,    2.2048777340e-1,    2.3192068930e-1,
        2.4411579970e-1,    2.5707310440e-1,    2.7003040910e-1,    2.8298771380e-1,
        2.9518282410e-1,    3.0814012890e-1,    3.1957304480e-1,    3.3100596070e-1,
        3.4015229340e-1,    3.4701204300e-1,    3.4999999400e-1,    3.8999998570e-1,
        3.9304876330e-1,    4.0143290160e-1,    4.1210362320e-1,    4.2429873350e-1,
        4.3649384380e-1,    4.4868895410e-1,    4.6164625880e-1,    4.7460356350e-1,
        4.8756086830e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,1,3,15
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n82_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        2.8614863750e-1,    3.4041404720e-3,    1.6632524130e-1,   -7.9737603660e-2,
        -1.7248192430e-1,    6.1409115790e-2,    9.4227790830e-3,    9.3253940340e-2,
        1.6842064260e-1,   -2.3573789000e-1,   -6.0760259630e-2,   -1.4159679410e-2,
        -3.8372337820e-2,    3.0130684380e-1,    1.9991576670e-2,   -2.0686031880e-1,
        2.0496517420e-2,   -2.4375110860e-1,    1.4708507060e-1,    3.0117571350e-1,
        -1.0720869900e-1,    2.3728579280e-2,   -1.5701508520e-1,   -2.8567618130e-1,
        3.2943701740e-1,    6.0640096660e-2,    5.9703886510e-2,    1.0485351090e-1,
        -3.6815273760e-1,   -9.1093122960e-2,    2.1457350250e-1,   -8.2264721390e-2,
        3.4748905900e-1,   -5.2387356760e-2,   -3.0297023060e-1,    5.2298903470e-3,
        -5.5031657220e-2,   -8.2062959670e-2,    6.0251021390e-1,   -1.7305946350e-1,
        3.0631327630e-1,
    ];
    let expected_deviations = vec![
        6.4827126260e-1,    1.6206781570e-1,    6.4827126260e-1,
    ];
    let expected_deviation_dbs = vec![
        4.3405742650e0,    -1.5806064610e1,     4.3405742650e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.2195123360e-2,    2.4390237410e-2,    3.7398356940e-2,
        4.9593467270e-2,    6.1788577590e-2,    7.3983691630e-2,    8.6178801950e-2,
        9.8373912270e-2,    1.1056902260e-1,    1.2276413290e-1,    1.3414624330e-1,
        1.4471533890e-1,    1.5000000600e-1,    1.8999999760e-1,    1.9487804170e-1,
        2.0382112260e-1,    2.1439021830e-1,    2.2495931390e-1,    2.3634141680e-1,
        2.4853652720e-1,    2.5991863010e-1,    2.7130073310e-1,    2.8268283610e-1,
        2.9487794640e-1,    3.0626004930e-1,    3.1764215230e-1,    3.2821124790e-1,
        3.3796733620e-1,    3.4609740970e-1,    3.4999999400e-1,    3.8999998570e-1,
        3.9731705190e-1,    4.0788614750e-1,    4.2008125780e-1,    4.3227636810e-1,
        4.4447147850e-1,    4.5666658880e-1,    4.6886169910e-1,    4.8024380210e-1,
        4.9243891240e-1,    4.9894297120e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,1,3,16
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n82_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        2.9183444380e-1,    3.9253830910e-3,    1.6965040560e-1,   -8.0916106700e-2,
        -1.7638015750e-1,    6.2492221590e-2,    9.1042816640e-3,    9.5635905860e-2,
        1.7193751040e-1,   -2.3971226810e-1,   -6.2746524810e-2,   -1.4623343940e-2,
        -3.9617300030e-2,    3.0787062640e-1,    2.1158128980e-2,   -2.1049138900e-1,
        1.9936352970e-2,   -2.4923846130e-1,    1.4927852150e-1,    3.0853003260e-1,
        -1.0857540370e-1,    2.4752855300e-2,   -1.6156762840e-1,   -2.9242557290e-1,
        3.3568638560e-1,    6.3250839710e-2,    6.2077879910e-2,    1.0743707420e-1,
        -3.7786984440e-1,   -9.3428194520e-2,    2.1823376420e-1,   -8.1578969960e-2,
        3.5565692190e-1,   -5.3908586500e-2,   -3.1160819530e-1,    5.5159926410e-3,
        -5.7459890840e-2,   -7.7630102630e-2,    6.1153489350e-1,   -1.7792177200e-1,
        3.0324113370e-1,
    ];
    let expected_deviations = vec![
        6.6208291050e-1,    1.6552072760e-1,    6.6208291050e-1,
    ];
    let expected_deviation_dbs = vec![
        4.4130539890e0,    -1.5622951510e1,     4.4130539890e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    1.2195123360e-2,    2.4390237410e-2,    3.7347543980e-2,
        4.9542654310e-2,    6.1737764630e-2,    7.3932878670e-2,    8.6127988990e-2,
        9.8323099320e-2,    1.1051820960e-1,    1.2271332000e-1,    1.3414624330e-1,
        1.4481696490e-1,    1.5000000600e-1,    1.8999999760e-1,    1.9457316400e-1,
        2.0371949670e-1,    2.1439021830e-1,    2.2506093980e-1,    2.3649385570e-1,
        2.4792677160e-1,    2.5935968760e-1,    2.7155479790e-1,    2.8298771380e-1,
        2.9442062970e-1,    3.0585354570e-1,    3.1728646160e-1,    3.2795718310e-1,
        3.3786571030e-1,    3.4624984860e-1,    3.4999999400e-1,    3.8999998570e-1,
        3.9685973530e-1,    4.0829265120e-1,    4.1972556710e-1,    4.3192067740e-1,
        4.4411578770e-1,    4.5631089810e-1,    4.6850600840e-1,    4.8070111870e-1,
        4.9213403460e-1,    4.9899378420e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,1,3,15
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n128_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 15);

    let expected_impulse_response = vec![
        -4.3577432630e-1,    8.3935856820e-3,   -1.4122599360e-1,    1.8629151580e-1,
        5.5289286380e-1,   -3.1599557400e-1,   -2.2429190580e-1,   -2.7339190240e-1,
        -4.7910499570e-1,     1.0566604140e0,    4.3248325590e-1,   -2.7777141330e-1,
        -4.7063469890e-2,    -1.5694651600e0,    1.3704633710e-1,     1.5592989920e0,
        1.1759898070e-1,    9.6637237070e-1,    -1.4197785850e0,    -2.3461766240e0,
        1.3035814760e0,    6.0027474160e-1,     1.7772367000e0,     1.3865345720e0,
        -3.7849349980e0,    -1.1348199840e0,    2.4299275880e-1,    6.7413485050e-1,
        4.7736654280e0,   -9.2101359370e-1,    -3.6573920250e0,   -8.5961949830e-1,
        -2.4812006950e0,     4.1672754290e0,     4.9118595120e0,    -2.6157364850e0,
        -1.3462765220e0,    -4.5104408260e0,    -2.0783288480e0,     7.4214935300e0,
        2.1398260590e0,    5.3632020950e-2,    -2.1212625500e0,    -8.3688344960e0,
        1.9346244340e0,     5.8757867810e0,     2.1404600140e0,     3.6434092520e0,
        -7.0511264800e0,    -7.1731557850e0,     3.5596575740e0,     2.4506645200e0,
        6.7935819630e0,     2.3701906200e0,    -1.0089677810e1,    -3.1534535880e0,
        -7.3758125310e-2,     3.3471791740e0,     1.0343073840e1,    -2.3490388390e0,
        -7.4652285580e0,    -2.6680126190e0,    -3.8631856440e0,     8.5359287260e0,
    ];
    let expected_deviations = vec![
        6.6486495730e-1,    1.6621623930e-1,    6.6486495730e-1,
    ];
    let expected_deviation_dbs = vec![
        4.4275803570e0,    -1.5586530690e1,     4.4275803570e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    7.8125000000e-3,    1.5625005590e-2,    2.3437498140e-2,
        3.1249990690e-2,    3.9062485100e-2,    4.6874977650e-2,    5.4687470200e-2,
        6.2499962750e-2,    7.0312455300e-2,    7.8124947850e-2,    8.5416607560e-2,
        9.3229100110e-2,    1.0104159270e-1,    1.0885408520e-1,    1.1614574490e-1,
        1.2343740460e-1,    1.3072913890e-1,    1.3750006260e-1,    1.4375014600e-1,
        1.4843770860e-1,    1.5000000600e-1,    1.8999999760e-1,    1.9208335880e-1,
        1.9729176160e-1,    2.0302100480e-1,    2.0979192850e-1,    2.1656285230e-1,
        2.2385461630e-1,    2.3062554000e-1,    2.3791730400e-1,    2.4572990830e-1,
        2.5302159790e-1,    2.6031315330e-1,    2.6760470870e-1,    2.7541708950e-1,
        2.8270864490e-1,    2.9000020030e-1,    2.9729175570e-1,    3.0458331110e-1,
        3.1187486650e-1,    3.1916642190e-1,    3.2593715190e-1,    3.3270788190e-1,
        3.3895778660e-1,    3.4468686580e-1,    3.4833264350e-1,    3.4999999400e-1,
        3.8999998570e-1,    3.9260411260e-1,    3.9781236650e-1,    4.0458309650e-1,
        4.1187465190e-1,    4.1916620730e-1,    4.2645776270e-1,    4.3427014350e-1,
        4.4156169890e-1,    4.4937407970e-1,    4.5718646050e-1,    4.6499884130e-1,
        4.7229039670e-1,    4.8010277750e-1,    4.8791515830e-1,    4.9520671370e-1,
        4.9937331680e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,1,3,16
0.000000000000,0.150000005960,0.189999997616,0.349999994040
0.389999985695,0.500000000000
1.000000000000,0.000000000000,1.000000000000
1.000000000000,4.000000000000,1.000000000000
*/
#[test]
fn multipleband_notch_n128_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.15f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.19f32,
        upper_edge: 0.35f32,
        desired_value: 0f32,
        weight: 4f32,
    });
    bands.push(Band {
        lower_edge: 0.39f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::MultipleBand, &bands, 16);

    let expected_impulse_response = vec![
        -4.3641445040e-1,    8.4998905660e-3,   -1.4148098230e-1,    1.8627345560e-1,
        5.5398875470e-1,   -3.1685459610e-1,   -2.2436386350e-1,   -2.7319908140e-1,
        -4.8048061130e-1,     1.0585031510e0,    4.3287295100e-1,   -2.7948957680e-1,
        -4.5154392720e-2,    -1.5718050000e0,    1.3703203200e-1,     1.5628697870e0,
        1.1414092780e-1,    9.6821492910e-1,    -1.4194850920e0,    -2.3511164190e0,
        1.3097027540e0,    5.9888434410e-1,     1.7744991780e0,     1.3926482200e0,
        -3.7931482790e0,    -1.1325008870e0,    2.4946856500e-1,    6.6537833210e-1,
        4.7816066740e0,   -9.2391681670e-1,    -3.6664514540e0,   -8.4591573480e-1,
        -2.4883761410e0,     4.1672134400e0,     4.9224100110e0,    -2.6341562270e0,
        -1.3371202950e0,    -4.5050449370e0,    -2.0916450020e0,     7.4404096600e0,
        2.1284039020e0,    4.3785333630e-2,    -2.1010956760e0,    -8.3865318300e0,
        1.9445304870e0,     5.8867359160e0,     2.1132063870e0,     3.6623337270e0,
        -7.0544672010e0,    -7.1855926510e0,     3.5890252590e0,     2.4278888700e0,
        6.7903947830e0,     2.3890130520e0,    -1.0116684910e1,    -3.1299853320e0,
        -6.8798303600e-2,     3.3194906710e0,     1.0369291310e1,    -2.3663656710e0,
        -7.4700517650e0,    -2.6360089780e0,    -3.8928277490e0,     8.5448369980e0,
    ];
    let expected_deviations = vec![
        6.6595274210e-1,    1.6648818550e-1,    6.6595274210e-1,
    ];
    let expected_deviation_dbs = vec![
        4.4332532880e0,    -1.5572332380e1,     4.4332532880e0,
    ];
    let expected_extremal_frequencies = vec![
        0.0000000000e0,    7.8125000000e-3,    1.5625000000e-2,    2.3437500000e-2,
        3.1250000000e-2,    3.9062500000e-2,    4.6875000000e-2,    5.4687500000e-2,
        6.2500000000e-2,    7.0312500000e-2,    7.8125000000e-2,    8.5449218750e-2,
        9.3261718750e-2,    1.0107421880e-1,    1.0839843750e-1,    1.1621093750e-1,
        1.2353515620e-1,    1.3085937500e-1,    1.3769531250e-1,    1.4355468750e-1,
        1.4794921880e-1,    1.5000000600e-1,    1.8999999760e-1,    1.9244140390e-1,
        1.9732421640e-1,    2.0318359140e-1,    2.0953124760e-1,    2.1636718510e-1,
        2.2369140390e-1,    2.3101562260e-1,    2.3833984140e-1,    2.4566406010e-1,
        2.5298827890e-1,    2.6031249760e-1,    2.6763671640e-1,    2.7544921640e-1,
        2.8277343510e-1,    2.9009765390e-1,    2.9742187260e-1,    3.0474609140e-1,
        3.1207031010e-1,    3.1939452890e-1,    3.2623046640e-1,    3.3306640390e-1,
        3.3892577890e-1,    3.4429687260e-1,    3.4869140390e-1,    3.4999999400e-1,
        3.8999998570e-1,    3.9244139190e-1,    3.9781248570e-1,    4.0464842320e-1,
        4.1148436070e-1,    4.1880857940e-1,    4.2662107940e-1,    4.3394529820e-1,
        4.4175779820e-1,    4.4957029820e-1,    4.5689451690e-1,    4.6470701690e-1,
        4.7251951690e-1,    4.8033201690e-1,    4.8765623570e-1,    4.9498045440e-1,
        4.9937498570e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_deviation_dbs(&pm_output, &expected_deviation_dbs);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n3_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        1.2830007080e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        3.3333337310e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.6666667160e-1,    3.3333334330e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n3_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        1.4963617920e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.6666667460e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.3333335070e-2,    4.6666666870e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n3_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        1.5017083290e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.7500005960e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.1250000000e-2,    4.6875000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n3_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        1.5053582190e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.8928890230e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3888888990e-2,    4.7222238780e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n3_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        1.5503965320e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        9.4594556090e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3513513840e-2,    4.8648640510e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n4_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -2.3522971200e-2,    2.0811291040e-1,
    ];
    let expected_deviations = vec![
        7.3456518350e-2,
    ];
    let expected_extremal_frequencies = vec![
        1.2500000000e-1,    3.7500000000e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n4_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -2.1930694580e-2,    2.0794318620e-1,
    ];
    let expected_deviations = vec![
        8.0504491930e-2,
    ];
    let expected_extremal_frequencies = vec![
        8.3333335820e-2,    3.3333334330e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n4_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -2.0865017550e-2,    2.0836640890e-1,
    ];
    let expected_deviations = vec![
        8.3074323830e-2,
    ];
    let expected_extremal_frequencies = vec![
        1.6666667540e-2,    3.1666672230e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n4_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -2.0859435200e-2,    2.0836457610e-1,
    ];
    let expected_deviations = vec![
        8.3103932440e-2,
    ];
    let expected_extremal_frequencies = vec![
        1.5625000000e-2,    3.2812500000e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n4_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        -2.0828133450e-2,    2.0836557450e-1,
    ];
    let expected_deviations = vec![
        8.3225160840e-2,
    ];
    let expected_extremal_frequencies = vec![
        6.9444444960e-3,    3.1944456700e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n4_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        -2.0827421920e-2,    2.0836259420e-1,
    ];
    let expected_deviations = vec![
        8.3239935340e-2,
    ];
    let expected_extremal_frequencies = vec![
        6.7567569200e-3,    3.2432425020e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n11_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        6.5518729390e-3,   -1.7719194290e-2,    3.4275103360e-2,   -6.5967135130e-2,
        1.5196879210e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.4619594220e-2,
    ];
    let expected_extremal_frequencies = vec![
        5.0000000750e-2,    1.0000000150e-1,    2.0000000300e-1,    3.0000001190e-1,
        3.5000002380e-1,    4.0000003580e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n11_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        3.5538535570e-2,   -5.7584494350e-2,    5.4814260450e-2,   -7.8577786680e-2,
        1.5824839470e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.6863399150e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.3333335070e-2,    1.0000000890e-1,    2.0000000300e-1,    3.0000001190e-1,
        4.0000000600e-1,    4.6666666870e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n11_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        7.9357676210e-2,   -1.1315567050e-1,    6.9273643200e-2,   -8.6063466970e-2,
        1.6174627840e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.6474267240e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.6666668280e-3,    1.0666667670e-1,    2.1333326400e-1,    3.1999984380e-1,
        4.2666640880e-1,    4.9333301190e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n11_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        8.0563113090e-2,   -1.1462855340e-1,    6.9540604950e-2,   -8.6071655150e-2,
        1.6171675920e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.7750891450e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.2500000930e-3,    1.0625001790e-1,    2.1249993150e-1,    3.1874984500e-1,
        4.2499974370e-1,    4.9374967810e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n11_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        8.2461714740e-2,   -1.1710360650e-1,    7.0155858990e-2,   -8.6404532190e-2,
        1.6184318070e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.9930156470e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.7777778450e-3,    1.0833333430e-1,    2.1388912200e-1,    3.2222273950e-1,
        4.2500078680e-1,    4.9444541340e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n11_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        8.2901775840e-2,   -1.1765115710e-1,    7.0320308210e-2,   -8.6509831250e-2,
        1.6191738840e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.0404168370e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.7027027680e-3,    1.0810814050e-1,    2.1351341900e-1,    3.2162174580e-1,
        4.2432484030e-1,    4.9459537860e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n12_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -2.6719036980e-3,    4.0246699940e-3,   -4.5096068640e-3,    8.4876883770e-3,
        -2.2726668040e-2,    2.0286192000e-1,
    ];
    let expected_deviations = vec![
        1.8870098520e-2,
    ];
    let expected_extremal_frequencies = vec![
        4.1666667910e-2,    1.2500000000e-1,    2.0833334330e-1,    2.9166665670e-1,
        3.7499997020e-1,    4.5833328370e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n12_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -2.3337167220e-3,    4.1035078470e-3,   -4.8021022230e-3,    8.6178928610e-3,
        -2.2864554080e-2,    2.0295949280e-1,
    ];
    let expected_deviations = vec![
        1.7274897550e-2,
    ];
    let expected_extremal_frequencies = vec![
        2.7777777980e-2,    1.1111111190e-1,    1.9444444780e-1,    2.7777779100e-1,
        3.6111116410e-1,    4.4444453720e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n12_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -2.3322761990e-3,    4.1650142520e-3,   -4.6351738270e-3,    8.4210420030e-3,
        -2.2772258150e-2,    2.0287729800e-1,
    ];
    let expected_deviations = vec![
        1.9187746570e-2,
    ];
    let expected_extremal_frequencies = vec![
        5.5555556900e-3,    9.4444446270e-2,    1.8888889250e-1,    2.8333342080e-1,
        3.7777811290e-1,    4.6111166480e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n12_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -2.3362378120e-3,    4.1693425740e-3,   -4.6319663520e-3,    8.4167057650e-3,
        -2.2768912840e-2,    2.0287269350e-1,
    ];
    let expected_deviations = vec![
        1.9216619430e-2,
    ];
    let expected_extremal_frequencies = vec![
        5.2083334890e-3,    9.3750007450e-2,    1.9270828370e-1,    2.8645828370e-1,
        3.7500011920e-1,    4.5833361150e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n12_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        -2.3314754940e-3,    4.1676256810e-3,   -4.6337908130e-3,    8.4185358140e-3,
        -2.2768648340e-2,    2.0287363230e-1,
    ];
    let expected_deviations = vec![
        1.9225206230e-2,
    ];
    let expected_extremal_frequencies = vec![
        2.3148148320e-3,    9.4907373190e-2,    1.8981492520e-1,    2.8472235800e-1,
        3.7731459740e-1,    4.5833280680e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n12_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        -2.3318585010e-3,    4.1678091510e-3,   -4.6329833570e-3,    8.4179621190e-3,
        -2.2768743340e-2,    2.0287379620e-1,
    ];
    let expected_deviations = vec![
        1.9227378070e-2,
    ];
    let expected_extremal_frequencies = vec![
        2.2522523070e-3,    9.4594560560e-2,    1.8918910620e-1,    2.8378364440e-1,
        3.7612593170e-1,    4.5945921540e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n16_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -1.7198007550e-3,    2.3658780850e-3,   -1.9610309970e-3,    2.7671703140e-3,
        -4.2860596440e-3,    8.2126623020e-3,   -2.2627817470e-2,    2.0272888240e-1,
    ];
    let expected_deviations = vec![
        1.3322710990e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.1250000000e-2,    6.2500000000e-2,    1.5625000000e-1,    2.1875000000e-1,
        2.8125000000e-1,    3.4375000000e-1,    4.0625000000e-1,    4.6875000000e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n16_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -1.4829445860e-3,    2.3435410110e-3,   -2.0808158440e-3,    2.7776712090e-3,
        -4.3269977900e-3,    8.3233900370e-3,   -2.2763034330e-2,    2.0291307570e-1,
    ];
    let expected_deviations = vec![
        1.1954023500e-2,
    ];
    let expected_extremal_frequencies = vec![
        2.0833333950e-2,    6.2500000000e-2,    1.4583332840e-1,    2.0833331350e-1,
        2.7083331350e-1,    3.5416668650e-1,    4.1666671630e-1,    4.5833340290e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n16_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -1.5309153820e-3,    2.4419734250e-3,   -2.0137743560e-3,    2.6988824830e-3,
        -4.2859390380e-3,    8.2362722610e-3,   -2.2634979340e-2,    2.0275488500e-1,
    ];
    let expected_deviations = vec![
        1.3609599320e-2,
    ];
    let expected_extremal_frequencies = vec![
        4.1666668840e-3,    7.0833340290e-2,    1.3750003280e-1,    2.0833329860e-1,
        2.7499991660e-1,    3.4583318230e-1,    4.1249978540e-1,    4.7083306310e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n16_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -1.5278203650e-3,    2.4397606030e-3,   -2.0158723930e-3,    2.7029989290e-3,
        -4.2884349820e-3,    8.2361511890e-3,   -2.2637169810e-2,    2.0275852080e-1,
    ];
    let expected_deviations = vec![
        1.3573216270e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    7.0312500000e-2,    1.3671875000e-1,    2.0703125000e-1,
        2.7734375000e-1,    3.4375000000e-1,    4.1015625000e-1,    4.7265625000e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n16_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        -1.5281887030e-3,    2.4415974040e-3,   -2.0141983400e-3,    2.6997313830e-3,
        -4.2867413720e-3,    8.2354824990e-3,   -2.2634726020e-2,    2.0275656880e-1,
    ];
    let expected_deviations = vec![
        1.3611136940e-2,
    ];
    let expected_extremal_frequencies = vec![
        1.7361111240e-3,    6.9444470110e-2,    1.3888888060e-1,    2.0833306010e-1,
        2.7604115010e-1,    3.4548532960e-1,    4.1145730020e-1,    4.7048485280e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n16_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        -1.5281020900e-3,    2.4416274390e-3,   -2.0145839080e-3,    2.7000124100e-3,
        -4.2865364810e-3,    8.2351723690e-3,   -2.2634565830e-2,    2.0275649430e-1,
    ];
    let expected_deviations = vec![
        1.3611694800e-2,
    ];
    let expected_extremal_frequencies = vec![
        1.6891892300e-3,    6.9256745280e-2,    1.3851352040e-1,    2.0777054130e-1,
        2.7702754740e-1,    3.4459537270e-1,    4.1216319800e-1,    4.7128504510e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n47_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        1.2864322400e-2,   -1.8106970940e-2,    1.0002977210e-2,   -7.7834427360e-3,
        7.0907678460e-3,   -7.0561431350e-3,    7.3854532090e-3,   -7.9465694730e-3,
        8.6706280710e-3,   -9.5232129100e-3,    1.0495759550e-2,   -1.1603213850e-2,
        1.2883830820e-2,   -1.4400359240e-2,    1.6244262460e-2,   -1.8546834590e-2,
        2.1504588430e-2,   -2.5435090070e-2,    3.0903466050e-2,   -3.9039947090e-2,
        5.2489623430e-2,   -7.9204142090e-2,    1.5896905960e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.0991412400e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0869565420e-2,    2.1739130840e-2,    4.3478261680e-2,    6.5217390660e-2,
        8.6956515910e-2,    1.0869564120e-1,    1.3043476640e-1,    1.5217389170e-1,
        1.7391301690e-1,    1.9565214220e-1,    2.1739126740e-1,    2.3913039270e-1,
        2.6086953280e-1,    2.8260865810e-1,    3.0434778330e-1,    3.2608690860e-1,
        3.4782603380e-1,    3.6956515910e-1,    3.9130428430e-1,    4.2391297220e-1,
        4.4565209750e-1,    4.5652166010e-1,    4.7826078530e-1,    4.8913034800e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n47_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        2.2639248520e-2,   -2.8685279190e-2,    1.0387834160e-2,   -8.4364339710e-3,
        9.3237832190e-3,   -9.3021243810e-3,    9.2745088040e-3,   -1.0039731860e-2,
        1.0874196890e-2,   -1.1432327330e-2,    1.2235086410e-2,   -1.3410419230e-2,
        1.4620944860e-2,   -1.5950575470e-2,    1.7748773100e-2,   -2.0033322270e-2,
        2.2815465930e-2,   -2.6551321150e-2,    3.1907692550e-2,   -3.9869651200e-2,
        5.3062081340e-2,   -7.9586997630e-2,    1.5920087700e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.2017829120e-1,
    ];
    let expected_extremal_frequencies = vec![
        7.2463769470e-3,    2.1739130840e-2,    4.3478257950e-2,    6.5217383210e-2,
        8.6956508460e-2,    1.0869563370e-1,    1.3043476640e-1,    1.5217389170e-1,
        1.7391301690e-1,    1.9565214220e-1,    2.1739126740e-1,    2.4637676780e-1,
        2.6811590790e-1,    2.8985503320e-1,    3.1159415840e-1,    3.3333328370e-1,
        3.5507240890e-1,    3.7681153420e-1,    3.9855065940e-1,    4.2028978470e-1,
        4.4202890990e-1,    4.6376803520e-1,    4.8550716040e-1,    4.9275353550e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n47_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        5.1063861700e-2,   -6.3241124150e-2,    1.9275132570e-2,   -1.2789063160e-2,
        1.1074475940e-2,   -1.0621905330e-2,    1.0663993660e-2,   -1.0974474250e-2,
        1.1455446480e-2,   -1.2077108030e-2,    1.2829713520e-2,   -1.3823956250e-2,
        1.4905519780e-2,   -1.6189411280e-2,    1.7962679270e-2,   -2.0179599520e-2,
        2.2997513410e-2,   -2.6746541260e-2,    3.2011106610e-2,   -3.9918169380e-2,
        5.3130462770e-2,   -7.9622313380e-2,    1.5920534730e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        5.3746765850e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.4492754130e-3,    2.1739134560e-2,    4.4927544890e-2,    6.6666670140e-2,
        8.8405750690e-2,    1.1014483120e-1,    1.3188391920e-1,    1.5507227180e-1,
        1.7681135240e-1,    1.9855043290e-1,    2.2028951350e-1,    2.4347786610e-1,
        2.6521709560e-1,    2.8695639970e-1,    3.1014499070e-1,    3.3188429470e-1,
        3.5362359880e-1,    3.7536290290e-1,    3.9855149390e-1,    4.2029079790e-1,
        4.4203010200e-1,    4.6376940610e-1,    4.8550871010e-1,    4.9710300560e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n47_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        6.9869861010e-2,   -8.5616104300e-2,    2.4293608960e-2,   -1.5036284920e-2,
        1.2423738840e-2,   -1.1576481160e-2,    1.1403888460e-2,   -1.1566057800e-2,
        1.1938355860e-2,   -1.2491866950e-2,    1.3201907280e-2,   -1.4115989210e-2,
        1.5170499680e-2,   -1.6583219170e-2,    1.8398195510e-2,   -2.0484030250e-2,
        2.3170217870e-2,   -2.6835784320e-2,    3.2060876490e-2,   -3.9979070430e-2,
        5.3219497200e-2,   -7.9683959480e-2,    1.5920576450e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.4056428670e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3586956770e-3,    2.1739127110e-2,    4.4836949560e-2,    6.6576078530e-2,
        8.8315203790e-2,    1.1005432900e-1,    1.3315218690e-1,    1.5489143130e-1,
        1.7663067580e-1,    1.9836992030e-1,    2.2146786750e-1,    2.4320711200e-1,
        2.6494619250e-1,    2.8804388640e-1,    3.0978289250e-1,    3.3152189850e-1,
        3.5326090460e-1,    3.7635859850e-1,    3.9809760450e-1,    4.1983661060e-1,
        4.4293430450e-1,    4.6467331050e-1,    4.8505362870e-1,    4.9864050750e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n47_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        8.1952489910e-2,   -9.9996097390e-2,    2.7437612410e-2,   -1.6418181360e-2,
        1.3249926270e-2,   -1.2148126960e-2,    1.1834383010e-2,   -1.1907622220e-2,
        1.2223869560e-2,   -1.2727990750e-2,    1.3410776850e-2,   -1.4283671980e-2,
        1.5364468100e-2,   -1.6699925070e-2,    1.8356144430e-2,   -2.0451381800e-2,
        2.3205921050e-2,   -2.6937484740e-2,    3.2158866520e-2,   -4.0040314200e-2,
        5.3242936730e-2,   -7.9696625470e-2,    1.5921375160e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.7255674600e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.0386472610e-4,    2.2342989220e-2,    4.4082131240e-2,    6.6425144670e-2,
        8.8164180520e-2,    1.1050707850e-1,    1.3285006580e-1,    1.5458936990e-1,
        1.7693254350e-1,    1.9867184760e-1,    2.2101502120e-1,    2.4335819480e-1,
        2.6509711150e-1,    2.8743973370e-1,    3.0978235600e-1,    3.3152112360e-1,
        3.5386374590e-1,    3.7620636820e-1,    3.9854899050e-1,    4.2028775810e-1,
        4.4263038040e-1,    4.6497300270e-1,    4.8610791560e-1,    4.9939271810e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n47_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        7.2585925460e-2,   -8.8881954550e-2,    2.5001287460e-2,   -1.5354596080e-2,
        1.2626536190e-2,   -1.1719778180e-2,    1.1508643630e-2,   -1.1648677290e-2,
        1.2011267240e-2,   -1.2546628710e-2,    1.3249188660e-2,   -1.4128938320e-2,
        1.5224441890e-2,   -1.6589537260e-2,    1.8278032540e-2,   -2.0395353440e-2,
        2.3150280120e-2,   -2.6886671780e-2,    3.2126173380e-2,   -4.0011525150e-2,
        5.3222104910e-2,   -7.9690903430e-2,    1.5921285750e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.7117580180e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.8754405470e-4,    2.2326670590e-2,    4.4065818190e-2,    6.6392540930e-2,
        8.8131718340e-2,    1.1045844110e-1,    1.3278506700e-1,    1.5452396870e-1,
        1.7685040830e-1,    1.9917684790e-1,    2.2091574970e-1,    2.4324218930e-1,
        2.6556903120e-1,    2.8730848430e-1,    3.0963549020e-1,    3.3196249600e-1,
        3.5370194910e-1,    3.7602895500e-1,    3.9835596080e-1,    4.2068296670e-1,
        4.4242241980e-1,    4.6474942560e-1,    4.8590132590e-1,    4.9882748720e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n48_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -4.0149121200e-4,    5.0405156800e-4,   -2.3123521530e-4,    2.2859827730e-4,
        -2.0114796640e-4,    1.8576800360e-4,   -2.0695140120e-4,    2.0760379270e-4,
        -2.2932252610e-4,    2.5960677890e-4,   -2.8588157150e-4,    3.3542269380e-4,
        -3.9118347920e-4,    4.6311994080e-4,   -5.6868477260e-4,    7.0402881830e-4,
        -9.0331386310e-4,    1.2036906560e-3,   -1.6752092630e-3,    2.5046230290e-3,
        -4.1379891340e-3,    8.1056896600e-3,   -2.2518979390e-2,    2.0264315610e-1,
    ];
    let expected_deviations = vec![
        3.6127059720e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0416666980e-2,    2.0833333950e-2,    4.1666667910e-2,    6.2500000000e-2,
        8.3333328370e-2,    1.0416665670e-1,    1.2499998510e-1,    1.5625000000e-1,
        1.7708334330e-1,    1.9791668650e-1,    2.1875002980e-1,    2.3958337310e-1,
        2.6041668650e-1,    2.8125000000e-1,    3.0208331350e-1,    3.2291662690e-1,
        3.4374994040e-1,    3.6458325390e-1,    3.8541656730e-1,    4.0624988080e-1,
        4.2708319430e-1,    4.4791650770e-1,    4.6874982120e-1,    4.8958313470e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n48_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -3.3905246530e-4,    4.3620925860e-4,   -1.8424802690e-4,    1.5653661100e-4,
        -1.5913933750e-4,    1.8551429090e-4,   -1.9980435900e-4,    2.0815979220e-4,
        -2.3313824200e-4,    2.6505882850e-4,   -3.0806136780e-4,    3.5273889080e-4,
        -4.0916152650e-4,    4.8707932000e-4,   -5.9056433380e-4,    7.3327834250e-4,
        -9.3134114290e-4,    1.2296598870e-3,   -1.7063212580e-3,    2.5341969450e-3,
        -4.1690682990e-3,    8.1380335610e-3,   -2.2548912090e-2,    2.0267559590e-1,
    ];
    let expected_deviations = vec![
        3.2764342610e-3,
    ];
    let expected_extremal_frequencies = vec![
        6.9444444960e-3,    2.0833333950e-2,    4.1666667910e-2,    6.2500000000e-2,
        8.3333343270e-2,    1.1111113430e-1,    1.3194447760e-1,    1.5277782080e-1,
        1.7361116410e-1,    1.9444450740e-1,    2.1527785060e-1,    2.3611119390e-1,
        2.5694453720e-1,    2.7777788040e-1,    3.0555567150e-1,    3.2638901470e-1,
        3.4722235800e-1,    3.6805570130e-1,    3.8888904450e-1,    4.0972238780e-1,
        4.3055573110e-1,    4.5138907430e-1,    4.7222241760e-1,    4.8611131310e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n48_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -3.9166532220e-4,    5.0887651740e-4,   -2.0222045710e-4,    1.6453681750e-4,
        -1.6201339890e-4,    1.7003476390e-4,   -1.8391838240e-4,    2.0259080340e-4,
        -2.2591362360e-4,    2.5488849500e-4,   -2.9107919540e-4,    3.3667261600e-4,
        -3.9494797240e-4,    4.7092366730e-4,   -5.7239551100e-4,    7.1196863430e-4,
        -9.1110949870e-4,    1.2092568210e-3,   -1.6847185320e-3,    2.5114086460e-3,
        -4.1446206160e-3,    8.1152459610e-3,   -2.2526115180e-2,    2.0265260340e-1,
    ];
    let expected_deviations = vec![
        4.0011010130e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.3888889230e-3,    2.2222222760e-2,    4.3055556710e-2,    6.5277785060e-2,
        8.6111173030e-2,    1.0833345350e-1,    1.2916684150e-1,    1.5138912200e-1,
        1.7222251000e-1,    1.9444479050e-1,    2.1666707100e-1,    2.3750045900e-1,
        2.5972262020e-1,    2.8055578470e-1,    3.0277782680e-1,    3.2361099120e-1,
        3.4583303330e-1,    3.6805507540e-1,    3.8888823990e-1,    4.1111028190e-1,
        4.3194344640e-1,    4.5277661090e-1,    4.7360977530e-1,    4.9166518450e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n48_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -3.9077783000e-4,    5.0790153910e-4,   -2.0211632360e-4,    1.6467674870e-4,
        -1.6228968160e-4,    1.7049242160e-4,   -1.8435547830e-4,    2.0282683540e-4,
        -2.2620939130e-4,    2.5511125570e-4,   -2.9113510390e-4,    3.3685343810e-4,
        -3.9525720060e-4,    4.7123545660e-4,   -5.7269679380e-4,    7.1222975380e-4,
        -9.1137166600e-4,    1.2093492550e-3,   -1.6845050270e-3,    2.5119397320e-3,
        -4.1459091010e-3,    8.1154648210e-3,   -2.2525392470e-2,    2.0265181360e-1,
    ];
    let expected_deviations = vec![
        3.9922827850e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.3020833720e-3,    2.2135417910e-2,    4.2968742550e-2,    6.5104141830e-2,
        8.5937514900e-2,    1.0807297380e-1,    1.3020840290e-1,    1.5104165670e-1,
        1.7317698900e-1,    1.9401024280e-1,    2.1614557500e-1,    2.3828090730e-1,
        2.5911426540e-1,    2.8124985100e-1,    3.0208334330e-1,    3.2421892880e-1,
        3.4635451440e-1,    3.6718800660e-1,    3.8932359220e-1,    4.1015708450e-1,
        4.3229267000e-1,    4.5312616230e-1,    4.7395965460e-1,    4.9088686700e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n48_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        -3.9124983600e-4,    5.0866807580e-4,   -2.0224752370e-4,    1.6466504890e-4,
        -1.6218451490e-4,    1.7033341280e-4,   -1.8415585510e-4,    2.0259107990e-4,
        -2.2584814000e-4,    2.5481486230e-4,   -2.9097544030e-4,    3.3656263260e-4,
        -3.9488592300e-4,    4.7093071040e-4,   -5.7239597660e-4,    7.1194628250e-4,
        -9.1104902090e-4,    1.2091645040e-3,   -1.6847208610e-3,    2.5116619650e-3,
        -4.1453391310e-3,    8.1153661010e-3,   -2.2525411100e-2,    2.0265191790e-1,
    ];
    let expected_deviations = vec![
        4.0035708810e-3,
    ];
    let expected_extremal_frequencies = vec![
        5.7870370800e-4,    2.1412029860e-2,    4.3402794750e-2,    6.4814873040e-2,
        8.6226828400e-2,    1.0821748520e-1,    1.2962944810e-1,    1.5104140340e-1,
        1.7303206030e-1,    1.9444401560e-1,    2.1585597100e-1,    2.3784662780e-1,
        2.5925859810e-1,    2.8067055340e-1,    3.0266121030e-1,    3.2407316570e-1,
        3.4548512100e-1,    3.6747577790e-1,    3.8888773320e-1,    4.1029968860e-1,
        4.3171164390e-1,    4.5312359930e-1,    4.7337815170e-1,    4.9131789800e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n48_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        -3.9121959710e-4,    5.0862046190e-4,   -2.0221396700e-4,    1.6465166120e-4,
        -1.6218711970e-4,    1.7034442860e-4,   -1.8416909730e-4,    2.0260916790e-4,
        -2.2589587020e-4,    2.5488273240e-4,   -2.9101147080e-4,    3.3654266740e-4,
        -3.9483892030e-4,    4.7083268870e-4,   -5.7228805960e-4,    7.1185606070e-4,
        -9.1104593590e-4,    1.2093440160e-3,   -1.6848018860e-3,    2.5115618480e-3,
        -4.1452921000e-3,    8.1154154610e-3,   -2.2525504230e-2,    2.0265203710e-1,
    ];
    let expected_deviations = vec![
        4.0033212860e-3,
    ];
    let expected_extremal_frequencies = vec![
        5.6306307670e-4,    2.1396389230e-2,    4.3355837460e-2,    6.4752221110e-2,
        8.6148604750e-2,    1.0810805110e-1,    1.2950448690e-1,    1.5146422390e-1,
        1.7286089060e-1,    1.9425755740e-1,    2.1621729430e-1,    2.3761396110e-1,
        2.5901037450e-1,    2.8096953030e-1,    3.0236563090e-1,    3.2432478670e-1,
        3.4572088720e-1,    3.6711698770e-1,    3.8907614350e-1,    4.1047224400e-1,
        4.3186834450e-1,    4.5326444510e-1,    4.7353443500e-1,    4.9155220390e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,2,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n81_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -1.1761515400e-2,    1.5921104700e-2,   -6.9101760160e-3,    4.1372757410e-3,
        -3.3569689840e-3,    3.4866761420e-3,   -3.9229877290e-3,    4.2912065980e-3,
        -4.4412482530e-3,    4.4156201180e-3,   -4.3569915000e-3,    4.3974705040e-3,
        -4.5852325860e-3,    4.8779249190e-3,   -5.1911771300e-3,    5.4617561400e-3,
        -5.6842491030e-3,    5.9039630000e-3,   -6.1789713800e-3,    6.5413974230e-3,
        -6.9832652810e-3,    7.4721388520e-3,   -7.9812705520e-3,    8.5118636490e-3,
        -9.0939328070e-3,    9.7687691450e-3,   -1.0568700730e-2,    1.1509418490e-2,
        -1.2599751350e-2,    1.3861551880e-2,   -1.5347346660e-2,    1.7149001360e-2,
        -1.9402287900e-2,    2.2303506730e-2,   -2.6162922380e-2,    3.1537055970e-2,
        -3.9554581050e-2,    5.2869573240e-2,   -7.9448297620e-2,    1.5908662970e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.0606078800e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.2500000930e-3,    1.2500000190e-2,    2.5000000370e-2,    3.7500001490e-2,
        5.0000004470e-2,    6.2500007450e-2,    7.5000010430e-2,    8.7500013410e-2,
        1.0000001640e-1,    1.1250001940e-1,    1.2500001490e-1,    1.3750000300e-1,
        1.4999999110e-1,    1.6249997910e-1,    1.7499996720e-1,    1.8749995530e-1,
        1.9999994340e-1,    2.1249993150e-1,    2.2499991950e-1,    2.3749990760e-1,
        2.4999989570e-1,    2.6249989870e-1,    2.7499988680e-1,    2.8749987480e-1,
        2.9999986290e-1,    3.1249985100e-1,    3.2499983910e-1,    3.3749982710e-1,
        3.4999981520e-1,    3.6874979730e-1,    3.8124978540e-1,    3.9374977350e-1,
        4.0624976160e-1,    4.1874974970e-1,    4.3124973770e-1,    4.4374972580e-1,
        4.5624971390e-1,    4.6874970200e-1,    4.8124969010e-1,    4.8749968410e-1,
        4.9374967810e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,2,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n81_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -2.1384071560e-2,    2.6087190960e-2,   -7.2077363730e-3,    5.0335861740e-3,
        -5.4889768360e-3,    4.9112178390e-3,   -4.5052915810e-3,    4.8723332580e-3,
        -5.0937756900e-3,    5.0053708260e-3,   -5.1781423390e-3,    5.5204220120e-3,
        -5.6605860590e-3,    5.7843625550e-3,   -6.1001405120e-3,    6.3960030670e-3,
        -6.5897032620e-3,    6.8811997770e-3,   -7.2680115700e-3,    7.5965002180e-3,
        -7.9347789290e-3,    8.3938464520e-3,   -8.8867619630e-3,    9.3696042900e-3,
        -9.9524706600e-3,    1.0653071110e-2,   -1.1400640010e-2,    1.2248896060e-2,
        -1.3289935890e-2,    1.4509268110e-2,   -1.5930369500e-2,    1.7695724960e-2,
        -1.9927397370e-2,    2.2758439180e-2,   -2.6529058810e-2,    3.1848013400e-2,
        -3.9810702200e-2,    5.3051963450e-2,   -7.9577341680e-2,    1.5916895870e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.1571552750e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.1666668840e-3,    1.2500001120e-2,    2.5000000370e-2,    3.7500001490e-2,
        5.0000000750e-2,    6.2500000000e-2,    7.5000010430e-2,    8.7500020860e-2,
        1.0000003130e-1,    1.1250004170e-1,    1.2500004470e-1,    1.3750003280e-1,
        1.5000002090e-1,    1.6250000890e-1,    1.7499999700e-1,    1.8749998510e-1,
        1.9999997320e-1,    2.1249996130e-1,    2.2916661200e-1,    2.4166660010e-1,
        2.5416660310e-1,    2.6666659120e-1,    2.7916657920e-1,    2.9166656730e-1,
        3.0416655540e-1,    3.1666654350e-1,    3.2916653160e-1,    3.4166651960e-1,
        3.5416650770e-1,    3.6666649580e-1,    3.7916648390e-1,    3.9166647200e-1,
        4.0416646000e-1,    4.1666644810e-1,    4.2916643620e-1,    4.4166642430e-1,
        4.5416641240e-1,    4.6666640040e-1,    4.7916638850e-1,    4.9166637660e-1,
        4.9583303930e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,2,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n81_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -6.7338176070e-2,    8.1169933080e-2,   -2.0381845530e-2,    1.0975174610e-2,
        -8.1223845480e-3,    6.9717541340e-3,   -6.4458772540e-3,    6.2045156960e-3,
        -6.1072483660e-3,    6.0941576960e-3,   -6.1252266170e-3,    6.2606781720e-3,
        -6.3136890530e-3,    6.4225420360e-3,   -6.7452862860e-3,    7.0042535660e-3,
        -7.2010606530e-3,    7.4086934330e-3,   -7.6586455110e-3,    7.9531371590e-3,
        -8.2983225580e-3,    8.6943954230e-3,   -9.1495066880e-3,    9.6372365950e-3,
        -1.0206893090e-2,    1.0868340730e-2,   -1.1581912640e-2,    1.2424245480e-2,
        -1.3446792960e-2,    1.4658331870e-2,   -1.6094401480e-2,    1.7843991520e-2,
        -2.0030140880e-2,    2.2841960190e-2,   -2.6601508260e-2,    3.1894668940e-2,
        -3.9843559270e-2,    5.3090542550e-2,   -7.9611599450e-2,    1.5918165450e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.2092562910e-1,
    ];
    let expected_extremal_frequencies = vec![
        8.3333335350e-4,    1.2500001120e-2,    2.4999992920e-2,    3.7499982860e-2,
        5.0833303480e-2,    6.3333295290e-2,    7.5833283360e-2,    8.8333271440e-2,
        1.0083325950e-1,    1.1333324760e-1,    1.2583324310e-1,    1.3916656370e-1,
        1.5166655180e-1,    1.6416653990e-1,    1.7666652800e-1,    1.8916651610e-1,
        2.0166650410e-1,    2.1416649220e-1,    2.2749981280e-1,    2.3999980090e-1,
        2.5249978900e-1,    2.6499977710e-1,    2.7749976520e-1,    2.8999975320e-1,
        3.0333307390e-1,    3.1583306190e-1,    3.2833305000e-1,    3.4083303810e-1,
        3.5333302620e-1,    3.6583301420e-1,    3.7916633490e-1,    3.9166632290e-1,
        4.0416631100e-1,    4.1666629910e-1,    4.2916628720e-1,    4.4249960780e-1,
        4.5499959590e-1,    4.6749958400e-1,    4.7999957200e-1,    4.9166622760e-1,
        4.9916622040e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,2,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n81_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -6.8560376760e-2,    8.2631506030e-2,   -2.0718932150e-2,    1.1142469940e-2,
        -8.2421079280e-3,    7.0742070670e-3,   -6.5380930900e-3,    6.2862634660e-3,
        -6.1868429180e-3,    6.1782300470e-3,   -6.2225088480e-3,    6.3260048630e-3,
        -6.3895359640e-3,    6.5874159340e-3,   -6.8375095730e-3,    6.9274604320e-3,
        -7.0752650500e-3,    7.3282271620e-3,   -7.6385587450e-3,    7.9790353780e-3,
        -8.3438754080e-3,    8.7415575980e-3,   -9.1888904570e-3,    9.6830427650e-3,
        -1.0253816840e-2,    1.0878041390e-2,   -1.1605650190e-2,    1.2490227820e-2,
        -1.3484045860e-2,    1.4644533400e-2,   -1.6072735190e-2,    1.7836138610e-2,
        -2.0038872960e-2,    2.2866278890e-2,   -2.6631399990e-2,    3.1919956210e-2,
        -3.9859622720e-2,    5.3110420700e-2,   -7.9613268380e-2,    1.5915825960e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.3425561190e-1,
    ];
    let expected_extremal_frequencies = vec![
        7.8125001160e-4,    1.2500002050e-2,    2.4999992920e-2,    3.7499982860e-2,
        5.0781220200e-2,    6.3281208280e-2,    7.5781255960e-2,    8.8281303640e-2,
        1.0078135130e-1,    1.1328139900e-1,    1.2656269970e-1,    1.3906274740e-1,
        1.5156279500e-1,    1.6406284270e-1,    1.7656289040e-1,    1.8906293810e-1,
        2.0156298580e-1,    2.1484428640e-1,    2.2734433410e-1,    2.3984438180e-1,
        2.5234436990e-1,    2.6484417920e-1,    2.7734398840e-1,    2.9062503580e-1,
        3.0312484500e-1,    3.1562465430e-1,    3.2812446360e-1,    3.4062427280e-1,
        3.5390532020e-1,    3.6640512940e-1,    3.7890493870e-1,    3.9140474800e-1,
        4.0390455720e-1,    4.1718560460e-1,    4.2968541380e-1,    4.4218522310e-1,
        4.5468503240e-1,    4.6718484160e-1,    4.7968465090e-1,    4.9218446020e-1,
        4.9921560290e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,2,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n82_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        2.1677627230e-4,   -2.7002239950e-4,    1.0620077960e-4,   -9.6150288300e-5,
        7.4083211080e-5,   -5.5207776310e-5,    5.9346173660e-5,   -5.1139246350e-5,
        5.1405626440e-5,   -5.7284683860e-5,    5.6197695810e-5,   -6.2945488030e-5,
        6.8847097280e-5,   -7.1824280890e-5,    8.0452970000e-5,   -8.5871026390e-5,
        9.1451336630e-5,   -1.0058869520e-4,    1.0670483610e-4,   -1.1548111800e-4,
        1.2652852460e-4,   -1.3645940630e-4,    1.5110048120e-4,   -1.6786009660e-4,
        1.8662527150e-4,   -2.1220321650e-4,    2.4165905780e-4,   -2.7813826450e-4,
        3.2572646160e-4,   -3.8394116560e-4,    4.6073237900e-4,   -5.6354491970e-4,
        7.0245750250e-4,   -9.0245268070e-4,    1.2009351050e-3,   -1.6753978560e-3,
        2.5029168460e-3,   -4.1359975000e-3,    8.1052454190e-3,   -2.2516006600e-2,
        2.0264166590e-1,
    ];
    let expected_deviations = vec![
        2.0176733380e-3,
    ];
    let expected_extremal_frequencies = vec![
        6.0975607480e-3,    1.2195121500e-2,    2.4390242990e-2,    3.6585364490e-2,
        4.8780485990e-2,    6.0975611210e-2,    7.3170736430e-2,    8.5365861650e-2,
        9.7560986880e-2,    1.0975611210e-1,    1.2195123730e-1,    1.4024390280e-1,
        1.5243901310e-1,    1.6463412340e-1,    1.7682923380e-1,    1.8902434410e-1,
        2.0121945440e-1,    2.1341456470e-1,    2.2560967500e-1,    2.3780478540e-1,
        2.4999989570e-1,    2.6219502090e-1,    2.7439013120e-1,    2.8658524160e-1,
        2.9878035190e-1,    3.1097546220e-1,    3.2317057250e-1,    3.3536568280e-1,
        3.4756079320e-1,    3.5975590350e-1,    3.7195101380e-1,    3.8414612410e-1,
        3.9634123440e-1,    4.0853634480e-1,    4.2073145510e-1,    4.3292656540e-1,
        4.5121923090e-1,    4.6341434120e-1,    4.7560945150e-1,    4.8780456190e-1,
        4.9390211700e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,2,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n82_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        2.1998726880e-4,   -2.7115351990e-4,    8.7336906290e-5,   -6.1103870390e-5,
        5.5238058850e-5,   -6.2000734030e-5,    5.8433193770e-5,   -5.0059505160e-5,
        5.2811272330e-5,   -5.6663164290e-5,    6.2105485990e-5,   -6.2833598350e-5,
        6.4228894190e-5,   -7.0122172470e-5,    7.6084463220e-5,   -8.2491300420e-5,
        8.7055304900e-5,   -9.3770839160e-5,    1.0327067870e-4,   -1.1323906070e-4,
        1.2421196150e-4,   -1.3567684800e-4,    1.5063931640e-4,   -1.6887221140e-4,
        1.8969154920e-4,   -2.1425468730e-4,    2.4362523980e-4,   -2.8108648260e-4,
        3.2785249640e-4,   -3.8669450440e-4,    4.6283338450e-4,   -5.6428788230e-4,
        7.0462247820e-4,   -9.0434553570e-4,    1.2027232440e-3,   -1.6780572480e-3,
        2.5049271060e-3,   -4.1391602720e-3,    8.1094373020e-3,   -2.2519480440e-2,
        2.0264565940e-1,
    ];
    let expected_deviations = vec![
        2.2075085440e-3,
    ];
    let expected_extremal_frequencies = vec![
        4.0650404990e-3,    1.2195121500e-2,    2.4390242990e-2,    3.6585364490e-2,
        4.8780485990e-2,    6.0975607480e-2,    7.3170721530e-2,    8.5365831850e-2,
        1.0162597890e-1,    1.1382108930e-1,    1.2601619960e-1,    1.3821130990e-1,
        1.5040642020e-1,    1.6260153060e-1,    1.7479664090e-1,    1.8699175120e-1,
        1.9918686150e-1,    2.1138197180e-1,    2.2357708220e-1,    2.3577219250e-1,
        2.4796730280e-1,    2.6016241310e-1,    2.7235752340e-1,    2.8861767050e-1,
        3.0081278090e-1,    3.1300789120e-1,    3.2520300150e-1,    3.3739811180e-1,
        3.4959322210e-1,    3.6178833250e-1,    3.7398344280e-1,    3.8617855310e-1,
        3.9837366340e-1,    4.1056877370e-1,    4.2276388410e-1,    4.3495899440e-1,
        4.4715410470e-1,    4.5934921500e-1,    4.7154432540e-1,    4.8373943570e-1,
        4.9593454600e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,2,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n82_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        2.1714779720e-4,   -2.7087601480e-4,    8.6070023830e-5,   -5.9119443900e-5,
        5.2234150640e-5,   -5.0575814380e-5,    5.0942297090e-5,   -5.2339826650e-5,
        5.4391486630e-5,   -5.6947174020e-5,    5.9954847530e-5,   -6.3407183920e-5,
        6.7321496320e-5,   -7.1740243580e-5,    7.6726093540e-5,   -8.2351820310e-5,
        8.8711138230e-5,   -9.5925468490e-5,    1.0414556890e-4,   -1.1355342580e-4,
        1.2438747220e-4,   -1.3695011150e-4,    1.5162487400e-4,   -1.6891624550e-4,
        1.8949601510e-4,   -2.1427412870e-4,    2.4402035340e-4,   -2.8108383410e-4,
        3.2780756010e-4,   -3.8652599320e-4,    4.6271248720e-4,   -5.6447123640e-4,
        7.0429279000e-4,   -9.0372050180e-4,    1.2021360450e-3,   -1.6777780840e-3,
        2.5047813540e-3,   -4.1385623630e-3,    8.1086829300e-3,   -2.2518793120e-2,
        2.0264533160e-1,
    ];
    let expected_deviations = vec![
        2.2766741460e-3,
    ];
    let expected_extremal_frequencies = vec![
        8.1300811140e-4,    1.2195123360e-2,    2.5203244760e-2,    3.7398356940e-2,
        4.9593467270e-2,    6.2601588670e-2,    7.4796698990e-2,    8.6991809310e-2,
        9.9999926980e-2,    1.1219503730e-1,    1.2439014760e-1,    1.3739827280e-1,
        1.4959338310e-1,    1.6178849340e-1,    1.7479661110e-1,    1.8699172140e-1,
        1.9918683170e-1,    2.1219494940e-1,    2.2439005970e-1,    2.3658517000e-1,
        2.4959328770e-1,    2.6178839800e-1,    2.7398350830e-1,    2.8699162600e-1,
        2.9918673630e-1,    3.1138184670e-1,    3.2438996430e-1,    3.3658507470e-1,
        3.4878018500e-1,    3.6178830270e-1,    3.7398341300e-1,    3.8617852330e-1,
        3.9918664100e-1,    4.1138175130e-1,    4.2357686160e-1,    4.3658497930e-1,
        4.4878008960e-1,    4.6097519990e-1,    4.7317031030e-1,    4.8455241320e-1,
        4.9512150880e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,2,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n82_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        2.1656911120e-4,   -2.7021823920e-4,    8.5973413660e-5,   -5.9114834580e-5,
        5.2211224100e-5,   -5.0555248890e-5,    5.0972186730e-5,   -5.2367278840e-5,
        5.4376141630e-5,   -5.6912351280e-5,    5.9908175900e-5,   -6.3350446000e-5,
        6.7270666480e-5,   -7.1686146840e-5,    7.6656360760e-5,   -8.2280355850e-5,
        8.8639440950e-5,   -9.5872688690e-5,    1.0414799910e-4,   -1.1360659850e-4,
        1.2447875630e-4,   -1.3707217290e-4,    1.5172973510e-4,   -1.6903242790e-4,
        1.8971515240e-4,   -2.1441976420e-4,    2.4436495730e-4,   -2.8140703220e-4,
        3.2734952400e-4,   -3.8616277740e-4,    4.6301132530e-4,   -5.6483957450e-4,
        7.0455262900e-4,   -9.0391666160e-4,    1.2022892480e-3,   -1.6779412980e-3,
        2.5049927640e-3,   -4.1387933310e-3,    8.1089213490e-3,   -2.2519042720e-2,
        2.0264558490e-1,
    ];
    let expected_deviations = vec![
        2.2708117030e-3,
    ];
    let expected_extremal_frequencies = vec![
        7.6219509360e-4,    1.2195123360e-2,    2.5152431800e-2,    3.7347543980e-2,
        4.9542654310e-2,    6.2499959020e-2,    7.4695073070e-2,    8.6890183390e-2,
        9.9847488110e-2,    1.1204259840e-1,    1.2423770870e-1,    1.3719502090e-1,
        1.4939013120e-1,    1.6158524160e-1,    1.7454254630e-1,    1.8673765660e-1,
        1.9969496130e-1,    2.1189007160e-1,    2.2408518200e-1,    2.3704248670e-1,
        2.4923759700e-1,    2.6143270730e-1,    2.7439001200e-1,    2.8658512230e-1,
        2.9878023270e-1,    3.1173753740e-1,    3.2393264770e-1,    3.3688995240e-1,
        3.4908506270e-1,    3.6128017310e-1,    3.7423747780e-1,    3.8643258810e-1,
        3.9862769840e-1,    4.1158500310e-1,    4.2378011350e-1,    4.3597522380e-1,
        4.4893252850e-1,    4.6112763880e-1,    4.7332274910e-1,    4.8475566510e-1,
        4.9542638660e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,2,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n128_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -1.3384570780e-4,    1.6450791740e-4,   -5.6053067960e-5,    4.6555956940e-5,
        -3.0926861650e-5,    1.7389851560e-5,   -1.9907944080e-5,    1.5348818120e-5,
        -1.5805908330e-5,    2.0373841830e-5,   -1.9971408620e-5,    2.3472472090e-5,
        -2.6209461790e-5,    2.5863429980e-5,   -2.8126796680e-5,    2.8162266970e-5,
        -2.7289650460e-5,    2.8372318410e-5,   -2.7595595380e-5,    2.7692447470e-5,
        -2.9117349190e-5,    2.9355371230e-5,   -3.1119358030e-5,    3.3225765950e-5,
        -3.4524397050e-5,    3.7100588090e-5,   -3.9138678400e-5,    4.0773684300e-5,
        -4.3310312320e-5,    4.5079803390e-5,   -4.7109628210e-5,    4.9880480220e-5,
        -5.2212613810e-5,    5.5423013690e-5,   -5.9233654610e-5,    6.3031351600e-5,
        -6.7959532320e-5,    7.3280119980e-5,   -7.8962548290e-5,    8.5815991040e-5,
        -9.3078160720e-5,    1.0126164120e-4,   -1.1086923040e-4,    1.2146583320e-4,
        -1.3396152640e-4,    1.4869213920e-4,   -1.6577276980e-4,    1.8643241490e-4,
        -2.1122477480e-4,    2.4120617310e-4,   -2.7838951790e-4,    3.2467336860e-4,
        -3.8349654640e-4,    4.6005629700e-4,   -5.6182080880e-4,    7.0162629710e-4,
        -9.0107874710e-4,    1.1993958610e-3,   -1.6750202050e-3,    2.5020218450e-3,
        -4.1357548910e-3,    8.1058889630e-3,   -2.2516006600e-2,    2.0264256000e-1,
    ];
    let expected_deviations = vec![
        1.2765545400e-3,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    7.8125000000e-3,    1.5625000000e-2,    2.3437500000e-2,
        3.1250000000e-2,    3.9062500000e-2,    4.6875000000e-2,    5.4687500000e-2,
        6.2500000000e-2,    7.0312500000e-2,    7.8125000000e-2,    8.5937500000e-2,
        9.3750000000e-2,    1.0156250000e-1,    1.0937500000e-1,    1.1718750000e-1,
        1.2500000000e-1,    1.3281250000e-1,    1.4453125000e-1,    1.5234375000e-1,
        1.6015625000e-1,    1.6796875000e-1,    1.7578125000e-1,    1.8359375000e-1,
        1.9140625000e-1,    1.9921875000e-1,    2.0703125000e-1,    2.1484375000e-1,
        2.2265625000e-1,    2.3046875000e-1,    2.3828125000e-1,    2.4609375000e-1,
        2.5390625000e-1,    2.6171875000e-1,    2.6953125000e-1,    2.7734375000e-1,
        2.8515625000e-1,    2.9296875000e-1,    3.0078125000e-1,    3.0859375000e-1,
        3.1640625000e-1,    3.2421875000e-1,    3.3203125000e-1,    3.3984375000e-1,
        3.4765625000e-1,    3.5546875000e-1,    3.6328125000e-1,    3.7109375000e-1,
        3.7890625000e-1,    3.8671875000e-1,    3.9453125000e-1,    4.0234375000e-1,
        4.1015625000e-1,    4.1796875000e-1,    4.2968750000e-1,    4.3750000000e-1,
        4.4531250000e-1,    4.5312500000e-1,    4.6093750000e-1,    4.6875000000e-1,
        4.7656250000e-1,    4.8437500000e-1,    4.9218750000e-1,    4.9609375000e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,2,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n128_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -1.3669385230e-4,    1.6701171990e-4,   -4.9292346380e-5,    3.1809311620e-5,
        -2.7089918150e-5,    2.9863189410e-5,   -2.6161036660e-5,    1.9389593940e-5,
        -1.9950879500e-5,    2.0865030820e-5,   -2.2415133570e-5,    2.0967254390e-5,
        -1.9804601830e-5,    2.1475243560e-5,   -2.2620575690e-5,    2.3631073420e-5,
        -2.3257121940e-5,    2.3828431950e-5,   -2.5706483940e-5,    2.6895540940e-5,
        -2.7945909100e-5,    2.8474303690e-5,   -2.9996434020e-5,    3.2055504560e-5,
        -3.3516080290e-5,    3.5032157030e-5,   -3.6508790800e-5,    3.8849168050e-5,
        -4.1353305280e-5,    4.3494044800e-5,   -4.5932338250e-5,    4.8630805400e-5,
        -5.2073955890e-5,    5.5597272880e-5,   -5.9116471680e-5,    6.3224732000e-5,
        -6.7856613900e-5,    7.3283037640e-5,   -7.9032717620e-5,    8.5329797000e-5,
        -9.2695445350e-5,    1.0109707360e-4,   -1.1080381230e-4,    1.2166592930e-4,
        -1.3422060870e-4,    1.4909214220e-4,   -1.6655726360e-4,    1.8731618180e-4,
        -2.1199331970e-4,    2.4201902850e-4,   -2.7913396480e-4,    3.2542797270e-4,
        -3.8431197750e-4,    4.6065473000e-4,   -5.6246493480e-4,    7.0240220520e-4,
        -9.0186262970e-4,    1.2003167070e-3,   -1.6759017020e-3,    2.5029275570e-3,
        -4.1368044910e-3,    8.1069376320e-3,   -2.2517070170e-2,    2.0264355840e-1,
    ];
    let expected_deviations = vec![
        1.3869783140e-3,
    ];
    let expected_extremal_frequencies = vec![
        2.6041667440e-3,    7.8125000000e-3,    1.5625000000e-2,    2.3437498140e-2,
        3.1249996270e-2,    3.9062500000e-2,    4.6875003730e-2,    5.4687507450e-2,
        6.2500007450e-2,    7.0312500000e-2,    7.8124992550e-2,    8.5937485100e-2,
        9.6354141830e-2,    1.0416663440e-1,    1.1197912690e-1,    1.1979161950e-1,
        1.2760411200e-1,    1.3541662690e-1,    1.4322914180e-1,    1.5104165670e-1,
        1.5885417160e-1,    1.6666668650e-1,    1.7447920140e-1,    1.8229171630e-1,
        1.9010423120e-1,    1.9791674610e-1,    2.0572926100e-1,    2.1354177590e-1,
        2.2135429080e-1,    2.2916680570e-1,    2.3697932060e-1,    2.4479183550e-1,
        2.5260433550e-1,    2.6041680570e-1,    2.6822927590e-1,    2.7604174610e-1,
        2.8645837310e-1,    2.9427084330e-1,    3.0208331350e-1,    3.0989578370e-1,
        3.1770825390e-1,    3.2552072410e-1,    3.3333319430e-1,    3.4114566450e-1,
        3.4895813470e-1,    3.5677060480e-1,    3.6458307500e-1,    3.7239554520e-1,
        3.8020801540e-1,    3.8802048560e-1,    3.9583295580e-1,    4.0364542600e-1,
        4.1145789620e-1,    4.1927036640e-1,    4.2708283660e-1,    4.3489530680e-1,
        4.4270777700e-1,    4.5052024720e-1,    4.5833271740e-1,    4.6614518760e-1,
        4.7395765780e-1,    4.8177012800e-1,    4.8958259820e-1,    4.9739506840e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,2,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n128_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -1.3524550010e-4,    1.6534118910e-4,   -4.6008652130e-5,    2.7858068280e-5,
        -2.2561249350e-5,    2.0600793500e-5,   -1.9851193430e-5,    1.9652361520e-5,
        -1.9758041160e-5,    2.0033568940e-5,   -2.0428891730e-5,    2.0924158890e-5,
        -2.1495112380e-5,    2.2120861100e-5,   -2.2805876140e-5,    2.3574515580e-5,
        -2.4424309230e-5,    2.5335559260e-5,   -2.6303056070e-5,    2.7354555640e-5,
        -2.8509599360e-5,    2.9740323950e-5,   -3.1073908760e-5,    3.2529576860e-5,
        -3.4095246520e-5,    3.5831704740e-5,   -3.7563673690e-5,    3.9702514190e-5,
        -4.2045692680e-5,    4.4163389250e-5,   -4.6562508940e-5,    4.9339054390e-5,
        -5.2450584920e-5,    5.5871081710e-5,   -5.9626159780e-5,    6.3784245870e-5,
        -6.8396489950e-5,    7.3533337850e-5,   -7.9288576670e-5,    8.5764273540e-5,
        -9.3076203480e-5,    1.0137583010e-4,   -1.1088080650e-4,    1.2184161460e-4,
        -1.3453551220e-4,    1.4931110490e-4,   -1.6668101310e-4,    1.8735261980e-4,
        -2.1213534640e-4,    2.4220583150e-4,   -2.7923041490e-4,    3.2547302540e-4,
        -3.8433275770e-4,    4.6074105190e-4,   -5.6256190870e-4,    7.0246879480e-4,
        -9.0188544710e-4,    1.2002564040e-3,   -1.6758944840e-3,    2.5029375680e-3,
        -4.1367644440e-3,    8.1069106240e-3,   -2.2517038510e-2,    2.0264358820e-1,
    ];
    let expected_deviations = vec![
        1.4358713520e-3,
    ];
    let expected_extremal_frequencies = vec![
        5.2083336050e-4,    7.8125000000e-3,    1.5625005590e-2,    2.3958330970e-2,
        3.1770825390e-2,    3.9583317940e-2,    4.7395810480e-2,    5.5208303030e-2,
        6.3541628420e-2,    7.1354120970e-2,    7.9166613520e-2,    8.6979106070e-2,
        9.4791598620e-2,    1.0312492400e-1,    1.1093741660e-1,    1.1874990910e-1,
        1.2656241660e-1,    1.3489586110e-1,    1.4270846550e-1,    1.5052106980e-1,
        1.5833367410e-1,    1.6614627840e-1,    1.7447972300e-1,    1.8229232730e-1,
        1.9010493160e-1,    1.9791753590e-1,    2.0573014020e-1,    2.1406358480e-1,
        2.2187618910e-1,    2.2968879340e-1,    2.3750139770e-1,    2.4531400200e-1,
        2.5364732740e-1,    2.6145970820e-1,    2.6927208900e-1,    2.7708446980e-1,
        2.8541767600e-1,    2.9323005680e-1,    3.0104243760e-1,    3.0885481830e-1,
        3.1666719910e-1,    3.2500040530e-1,    3.3281278610e-1,    3.4062516690e-1,
        3.4843754770e-1,    3.5677075390e-1,    3.6458313470e-1,    3.7239551540e-1,
        3.8020789620e-1,    3.8854110240e-1,    3.9635348320e-1,    4.0416586400e-1,
        4.1197824480e-1,    4.1979062560e-1,    4.2812383170e-1,    4.3593621250e-1,
        4.4374859330e-1,    4.5156097410e-1,    4.5937335490e-1,    4.6770656110e-1,
        4.7551894190e-1,    4.8281049730e-1,    4.9062287810e-1,    4.9687278270e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,2,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn differentiator_allpass_n128_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -1.3501262580e-4,    1.6505343960e-4,   -4.5909757320e-5,    2.7781969040e-5,
        -2.2492044080e-5,    2.0534898790e-5,   -1.9801889720e-5,    1.9627588700e-5,
        -1.9740447900e-5,    2.0020261220e-5,   -2.0424875400e-5,    2.0933832270e-5,
        -2.1513295000e-5,    2.2153300960e-5,   -2.2868596720e-5,    2.3646487530e-5,
        -2.4484299500e-5,    2.5391156670e-5,   -2.6369692930e-5,    2.7426151060e-5,
        -2.8560029020e-5,    2.9786227970e-5,   -3.1117051550e-5,    3.2554569770e-5,
        -3.4103461080e-5,    3.5781617040e-5,   -3.7600759240e-5,    3.9633378040e-5,
        -4.1694765970e-5,    4.4147265730e-5,   -4.6942044720e-5,    4.9629903510e-5,
        -5.2552244600e-5,    5.5869946660e-5,   -5.9604342820e-5,    6.3769635740e-5,
        -6.8408815420e-5,    7.3575938590e-5,   -7.9336437920e-5,    8.5795691120e-5,
        -9.3114300400e-5,    1.0145182020e-4,   -1.1095446830e-4,    1.2189036350e-4,
        -1.3456639140e-4,    1.4934127100e-4,   -1.6672819040e-4,    1.8737293430e-4,
        -2.1214934530e-4,    2.4222736830e-4,   -2.7923303420e-4,    3.2547820590e-4,
        -3.8431215220e-4,    4.6074378770e-4,   -5.6256982500e-4,    7.0240104110e-4,
        -9.0186716990e-4,    1.2002896980e-3,   -1.6759294090e-3,    2.5030111430e-3,
        -4.1368287060e-3,    8.1069217990e-3,   -2.2516995670e-2,    2.0264351370e-1,
    ];
    let expected_deviations = vec![
        1.4337963660e-3,
    ];
    let expected_extremal_frequencies = vec![
        4.8828125000e-4,    7.8125000000e-3,    1.5625000000e-2,    2.3925781250e-2,
        3.1738281250e-2,    3.9550781250e-2,    4.7363281250e-2,    5.5664062500e-2,
        6.3476562500e-2,    7.1289062500e-2,    7.9101562500e-2,    8.6914062500e-2,
        9.5214843750e-2,    1.0302734380e-1,    1.1083984380e-1,    1.1865234380e-1,
        1.2695312500e-1,    1.3476562500e-1,    1.4257812500e-1,    1.5039062500e-1,
        1.5820312500e-1,    1.6650390620e-1,    1.7431640620e-1,    1.8212890620e-1,
        1.8994140620e-1,    1.9824218750e-1,    2.0605468750e-1,    2.1386718750e-1,
        2.2167968750e-1,    2.2949218750e-1,    2.3779296880e-1,    2.4560546880e-1,
        2.5341796880e-1,    2.6123046880e-1,    2.6953125000e-1,    2.7734375000e-1,
        2.8515625000e-1,    2.9296875000e-1,    3.0126953120e-1,    3.0908203120e-1,
        3.1689453120e-1,    3.2470703120e-1,    3.3300781250e-1,    3.4082031250e-1,
        3.4863281250e-1,    3.5644531250e-1,    3.6425781250e-1,    3.7255859380e-1,
        3.8037109380e-1,    3.8818359380e-1,    3.9599609380e-1,    4.0429687500e-1,
        4.1210937500e-1,    4.1992187500e-1,    4.2773437500e-1,    4.3603515620e-1,
        4.4384765620e-1,    4.5166015620e-1,    4.5947265620e-1,    4.6728515620e-1,
        4.7509765620e-1,    4.8291015620e-1,    4.9023437500e-1,    4.9707031250e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n3_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        7.3218137030e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.3926917310e-1,    1.3926917310e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n3_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        7.3218137030e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.3926917310e-1,    1.3926917310e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n3_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        7.2755061090e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.4471299950e-1,    1.4471299950e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.6666668060e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n3_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        7.2729706760e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.4501102270e-1,    1.4501102270e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.6249998810e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n3_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        7.2703547780e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.4531852300e-1,    1.4531852300e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.5555557010e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n3_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        7.2700537740e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.4535391330e-1,    1.4535391330e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.5405406950e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n4_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        5.4617356510e-2,    1.6835704450e-3,
    ];
    let expected_deviations = vec![
        1.0586754230e-1,    1.0586754230e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n4_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        5.4617356510e-2,    1.6835704450e-3,
    ];
    let expected_deviations = vec![
        1.0586754230e-1,    1.0586754230e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n4_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        5.4617356510e-2,    1.6835704450e-3,
    ];
    let expected_deviations = vec![
        1.0586754230e-1,    1.0586754230e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n4_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        5.4617356510e-2,    1.6835704450e-3,
    ];
    let expected_deviations = vec![
        1.0586754230e-1,    1.0586754230e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n4_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        5.4617356510e-2,    1.6835704450e-3,
    ];
    let expected_deviations = vec![
        1.0586754230e-1,    1.0586754230e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n4_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        5.4617356510e-2,    1.6835704450e-3,
    ];
    let expected_deviations = vec![
        1.0586754230e-1,    1.0586754230e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n11_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -1.2871638870e-2,    1.0790772740e-2,    1.1589542030e-2,    2.4629242720e-2,
        1.3869040650e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.1184291690e-2,    2.1184291690e-2,
    ];
    let expected_extremal_frequencies = vec![
        5.0000000750e-2,    1.0000000150e-1,    2.0000000300e-1,    2.5000000000e-1,
        3.5000002380e-1,    4.0000003580e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n11_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -1.1146811770e-2,    7.1909027170e-3,    1.6088150440e-2,    2.1278087050e-2,
        1.5543816610e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.1989312020e-2,    2.1989312020e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.3333335070e-2,    1.0000000150e-1,    2.0000000300e-1,    2.6666668060e-1,
        3.6666667460e-1,    4.3333333730e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n11_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -9.1291666030e-3,    1.5109289670e-3,    2.0531717690e-2,    1.9140319900e-2,
        1.6996905210e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.7820391580e-2,    2.7820391580e-2,
    ];
    let expected_extremal_frequencies = vec![
        6.6666668280e-3,    6.6666670140e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.5999996070e-1,    3.7333318590e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n11_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -9.1284401710e-3,    1.5254402530e-3,    2.0504768940e-2,    1.9144531340e-2,
        1.7015416180e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.7864629400e-2,    2.7864629400e-2,
    ];
    let expected_extremal_frequencies = vec![
        6.2500000930e-3,    6.8750008940e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.6249995830e-1,    3.6874985690e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n11_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        -9.0701393780e-3,    1.4136815440e-3,    2.0594067870e-2,    1.9109249110e-2,
        1.7031546680e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.7961436660e-2,    2.7961436660e-2,
    ];
    let expected_extremal_frequencies = vec![
        2.7777778450e-3,    6.6666670140e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.6388904450e-1,    3.7222266200e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n11_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        -9.0605141600e-3,    1.3903932410e-3,    2.0611181860e-2,    1.9106771800e-2,
        1.7029643060e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.7979077770e-2,    2.7979077770e-2,
    ];
    let expected_extremal_frequencies = vec![
        2.7027027680e-3,    6.7567557100e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.6216211920e-1,    3.7027063970e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n12_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -8.7079685180e-3,    1.1058207600e-3,    1.0644664060e-2,    1.9422199580e-2,
        1.9391840320e-2,    8.1296823920e-3,
    ];
    let expected_deviations = vec![
        1.4658345840e-2,    1.4658345840e-2,
    ];
    let expected_extremal_frequencies = vec![
        4.1666667910e-2,    1.0000000150e-1,    2.0000000300e-1,    2.4166667460e-1,
        3.2499998810e-1,    4.0833330150e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n12_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -4.8890295440e-3,   -6.6093639470e-3,    1.4909936120e-2,    1.8018078060e-2,
        2.0797237750e-2,    8.0609619620e-3,
    ];
    let expected_deviations = vec![
        2.2696949540e-2,    2.2696949540e-2,
    ];
    let expected_extremal_frequencies = vec![
        2.7777777980e-2,    5.5555555970e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.5555557010e-1,    3.6666673420e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n12_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -2.8328415940e-3,   -9.5318071540e-3,    1.4959503900e-2,    1.9381238150e-2,
        2.0431598650e-2,    8.5446089510e-3,
    ];
    let expected_deviations = vec![
        2.8328439220e-2,    2.8328439220e-2,
    ];
    let expected_extremal_frequencies = vec![
        5.5555556900e-3,    6.6666670140e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.6666671040e-1,    3.7777811290e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n12_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -2.8204461560e-3,   -9.5489565280e-3,    1.4952288940e-2,    1.9402571020e-2,
        2.0427877080e-2,    8.5298605260e-3,
    ];
    let expected_deviations = vec![
        2.8352476660e-2,    2.8352476660e-2,
    ];
    let expected_extremal_frequencies = vec![
        5.2083334890e-3,    6.7708328370e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.6249998810e-1,    3.7708354000e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n12_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        -2.8030222750e-3,   -9.5616094770e-3,    1.4949640260e-2,    1.9398109990e-2,
        2.0440077410e-2,    8.5375979540e-3,
    ];
    let expected_deviations = vec![
        2.8425181280e-2,    2.8425181280e-2,
    ];
    let expected_extremal_frequencies = vec![
        2.3148148320e-3,    6.7129611970e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.6481488350e-1,    3.7592557070e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n12_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        -2.8015940920e-3,   -9.5635568720e-3,    1.4949146660e-2,    1.9399059940e-2,
        2.0440479740e-2,    8.5372887550e-3,
    ];
    let expected_deviations = vec![
        2.8430484240e-2,    2.8430484240e-2,
    ];
    let expected_extremal_frequencies = vec![
        2.2522523070e-3,    6.7567549650e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.6531529430e-1,    3.7567558880e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n16_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        3.9568273350e-3,   -5.8662290680e-3,   -5.5477153510e-3,   -6.5453024580e-4,
        1.0071557020e-2,    2.0074304190e-2,    1.9965482880e-2,    8.7783746420e-3,
    ];
    let expected_deviations = vec![
        1.2228463780e-2,    1.2228463780e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.1250000000e-2,    6.2500000000e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.3125000300e-1,    2.9374998810e-1,    3.5624998810e-1,    4.1874998810e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n16_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        3.9711007850e-3,   -5.5465977640e-3,   -5.9567363930e-3,   -6.4344285060e-4,
        1.0238593440e-2,    1.9934264940e-2,    2.0218316470e-2,    8.5867755120e-3,
    ];
    let expected_deviations = vec![
        1.2280519120e-2,    1.2280519120e-2,
    ];
    let expected_extremal_frequencies = vec![
        2.0833333950e-2,    6.2500000000e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.4166665970e-1,    3.0416667460e-1,    3.6666670440e-1,    4.2916673420e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n16_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        4.1756555440e-3,   -5.6630545300e-3,   -5.9292213990e-3,   -8.0055976290e-4,
        1.0069048030e-2,    2.0017176870e-2,    2.0419236270e-2,    8.6522772910e-3,
    ];
    let expected_deviations = vec![
        1.3057715260e-2,    1.3057715260e-2,
    ];
    let expected_extremal_frequencies = vec![
        4.1666668840e-3,    7.5000010430e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.3333330450e-1,    2.9583325980e-1,    3.6249986290e-1,    4.3333312870e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n16_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        4.1769058440e-3,   -5.6633539500e-3,   -5.9263650330e-3,   -8.0783059820e-4,
        1.0070610790e-2,    2.0017150790e-2,    2.0419647920e-2,    8.6580142380e-3,
    ];
    let expected_deviations = vec![
        1.3073607350e-2,    1.3073607350e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    7.8125000000e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.3515625300e-1,    2.9765623810e-1,    3.6406248810e-1,    4.3046873810e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n16_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        4.1780886240e-3,   -5.6623909620e-3,   -5.9287217450e-3,   -8.0700404940e-4,
        1.0068306700e-2,    2.0016912370e-2,    2.0424464720e-2,    8.6542107160e-3,
    ];
    let expected_deviations = vec![
        1.3080768290e-2,    1.3080768290e-2,
    ];
    let expected_extremal_frequencies = vec![
        1.7361111240e-3,    7.6388917860e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.3472209270e-1,    2.9722186920e-1,    3.6492994430e-1,    4.3263801930e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n16_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        4.1776238940e-3,   -5.6607648730e-3,   -5.9308032510e-3,   -8.0699007960e-4,
        1.0070761670e-2,    2.0015483720e-2,    2.0423863080e-2,    8.6552761500e-3,
    ];
    let expected_deviations = vec![
        1.3076831590e-2,    1.3076831590e-2,
    ];
    let expected_extremal_frequencies = vec![
        1.6891892300e-3,    7.7702686190e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.3378391560e-1,    2.9628413920e-1,    3.6385196450e-1,    4.3141978980e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n47_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        3.3780422650e-6,    3.7776429960e-5,    9.0408957480e-6,   -7.7828895880e-5,
        -1.5329929010e-4,   -6.0599708380e-5,    2.3358636830e-4,    4.8572686500e-4,
        2.8997071790e-4,   -4.5687903180e-4,   -1.1932844060e-3,   -9.7370485310e-4,
        5.6714389940e-4,    2.4169646200e-3,    2.6015937330e-3,   -5.9663085270e-5,
        -4.2263213550e-3,   -6.2058456240e-3,   -2.5785781910e-3,    6.4932089300e-3,
        1.6207765790e-2,    1.9843494520e-2,    1.3631141740e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.5176841520e-5,    2.5176841520e-5,
    ];
    let expected_extremal_frequencies = vec![
        1.0869565420e-2,    2.1739130840e-2,    3.2608695330e-2,    5.4347828030e-2,
        6.5217390660e-2,    7.6086953280e-2,    8.6956515910e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.1086956560e-1,    2.2173912820e-1,    2.3260869090e-1,
        2.5434783100e-1,    2.7608695630e-1,    2.8695651890e-1,    3.0869564410e-1,
        3.3043476940e-1,    3.5217389460e-1,    3.7391301990e-1,    3.9565214510e-1,
        4.1739127040e-1,    4.3913039570e-1,    4.6086952090e-1,    4.8260864620e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n47_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        3.8638017940e-6,    4.1283110480e-5,    6.4404648580e-6,   -8.6618580100e-5,
        -1.5726806300e-4,   -5.2957049770e-5,    2.5300297420e-4,    4.9589737320e-4,
        2.7286467960e-4,   -4.9281888640e-4,   -1.2116474100e-3,   -9.4753835580e-4,
        6.2243628780e-4,    2.4488591590e-3,    2.5712242350e-3,   -1.3413256970e-4,
        -4.2773643510e-3,   -6.1789290050e-3,   -2.4896557440e-3,    6.5644914280e-3,
        1.6191784290e-2,    1.9749846310e-2,    1.3544066810e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        3.4387441700e-5,    3.4387441700e-5,
    ];
    let expected_extremal_frequencies = vec![
        7.2463769470e-3,    2.1739130840e-2,    3.6231882870e-2,    5.0724633040e-2,
        7.2463758290e-2,    7.9710133370e-2,    8.6956508460e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0724637810e-1,    2.1449275310e-1,    2.2898550330e-1,
        2.5072464350e-1,    2.7246376870e-1,    2.9420289400e-1,    3.0869564410e-1,
        3.3043476940e-1,    3.5217389460e-1,    3.8115939500e-1,    4.0289852020e-1,
        4.2463764550e-1,    4.4637677070e-1,    4.6811589600e-1,    4.8985502120e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n47_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -9.8024966060e-6,    5.3300449510e-5,    3.5062857930e-5,   -6.7469882200e-5,
        -1.8622213970e-4,   -1.2986674850e-4,    1.9243150020e-4,    5.3814315470e-4,
        4.2800031950e-4,   -3.4329880150e-4,   -1.2349708700e-3,   -1.1888371080e-3,
        3.3579554290e-4,    2.3930377790e-3,    2.8649617450e-3,    3.0806974970e-4,
        -4.0770778430e-3,   -6.4510968510e-3,   -3.0533580580e-3,    6.1866827310e-3,
        1.6356626530e-2,    2.0349450410e-2,    1.4073755590e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.4234292950e-5,    6.4234292950e-5,
    ];
    let expected_extremal_frequencies = vec![
        1.4492754130e-3,    1.8840583040e-2,    3.7681166080e-2,    5.6521750990e-2,
        7.2463758290e-2,    8.6956478660e-2,    9.7101382910e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0434781910e-1,    2.1594199540e-1,    2.3333325980e-1,
        2.5217381120e-1,    2.7246382830e-1,    2.9275384550e-1,    3.1449314950e-1,
        3.3623245360e-1,    3.5797175770e-1,    3.7971106170e-1,    4.0145036580e-1,
        4.2318966980e-1,    4.4492897390e-1,    4.6666827800e-1,    4.8840758200e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n47_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -1.0020128680e-5,    5.3373125410e-5,    3.5544755520e-5,   -6.6968830650e-5,
        -1.8649989220e-4,   -1.3113467140e-4,    1.9110032010e-4,    5.3839344760e-4,
        4.3053951230e-4,   -3.4034726560e-4,   -1.2347097510e-3,   -1.1926747390e-3,
        3.3051642820e-4,    2.3913076150e-3,    2.8695818040e-3,    3.1592836600e-4,
        -4.0729884060e-3,   -6.4553730190e-3,   -3.0631423000e-3,    6.1797928070e-3,
        1.6359191390e-2,    2.0359799270e-2,    1.4083035290e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.4505635240e-5,    6.4505635240e-5,
    ];
    let expected_extremal_frequencies = vec![
        1.3586956770e-3,    1.9021736460e-2,    3.8043472920e-2,    5.7065207510e-2,
        7.3369555180e-2,    8.6956508460e-2,    9.6467375760e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0407611130e-1,    2.1630443630e-1,    2.3260886970e-1,
        2.5163069370e-1,    2.7201101180e-1,    2.9239133000e-1,    3.1413033600e-1,
        3.3586934210e-1,    3.5760834810e-1,    3.7934735420e-1,    4.0108636020e-1,
        4.2282536630e-1,    4.4456437230e-1,    4.6766206620e-1,    4.8940107230e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n47_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        -1.0021499290e-5,    5.3422845670e-5,    3.5484223190e-5,   -6.7103967010e-5,
        -1.8649736010e-4,   -1.3094193130e-4,    1.9140265070e-4,    5.3843902420e-4,
        4.3010435180e-4,   -3.4095824230e-4,   -1.2348762250e-3,   -1.1920080290e-3,
        3.3155526030e-4,    2.3917737420e-3,    2.8687738810e-3,    3.1438900620e-4,
        -4.0738731620e-3,   -6.4546037470e-3,   -3.0612582340e-3,    6.1811618510e-3,
        1.6358744350e-2,    2.0357813690e-2,    1.4081222010e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.4682768420e-5,    6.4682768420e-5,
    ];
    let expected_extremal_frequencies = vec![
        6.0386472610e-4,    1.9323669370e-2,    3.8043472920e-2,    5.6159447880e-2,
        7.3067627850e-2,    8.6956456300e-2,    9.6618250010e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0483095940e-1,    2.1630448100e-1,    2.3321282860e-1,
        2.5193274020e-1,    2.7185994390e-1,    2.9299485680e-1,    3.1412976980e-1,
        3.3586853740e-1,    3.5760730500e-1,    3.7934607270e-1,    4.0108484030e-1,
        4.2282360790e-1,    4.4516623020e-1,    4.6690499780e-1,    4.8924762010e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n47_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        -1.0003680470e-5,    5.3440406190e-5,    3.5435114110e-5,   -6.7189357650e-5,
        -1.8650670240e-4,   -1.3080264030e-4,    1.9160113880e-4,    5.3847138770e-4,
        4.2984008910e-4,   -3.4134558520e-4,   -1.2349882160e-3,   -1.1916186190e-3,
        3.3218716270e-4,    2.3920456880e-3,    2.8683333660e-3,    3.1349738130e-4,
        -4.0744054130e-3,   -6.4542032780e-3,   -3.0601697510e-3,    6.1819767580e-3,
        1.6358500350e-2,    2.0356668160e-2,    1.4080177990e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.4694904720e-5,    6.4694904720e-5,
    ];
    let expected_extremal_frequencies = vec![
        5.8754405470e-4,    1.9388953220e-2,    3.8190364840e-2,    5.6404270230e-2,
        7.2855539620e-2,    8.6956627670e-2,    9.6357353030e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0470030610e-1,    2.1645106380e-1,    2.3290212450e-1,
        2.5170338150e-1,    2.7226772900e-1,    2.9283207650e-1,    3.1398397680e-1,
        3.3572342990e-1,    3.5746288300e-1,    3.7920233610e-1,    4.0094178920e-1,
        4.2326879500e-1,    4.4500824810e-1,    4.6674770120e-1,    4.8907470700e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n48_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -4.7372924430e-7,    2.7429397960e-5,    2.1437670510e-5,   -4.2706044040e-5,
        -1.2737193900e-4,   -1.0738690620e-4,    1.1231253300e-4,    4.0206403360e-4,
        4.0633173190e-4,   -1.3058981860e-4,   -9.3206122980e-4,   -1.1727912350e-3,
        -1.6896764280e-4,    1.6875182050e-3,    2.7631367560e-3,    1.3717541010e-3,
        -2.3678285070e-3,   -5.7336003520e-3,   -4.9209743740e-3,    1.7973994840e-3,
        1.1774261480e-2,    1.8892288210e-2,    1.7594654110e-2,    7.1716867390e-3,
    ];
    let expected_deviations = vec![
        1.7217884310e-5,    1.7217884310e-5,
    ];
    let expected_extremal_frequencies = vec![
        1.0416666980e-2,    2.0833333950e-2,    3.1250000000e-2,    5.2083335820e-2,
        6.2500000000e-2,    7.2916664180e-2,    8.3333328370e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.1041667460e-1,    2.2083334620e-1,    2.3125001790e-1,
        2.5208336110e-1,    2.7291667460e-1,    2.8333333130e-1,    3.0416664480e-1,
        3.2499995830e-1,    3.4583327170e-1,    3.6666658520e-1,    3.9791655540e-1,
        4.1874986890e-1,    4.3958318230e-1,    4.6041649580e-1,    4.8124980930e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n48_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -1.2381036870e-5,    3.4363823940e-5,    4.6917502910e-5,   -1.6924099330e-5,
        -1.4261719480e-4,   -1.7795184980e-4,    3.5082906830e-5,    4.0933344280e-4,
        5.3983536780e-4,    4.6477187420e-5,   -8.8145653720e-4,   -1.3545525730e-3,
        -4.8481428530e-4,    1.5079164880e-3,    2.9341231570e-3,    1.8190767150e-3,
        -2.0029440060e-3,   -5.8035776020e-3,   -5.4315831510e-3,    1.2453999370e-3,
        1.1664247140e-2,    1.9356213510e-2,    1.8265144900e-2,    7.4842832980e-3,
    ];
    let expected_deviations = vec![
        4.1005809180e-5,    4.1005809180e-5,
    ];
    let expected_extremal_frequencies = vec![
        6.9444444960e-3,    2.0833333950e-2,    3.4722223880e-2,    5.5555555970e-2,
        6.9444447760e-2,    8.3333343270e-2,    9.0277791020e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0694445070e-1,    2.1388889850e-1,    2.3472224180e-1,
        2.4861113730e-1,    2.6944446560e-1,    2.9027780890e-1,    3.1111115220e-1,
        3.3194449540e-1,    3.5277783870e-1,    3.7361118200e-1,    3.9444452520e-1,
        4.1527786850e-1,    4.3611121180e-1,    4.5694455500e-1,    4.7777789830e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n48_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -1.9411403630e-5,    3.8529124140e-5,    5.9359579610e-5,   -4.8269357650e-6,
        -1.4853305770e-4,   -2.0900717940e-4,   -2.8635258790e-7,    4.0923690540e-4,
        5.9567129940e-4,    1.2482755120e-4,   -8.5482839490e-4,   -1.4284276400e-3,
        -6.1988167000e-4,    1.4267512600e-3,    3.0012044590e-3,    2.0066956060e-3,
        -1.8455632960e-3,   -5.8278869840e-3,   -5.6439014150e-3,    1.0119285430e-3,
        1.1614916850e-2,    1.9549055020e-2,    1.8546413630e-2,    7.6157674190e-3,
    ];
    let expected_deviations = vec![
        5.4959335100e-5,    5.4959335100e-5,
    ];
    let expected_extremal_frequencies = vec![
        1.3888889230e-3,    1.9444445150e-2,    3.8888890300e-2,    5.6944444780e-2,
        7.3611140250e-2,    8.7500065570e-2,    9.5833420750e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0416668060e-1,    2.1527782080e-1,    2.3194453120e-1,
        2.5000011920e-1,    2.6944440600e-1,    2.8888869290e-1,    3.0972185730e-1,
        3.3055502180e-1,    3.5138818620e-1,    3.7222135070e-1,    3.9444339280e-1,
        4.1527655720e-1,    4.3610972170e-1,    4.5694288610e-1,    4.7916492820e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n48_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -1.9597717260e-5,    3.8517915530e-5,    5.9758349380e-5,   -4.2819519880e-6,
        -1.4847636340e-4,   -2.0999368280e-4,   -1.8360733520e-6,    4.0873314720e-4,
        5.9746671470e-4,    1.2801354750e-4,   -8.5308152480e-4,   -1.4306863300e-3,
        -6.2525714750e-4,    1.4228442450e-3,    3.0030780470e-3,    2.0140132400e-3,
        -1.8387769810e-3,   -5.8281579990e-3,   -5.6521021760e-3,    1.0023349900e-3,
        1.1612522420e-2,    1.9556581970e-2,    1.8557736650e-2,    7.6210834090e-3,
    ];
    let expected_deviations = vec![
        5.5136290030e-5,    5.5136290030e-5,
    ];
    let expected_extremal_frequencies = vec![
        1.3020833720e-3,    1.9531250000e-2,    3.9062496270e-2,    5.7291645560e-2,
        7.2916656730e-2,    8.7239600720e-2,    9.6354201440e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0390623810e-1,    2.1562494340e-1,    2.3124988380e-1,
        2.4947898090e-1,    2.6901036500e-1,    2.8984385730e-1,    3.0937525630e-1,
        3.3020874860e-1,    3.5234433410e-1,    3.7317782640e-1,    3.9401131870e-1,
        4.1484481100e-1,    4.3567830320e-1,    4.5781388880e-1,    4.7864738110e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n48_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        -1.9607539800e-5,    3.8557838710e-5,    5.9741738370e-5,   -4.3236141210e-6,
        -1.4854408800e-4,   -2.0996737290e-4,   -1.6775447880e-6,    4.0886708300e-4,
        5.9738586420e-4,    1.2774730570e-4,   -8.5336092160e-4,   -1.4306118250e-3,
        -6.2480068300e-4,    1.4233004770e-3,    3.0030566270e-3,    2.0133992660e-3,
        -1.8394505610e-3,   -5.8282492680e-3,   -5.6514292960e-3,    1.0032113640e-3,
        1.1612805540e-2,    1.9555950540e-2,    1.8556728960e-2,    7.6206177470e-3,
    ];
    let expected_deviations = vec![
        5.5293327020e-5,    5.5293327020e-5,
    ];
    let expected_extremal_frequencies = vec![
        5.7870370800e-4,    1.9675919790e-2,    3.8773152980e-2,    5.6713014840e-2,
        7.3495395480e-2,    8.7384231390e-2,    9.6643455330e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0462961490e-1,    2.1562494340e-1,    2.3182858530e-1,
        2.4976833160e-1,    2.6944419740e-1,    2.8969874980e-1,    3.0995330210e-1,
        3.3078655600e-1,    3.5161980990e-1,    3.7303176520e-1,    3.9386501910e-1,
        4.1527697440e-1,    4.3611022830e-1,    4.5752218370e-1,    4.7893413900e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n48_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        -1.9583525500e-5,    3.8576145020e-5,    5.9685986340e-5,   -4.4339612940e-6,
        -1.4858381470e-4,   -2.0980898990e-4,   -1.3837416190e-6,    4.0900381280e-4,
        5.9709628110e-4,    1.2714328480e-4,   -8.5372210010e-4,   -1.4302394120e-3,
        -6.2382477340e-4,    1.4240605520e-3,    3.0027679170e-3,    2.0120788830e-3,
        -1.8407322930e-3,   -5.8282502000e-3,   -5.6499475610e-3,    1.0049948470e-3,
        1.1613273060e-2,    1.9554585220e-2,    1.8554644660e-2,    7.6196230950e-3,
    ];
    let expected_deviations = vec![
        5.5293152400e-5,    5.5293152400e-5,
    ];
    let expected_extremal_frequencies = vec![
        5.6306307670e-4,    1.9707201050e-2,    3.8851335640e-2,    5.6869342920e-2,
        7.3198162020e-2,    8.7274730210e-2,    9.6846796570e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0450456440e-1,    2.1576596800e-1,    2.3153193300e-1,
        2.5011324880e-1,    2.6925712820e-1,    2.8952711820e-1,    3.0979710820e-1,
        3.3063015340e-1,    3.5146319870e-1,    3.7285929920e-1,    3.9369234440e-1,
        4.1508844490e-1,    4.3648454550e-1,    4.5731759070e-1,    4.7871369120e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,2,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n81_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -1.3977419220e-8,   -1.5797246530e-7,   -1.8694585440e-7,    3.4349599790e-7,
        1.3681415110e-6,    1.4259968570e-6,   -1.4460648570e-6,   -6.3986522040e-6,
        -7.3840919870e-6,    2.6470968350e-6,    2.0874995240e-5,    2.8297992690e-5,
        2.6787611200e-6,   -5.1370341680e-5,   -8.5073203080e-5,   -3.6427038140e-5,
        9.6496543850e-5,    2.0853633760e-4,    1.4714139980e-4,   -1.2832711220e-4,
        -4.2671099070e-4,   -4.1590799810e-4,    6.6808934210e-5,    7.3536334090e-4,
        9.4746972900e-4,    2.4535370180e-4,   -1.0532360760e-3,   -1.8439613050e-3,
        -1.0664198780e-3,    1.1617328270e-3,    3.1741899440e-3,    2.8210629240e-3,
        -5.7383114470e-4,   -4.9658883360e-3,   -6.4643267540e-3,   -2.1551474930e-3,
        7.1133146990e-3,    1.6381399710e-2,    1.9443271680e-2,    1.3152945790e-2,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        3.8582960830e-8,    3.8582960830e-8,
    ];
    let expected_extremal_frequencies = vec![
        6.2500000930e-3,    1.2500000190e-2,    2.5000000370e-2,    3.1250000000e-2,
        4.3750002980e-2,    5.0000004470e-2,    5.6250005960e-2,    6.2500007450e-2,
        6.8750008940e-2,    7.5000010430e-2,    8.1250011920e-2,    8.7500013410e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0624999700e-1,    2.1249999110e-1,
        2.1874998510e-1,    2.2499997910e-1,    2.3124997320e-1,    2.4374996130e-1,
        2.4999995530e-1,    2.6249995830e-1,    2.7499994640e-1,    2.8749993440e-1,
        2.9999992250e-1,    3.0624991660e-1,    3.1874990460e-1,    3.3124989270e-1,
        3.4374988080e-1,    3.5624986890e-1,    3.6874985690e-1,    3.8124984500e-1,
        3.9374983310e-1,    4.0624982120e-1,    4.1874980930e-1,    4.3124979730e-1,
        4.4374978540e-1,    4.5624977350e-1,    4.6874976160e-1,    4.8124974970e-1,
        4.9374973770e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
81,2,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n81_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -3.0404567750e-8,   -2.7876092190e-7,   -2.6505875890e-7,    5.9299549090e-7,
        2.0542952370e-6,    1.9568396970e-6,   -2.1944035780e-6,   -8.7903936220e-6,
        -9.6697867780e-6,    3.6870526400e-6,    2.6667708880e-5,    3.5230827050e-5,
        3.2654697860e-6,   -6.1693332100e-5,   -1.0087911510e-4,   -4.3411066140e-5,
        1.0969155120e-4,    2.3672735550e-4,    1.6772783420e-4,   -1.3802309700e-4,
        -4.6650628790e-4,   -4.5654035060e-4,    6.2481864010e-5,    7.7872304250e-4,
        1.0086419060e-3,    2.7395217330e-4,   -1.0855493600e-3,   -1.9161519590e-3,
        -1.1223531330e-3,    1.1675813000e-3,    3.2386928800e-3,    2.8949792030e-3,
        -5.4568436460e-4,   -5.0018867480e-3,   -6.5353186800e-3,   -2.2104692180e-3,
        7.1082320060e-3,    1.6425075010e-2,    1.9505236300e-2,    1.3195858340e-2,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.1637544720e-8,    8.1637544720e-8,
    ];
    let expected_extremal_frequencies = vec![
        4.1666668840e-3,    1.2500001120e-2,    2.5000000370e-2,    3.7500001490e-2,
        4.1666667910e-2,    5.4166667160e-2,    5.8333333580e-2,    6.6666670140e-2,
        7.5000010430e-2,    8.3333350720e-2,    8.7500020860e-2,    9.1666691010e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0416666570e-1,    2.0833332840e-1,
        2.1249999110e-1,    2.2083331640e-1,    2.2916664180e-1,    2.4166662990e-1,
        2.4999995530e-1,    2.6249995830e-1,    2.7083328370e-1,    2.8333327170e-1,
        2.9583325980e-1,    3.0833324790e-1,    3.2083323600e-1,    3.2916656140e-1,
        3.4166654940e-1,    3.5416653750e-1,    3.6666652560e-1,    3.7916651370e-1,
        3.9166650180e-1,    4.0416648980e-1,    4.1666647790e-1,    4.2916646600e-1,
        4.4166645410e-1,    4.5416644220e-1,    4.7083309290e-1,    4.8333308100e-1,
        4.9583306910e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
81,2,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n81_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -1.1582053360e-7,   -1.7000526500e-7,   -3.0099194250e-8,    4.4776754750e-7,
        1.1987876860e-6,    1.4489738760e-6,   -4.0692577840e-7,   -5.0070684670e-6,
        -8.7341140900e-6,   -3.7196184620e-6,    1.3874301660e-5,    3.1930358090e-5,
        2.3552491260e-5,   -2.6494350100e-5,   -8.7357155280e-5,   -8.5190491520e-5,
        2.9049626390e-5,    1.9308537590e-4,    2.3443752430e-4,    1.8462647860e-5,
        -3.5597258830e-4,   -5.3601484980e-4,   -1.9768948550e-4,    5.4961279970e-4,
        1.0636926160e-3,    6.4662902150e-4,   -6.8416428990e-4,   -1.8838160900e-3,
        -1.5811605840e-3,    5.6058238260e-4,    3.0438685790e-3,    3.3727025150e-3,
        2.5679171090e-4,   -4.5839212830e-3,   -6.9364197550e-3,   -3.1437827270e-3,
        6.4500123260e-3,    1.6655061390e-2,    2.0460641010e-2,    1.4051191510e-2,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.1381889790e-8,    6.1381889790e-8,
    ];
    let expected_extremal_frequencies = vec![
        2.4999999440e-3,    3.3333334140e-3,    7.5000007640e-3,    1.5000001530e-2,
        2.2499995310e-2,    3.0833320690e-2,    3.7499982860e-2,    5.7499963790e-2,
        6.9999955590e-2,    8.1666611140e-2,    8.9999936520e-2,    9.4999931750e-2,
        9.7499929370e-2,    1.0000000150e-1,    2.0250000060e-1,    2.0749999580e-1,
        2.1499998870e-1,    2.2499997910e-1,    2.3583330210e-1,    2.4666662510e-1,
        2.5749996300e-1,    2.6916661860e-1,    2.8083327410e-1,    2.9166659710e-1,
        3.0333325270e-1,    3.1499990820e-1,    3.2666656370e-1,    3.3833321930e-1,
        3.4999987480e-1,    3.6249986290e-1,    3.7416651840e-1,    3.8583317400e-1,
        3.9833316210e-1,    4.0999981760e-1,    4.2166647320e-1,    4.3416646120e-1,
        4.4583311680e-1,    4.5749977230e-1,    4.6999976040e-1,    4.8166641590e-1,
        4.9416640400e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
81,2,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n81_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -1.9891747630e-7,   -3.6506995120e-7,    8.2796269400e-8,    1.4761961890e-6,
        2.5790429850e-6,    6.4272467170e-7,   -5.4804213510e-6,   -1.1125903260e-5,
        -6.6280249480e-6,    1.2783100830e-5,    3.4604549000e-5,    3.0740386140e-5,
        -1.6845144270e-5,   -8.3029626690e-5,   -9.8671611340e-5,   -6.8900008050e-6,
        1.5688681740e-4,    2.4699984350e-4,    1.1296512090e-4,   -2.2629700830e-4,
        -5.0778570580e-4,   -3.9042747810e-4,    2.0444639090e-4,    8.7620737030e-4,
        9.5083180350e-4,    7.6184864160e-5,   -1.2644130040e-3,   -1.8985919890e-3,
        -8.8420079560e-4,    1.4423964310e-3,    3.2984889110e-3,    2.6510825850e-3,
        -9.1011985210e-4,   -5.1686055960e-3,   -6.3328719700e-3,   -1.7880480740e-3,
        7.3905056340e-3,    1.6309767960e-2,    1.9076405090e-2,    1.2817780490e-2,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.1385527770e-8,    6.1385527770e-8,
    ];
    let expected_extremal_frequencies = vec![
        3.1250000470e-3,    6.2500005590e-3,    8.5937511180e-3,    1.4843752610e-2,
        2.1093746650e-2,    2.8124989940e-2,    4.7656223180e-2,    6.0937460510e-2,
        7.1093738080e-2,    8.0468773840e-2,    8.9062556620e-2,    9.2968821530e-2,
        9.8437592390e-2,    2.0078125600e-1,    2.0468752090e-1,    2.1093754470e-1,
        2.1796882150e-1,    2.2656260430e-1,    2.3593764010e-1,    2.4609392880e-1,
        2.5703132150e-1,    2.6796865460e-1,    2.7890598770e-1,    2.9062455890e-1,
        3.0234313010e-1,    3.1406170130e-1,    3.2578027250e-1,    3.3749884370e-1,
        3.4921741490e-1,    3.6093598600e-1,    3.7343579530e-1,    3.8515436650e-1,
        3.9687293770e-1,    4.0937274690e-1,    4.2187255620e-1,    4.3359112740e-1,
        4.4530969860e-1,    4.5780950780e-1,    4.7030931710e-1,    4.8202788830e-1,
        4.9374645950e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,2,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n82_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        6.0845593450e-8,   -6.8900341430e-8,   -3.8228799330e-7,   -3.4852246240e-7,
        7.6823386050e-7,    2.4391972600e-6,    1.9441440600e-6,   -3.2055138490e-6,
        -1.0097067390e-5,   -8.9512113850e-6,    7.9648261820e-6,    3.1459479940e-5,
        3.2926367570e-5,   -1.0622148690e-5,   -7.8199664130e-5,   -9.8321252150e-5,
        -7.7032746050e-6,    1.5862572760e-4,    2.4503067830e-4,    9.3695416580e-5,
        -2.6115684890e-4,   -5.2178383340e-4,   -3.3399759560e-4,    3.2836454920e-4,
        9.6544891130e-4,    8.5986993510e-4,   -2.2499542680e-4,   -1.5664760720e-3,
        -1.8457232510e-3,   -3.0472222720e-4,    2.2221722170e-3,    3.5299872980e-3,
        1.7507562880e-3,   -2.6412985750e-3,   -6.3806306570e-3,   -5.4190107620e-3,
        1.7560585400e-3,    1.2098768730e-2,    1.9286893310e-2,    1.7849616710e-2,
        7.2509087620e-3,
    ];
    let expected_deviations = vec![
        5.2750685600e-8,    5.2750685600e-8,
    ];
    let expected_extremal_frequencies = vec![
        6.0975607480e-3,    1.2195121500e-2,    2.4390242990e-2,    3.0487803740e-2,
        4.2682923380e-2,    5.4878048600e-2,    6.0975611210e-2,    6.7073173820e-2,
        7.3170736430e-2,    7.9268299040e-2,    8.5365861650e-2,    9.1463424270e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0609755810e-1,    2.1219511330e-1,
        2.1829266850e-1,    2.2439022360e-1,    2.3048777880e-1,    2.4268288910e-1,
        2.4878044430e-1,    2.6097556950e-1,    2.7317067980e-1,    2.8536579010e-1,
        2.9756090040e-1,    3.0975601080e-1,    3.2195112110e-1,    3.3414623140e-1,
        3.4024378660e-1,    3.5243889690e-1,    3.6463400720e-1,    3.7682911750e-1,
        3.8902422790e-1,    4.0121933820e-1,    4.1341444850e-1,    4.2560955880e-1,
        4.3780466910e-1,    4.4999977950e-1,    4.6219488980e-1,    4.7439000010e-1,
        4.8658511040e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
82,2,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n82_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        4.9786962110e-8,   -1.5759169970e-7,   -4.5078144240e-7,   -1.4910955310e-7,
        1.3401318030e-6,    2.8493163880e-6,    1.1982263000e-6,   -5.3253074840e-6,
        -1.1815620380e-5,   -7.3324072220e-6,    1.3646693330e-5,    3.6857101800e-5,
        3.0829014580e-5,   -2.2631269530e-5,   -9.1590525700e-5,   -9.7961106800e-5,
        1.3161203240e-5,    1.8611502310e-4,    2.5139501670e-4,    6.3306040830e-5,
        -3.0935564430e-4,   -5.4268195530e-4,   -2.9686454220e-4,    4.0217407510e-4,
        1.0101955850e-3,    8.2286115500e-4,   -3.2528617880e-4,   -1.6434079730e-3,
        -1.8188937100e-3,   -1.8249778080e-4,    2.3358666800e-3,    3.5242915620e-3,
        1.6163394320e-3,   -2.7905695610e-3,   -6.4050490040e-3,   -5.2855378020e-3,
        1.9334577960e-3,    1.2157833200e-2,    1.9168151540e-2,    1.7656710000e-2,
        7.1584731340e-3,
    ];
    let expected_deviations = vec![
        8.8689894540e-8,    8.8689894540e-8,
    ];
    let expected_extremal_frequencies = vec![
        4.0650404990e-3,    1.2195121500e-2,    2.4390242990e-2,    3.2520323990e-2,
        4.0650404990e-2,    5.2845526490e-2,    6.0975607480e-2,    7.3170721530e-2,
        7.7235758300e-2,    8.5365831850e-2,    8.9430868630e-2,    9.3495905400e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0406503980e-1,    2.0813007650e-1,
        2.1219511330e-1,    2.2032518680e-1,    2.2845526040e-1,    2.4065037070e-1,
        2.4878044430e-1,    2.6097556950e-1,    2.7317067980e-1,    2.8536579010e-1,
        2.9349586370e-1,    3.0569097400e-1,    3.1788608430e-1,    3.3008119460e-1,
        3.4227630500e-1,    3.5447141530e-1,    3.6666652560e-1,    3.7886163590e-1,
        3.9105674620e-1,    4.0325185660e-1,    4.1544696690e-1,    4.2764207720e-1,
        4.3983718750e-1,    4.5203229780e-1,    4.6422740820e-1,    4.7642251850e-1,
        4.8861762880e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
82,2,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n82_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        2.8446385160e-7,    6.3681682150e-8,   -1.0551857490e-6,   -1.5357045410e-6,
        5.8918476500e-7,    5.0190124060e-6,    6.2020317270e-6,   -2.0864854380e-6,
        -1.6761507140e-5,   -2.0914079870e-5,    2.3078282540e-6,    4.3837877460e-5,
        6.0191094240e-5,    8.4189669000e-6,   -9.3813192510e-5,   -1.4919941900e-4,
        -5.5520169550e-5,    1.6597793730e-4,    3.2260423180e-4,    1.8927251220e-4,
        -2.3669176150e-4,   -6.1560241740e-4,   -4.9112015400e-4,    2.3967371090e-4,
        1.0456816530e-3,    1.0746095800e-3,   -4.2426749130e-5,   -1.5866318720e-3,
        -2.0869453440e-3,   -5.9088249690e-4,    2.1332472100e-3,    3.7406173070e-3,
        2.1163113420e-3,   -2.4135313000e-3,   -6.4934631810e-3,   -5.8046486230e-3,
        1.3978183270e-3,    1.2062057850e-2,    1.9614160060e-2,    1.8287193030e-2,
        7.4500776830e-3,
    ];
    let expected_deviations = vec![
        5.4456737790e-8,    5.4456737790e-8,
    ];
    let expected_extremal_frequencies = vec![
        4.8780487850e-3,    8.9430902150e-3,    1.3008131650e-2,    1.7073171210e-2,
        2.2764222700e-2,    2.8455274180e-2,    3.5772342230e-2,    4.9593467270e-2,
        6.3414596020e-2,    7.6422713700e-2,    8.6178801950e-2,    9.5934890210e-2,
        9.9186919630e-2,    1.0000000150e-1,    2.0000000300e-1,    2.0162601770e-1,
        2.0650406180e-1,    2.1382112800e-1,    2.2276420890e-1,    2.3333330450e-1,
        2.4390240010e-1,    2.5528451800e-1,    2.6666662100e-1,    2.7886173130e-1,
        2.9024383430e-1,    3.0162593720e-1,    3.1382104750e-1,    3.2601615790e-1,
        3.3739826080e-1,    3.4959337120e-1,    3.6260148880e-1,    3.7479659910e-1,
        3.8617870210e-1,    3.9918681980e-1,    4.1138193010e-1,    4.2357704040e-1,
        4.3658515810e-1,    4.4878026840e-1,    4.6178838610e-1,    4.7398349640e-1,
        4.8699161410e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
82,2,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n82_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -1.7106265200e-7,   -4.8710069220e-7,   -3.1194576880e-7,    1.0895910240e-6,
        3.2459524850e-6,    3.2128079970e-6,   -2.3225343280e-6,   -1.1791684300e-5,
        -1.5301653550e-5,   -1.0048061090e-6,    2.9165659730e-5,    4.9137721360e-5,
        2.3642074670e-5,   -5.1912978960e-5,   -1.2163483190e-4,   -9.6326250060e-5,
        5.8368692410e-5,    2.4443035360e-4,    2.6694740520e-4,    5.9188460000e-6,
        -4.0370645000e-4,   -5.8965344210e-4,   -2.3779459300e-4,    5.3355429550e-4,
        1.1005960400e-3,    7.7676429650e-4,   -4.8666330990e-4,   -1.7828084530e-3,
        -1.8002201100e-3,   -4.5613851400e-6,    2.5211586620e-3,    3.5425093960e-3,
        1.4372672890e-3,   -3.0113726390e-3,   -6.4621949570e-3,   -5.1178210410e-3,
        2.1767895670e-3,    1.2250510980e-2,    1.9019924100e-2,    1.7402796070e-2,
        7.0354230700e-3,
    ];
    let expected_deviations = vec![
        5.0181490960e-8,    5.0181490960e-8,
    ];
    let expected_extremal_frequencies = vec![
        7.6219509360e-4,    3.8109754680e-3,    1.2957318690e-2,    1.7530487850e-2,
        2.2103654220e-2,    2.5914626200e-2,    2.9725598170e-2,    8.0792628230e-2,
        8.8414572180e-2,    9.3749932940e-2,    9.6036516130e-2,    9.9085293710e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0152439180e-1,    2.0533536370e-1,
        2.1067072450e-1,    2.1829266850e-1,    2.2667680680e-1,    2.3582313950e-1,
        2.4573166670e-1,    2.5640240310e-1,    2.6707312460e-1,    2.7850604060e-1,
        2.8917676210e-1,    3.0060967800e-1,    3.1280478840e-1,    3.2423770430e-1,
        3.3567062020e-1,    3.4786573050e-1,    3.5929864650e-1,    3.7149375680e-1,
        3.8368886710e-1,    3.9588397740e-1,    4.0807908770e-1,    4.2027419810e-1,
        4.3246930840e-1,    4.4542661310e-1,    4.5838391780e-1,    4.7134122250e-1,
        4.8582291600e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,2,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n128_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -3.1513954470e-7,   -3.2277776540e-7,    3.6025380720e-7,    2.6645318480e-6,
        4.2684150690e-6,   -1.2785599210e-6,   -1.6559286450e-5,   -2.6655941840e-5,
        -1.1845149860e-6,    6.8386863860e-5,    1.2025843900e-4,    3.9469232430e-5,
        -2.0915143250e-4,   -4.2786917770e-4,   -2.4571767430e-4,    4.8571231310e-4,
        1.2495913540e-3,    1.0009646650e-3,   -8.1020896320e-4,   -3.0613308770e-3,
        -3.1535341400e-3,    6.5071275460e-4,    6.3485577700e-3,    8.1955585630e-3,
        1.4670258390e-3,   -1.1087709110e-2,   -1.8146442250e-2,   -8.4756575520e-3,
        1.5810746700e-2,    3.4818504010e-2,    2.4879552420e-2,   -1.6485434030e-2,
        -5.8286353950e-2,   -5.5896840990e-2,    5.8574602010e-3,    8.4837101400e-2,
        1.0525161770e-1,    2.5830678640e-2,   -1.0532234610e-1,   -1.7189446090e-1,
        -8.8051348920e-2,    1.0528305920e-1,    2.4687424300e-1,    1.8533326690e-1,
        -6.7990407350e-2,   -3.1228679420e-1,   -3.1224533920e-1,   -1.9599467520e-2,
        3.4403407570e-1,    4.5056310300e-1,    1.5971291070e-1,   -3.1873953340e-1,
        -5.7106882330e-1,   -3.3860400320e-1,    2.2300222520e-1,    6.4087498190e-1,
        5.2663415670e-1,   -6.1407506470e-2,   -6.3487035040e-1,   -6.8547320370e-1,
        -1.4000499250e-1,    5.5197918420e-1,    7.9534369710e-1,    3.6967259650e-1,
    ];
    let expected_deviations = vec![
        2.7755863390e-8,    2.7755863390e-8,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    7.8125000000e-3,    1.1718750000e-2,    1.5625000000e-2,
        1.9531250000e-2,    2.7343750000e-2,    3.1250000000e-2,    4.2968750000e-2,
        5.0781250000e-2,    5.8593750000e-2,    6.2500000000e-2,    6.6406250000e-2,
        7.0312500000e-2,    7.4218750000e-2,    7.8125000000e-2,    8.2031250000e-2,
        8.5937500000e-2,    8.9843750000e-2,    9.3750000000e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0390625300e-1,    2.0781250300e-1,    2.1171875300e-1,
        2.1562500300e-1,    2.1953125300e-1,    2.2343750300e-1,    2.2734375300e-1,
        2.3515625300e-1,    2.3906250300e-1,    2.4687500300e-1,    2.5468748810e-1,
        2.5859373810e-1,    2.6640623810e-1,    2.7421873810e-1,    2.8203123810e-1,
        2.8984373810e-1,    2.9374998810e-1,    3.0156248810e-1,    3.0937498810e-1,
        3.1718748810e-1,    3.2499998810e-1,    3.3281248810e-1,    3.4062498810e-1,
        3.4843748810e-1,    3.5624998810e-1,    3.6406248810e-1,    3.7187498810e-1,
        3.7968748810e-1,    3.8359373810e-1,    3.9140623810e-1,    3.9921873810e-1,
        4.0703123810e-1,    4.1484373810e-1,    4.2265623810e-1,    4.3046873810e-1,
        4.3828123810e-1,    4.4609373810e-1,    4.5390623810e-1,    4.6171873810e-1,
        4.6953123810e-1,    4.7734373810e-1,    4.8515623810e-1,    4.9296873810e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,2,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n128_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -2.0711274830e-7,   -1.9783118430e-7,    2.5412657580e-7,    1.5619959870e-6,
        2.2870858630e-6,   -8.6511295190e-7,   -8.6617237680e-6,   -1.3162325560e-5,
        -1.3791941460e-7,    3.2674761310e-5,    5.5315376810e-5,    1.7296901210e-5,
        -9.2387846960e-5,   -1.8452102090e-4,   -1.0453093270e-4,    1.9956800680e-4,
        5.0799298330e-4,    4.0496859580e-4,   -3.0867569150e-4,   -1.1788127010e-3,
        -1.2136772270e-3,    2.1333724730e-4,    2.3249795190e-3,    3.0112515670e-3,
        6.0033146290e-4,   -3.8727449720e-3,   -6.3919946550e-3,   -3.0779717490e-3,
        5.2666757260e-3,    1.1806070800e-2,    8.5719358180e-3,   -5.1783695820e-3,
        -1.9096322360e-2,   -1.8522985280e-2,    1.3978965580e-3,    2.6945708320e-2,
        3.3786106850e-2,    8.8571421800e-3,   -3.2511189580e-2,   -5.3730193530e-2,
        -2.8079196810e-2,    3.1599450860e-2,    7.5477361680e-2,    5.7177681480e-2,
        -1.9630230960e-2,   -9.3784660100e-2,   -9.4207450750e-2,   -6.6953301430e-3,
        1.0196445880e-1,    1.3378739360e-1,    4.7539666300e-2,   -9.3850106000e-2,
        -1.6779975590e-1,   -9.8721027370e-2,    6.6284567120e-2,    1.8756183980e-1,
        1.5197940170e-1,   -2.1274492140e-2,   -1.8733105060e-1,   -1.9745673240e-1,
        -3.2492443920e-2,    1.7286400500e-1,    2.4241515990e-1,    1.1196874080e-1,
    ];
    let expected_deviations = vec![
        2.7402039750e-8,    2.7402039750e-8,
    ];
    let expected_extremal_frequencies = vec![
        5.2083334890e-3,    1.0416666980e-2,    1.3020833950e-2,    1.5625000000e-2,
        2.3437498140e-2,    2.8645830230e-2,    3.3854164180e-2,    3.9062500000e-2,
        4.6875003730e-2,    5.2083339540e-2,    5.9895843270e-2,    6.7708335820e-2,
        7.2916664180e-2,    7.5520828370e-2,    8.0729156730e-2,    8.3333320920e-2,
        8.8541649280e-2,    9.3749977650e-2,    9.6354141830e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0260417460e-1,    2.0520834620e-1,    2.0781251790e-1,
        2.1041668950e-1,    2.1562503280e-1,    2.2083337600e-1,    2.2604171930e-1,
        2.3125006260e-1,    2.3906257750e-1,    2.4427092080e-1,    2.5208342080e-1,
        2.5989589100e-1,    2.6510420440e-1,    2.7291667460e-1,    2.8072914480e-1,
        2.8854161500e-1,    2.9374992850e-1,    3.0156239870e-1,    3.0937486890e-1,
        3.1718733910e-1,    3.2499980930e-1,    3.3281227950e-1,    3.4062474970e-1,
        3.4583306310e-1,    3.5364553330e-1,    3.6145800350e-1,    3.6927047370e-1,
        3.7708294390e-1,    3.8489541410e-1,    3.9270788430e-1,    4.0052035450e-1,
        4.0833282470e-1,    4.1614529490e-1,    4.2395776510e-1,    4.3177023530e-1,
        4.3958270550e-1,    4.4739517570e-1,    4.5260348920e-1,    4.6041595940e-1,
        4.6822842960e-1,    4.7604089980e-1,    4.8385337000e-1,    4.9166584010e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
128,2,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n128_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -2.1510792070e-5,    1.3608901100e-5,    1.9305622120e-5,    1.0594510970e-5,
        -1.0471035240e-5,   -2.5335986720e-5,   -7.9373512560e-6,    4.2479325200e-5,
        7.2044400440e-5,    1.0329313230e-6,   -1.7256593860e-4,   -2.8468383240e-4,
        -8.8461005360e-5,    4.5328194390e-4,    9.1296335450e-4,    5.8558932510e-4,
        -7.7871064420e-4,   -2.2762750740e-3,   -2.1483832970e-3,    6.0775014570e-4,
        4.4934181500e-3,    5.7245921340e-3,    1.3276953250e-3,   -6.9553153590e-3,
        -1.2114831250e-2,   -7.0347329600e-3,    7.7292118220e-3,    2.1063178780e-2,
        1.8623154610e-2,   -3.3738687630e-3,   -3.0166577550e-2,   -3.6873996260e-2,
        -1.0245978830e-2,    3.4336037930e-2,    5.9426847850e-2,    3.5892192270e-2,
        -2.6700187470e-2,   -7.9645290970e-2,   -7.2186492380e-2,    1.3662800190e-3,
        8.7443560360e-2,    1.1159774660e-1,    4.2655795810e-2,   -7.2716243570e-2,
        -1.4094690980e-1,   -9.8266609010e-2,    3.0552372340e-2,    1.4537110920e-1,
        1.5014937520e-1,    3.4142494200e-2,   -1.1484231050e-1,   -1.7902186510e-1,
        -1.0555842520e-1,    5.0396740440e-2,    1.6936515270e-1,    1.6124086080e-1,
        3.3520370720e-2,   -1.1749073120e-1,   -1.8121254440e-1,   -1.1326278750e-1,
        3.6496788260e-2,    1.6317869720e-1,    1.8048115070e-1,    7.7295675870e-2,
    ];
    let expected_deviations = vec![
        5.3452023250e-8,    5.3452023250e-8,
    ];
    let expected_extremal_frequencies = vec![
        1.5625001400e-3,    6.2500000930e-3,    7.2916666980e-3,    1.3020837680e-2,
        1.4062505220e-2,    1.7187504100e-2,    2.6041662320e-2,    3.4374989570e-2,
        4.3229147790e-2,    4.8437476160e-2,    5.5208303030e-2,    6.0416631400e-2,
        6.6145792600e-2,    7.2395786640e-2,    7.9166613520e-2,    8.5937440400e-2,
        9.1666601600e-2,    9.6354097130e-2,    9.8958261310e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0156252380e-1,    2.0468756560e-1,    2.0833344760e-1,
        2.1302101020e-1,    2.1875025330e-1,    2.2447949650e-1,    2.3020873960e-1,
        2.3593798280e-1,    2.4166722600e-1,    2.4843814970e-1,    2.5520890950e-1,
        2.6197963950e-1,    2.6927119490e-1,    2.7604192500e-1,    2.8333348040e-1,
        2.9010421040e-1,    2.9739576580e-1,    3.0468732120e-1,    3.1197887660e-1,
        3.1927043200e-1,    3.2656198740e-1,    3.3437436820e-1,    3.4166592360e-1,
        3.4895747900e-1,    3.5624903440e-1,    3.6406141520e-1,    3.7135297060e-1,
        3.7916535140e-1,    3.8645690680e-1,    3.9374846220e-1,    4.0156084300e-1,
        4.0885239840e-1,    4.1666477920e-1,    4.2395633460e-1,    4.3176871540e-1,
        4.3906027080e-1,    4.4687265160e-1,    4.5468503240e-1,    4.6197658780e-1,
        4.6978896860e-1,    4.7708052400e-1,    4.8489290480e-1,    4.9218446020e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
128,2,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn differentiator_lowpass_n128_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        3.2223226930e2,    -2.8080657960e2,    -3.0254644780e2,    -5.7499633790e1,
        2.1317926030e2,     2.6564596560e2,     8.5291351320e1,    -1.2632241820e2,
        -1.7555563350e2,    -6.3599296570e1,     5.2080970760e1,     5.4447658540e1,
        -1.0618235590e1,    -1.7072372440e1,     6.5621246340e1,     1.2217527010e2,
        3.3844024660e1,    -1.5282189940e2,    -2.4203286740e2,    -9.6426559450e1,
        1.8581002810e2,     3.3629479980e2,     1.8180535890e2,    -1.6087547300e2,
        -3.7710150150e2,    -2.5703771970e2,     9.3580810550e1,     3.5198883060e2,
        2.8994607540e2,    -1.4110168460e1,    -2.6857342530e2,    -2.6004138180e2,
        -4.2460769650e1,     1.5273097230e2,     1.6599667360e2,     4.7968315120e1,
        -4.0478439330e1,    -2.6847213750e1,     8.9635772710e0,    -3.3994552610e1,
        -1.2394583130e2,    -1.1867782590e2,     5.0041412350e1,     2.4986248780e2,
        2.5382322690e2,    -7.5940246580e0,    -3.2413220210e2,    -3.7944296260e2,
        -7.3482421880e1,     3.3902163700e2,     4.6613650510e2,     1.6181643680e2,
        -3.0766625980e2,    -5.0100216670e2,    -2.2752615360e2,     2.5741430660e2,
        4.9163793950e2,     2.5434211730e2,    -2.1733831790e2,    -4.6128579710e2,
        -2.4549951170e2,     2.0535134890e2,     4.3716537480e2,     2.2070025630e2,
    ];
    let expected_deviations = vec![
        6.3752246150e-8,    6.3752246150e-8,
    ];
    let expected_extremal_frequencies = vec![
        2.9296875000e-3,    3.4179687500e-3,    5.3710937500e-3,    6.3476562500e-3,
        1.1230468750e-2,    1.9531250000e-2,    2.3925781250e-2,    2.8808593750e-2,
        3.5156250000e-2,    4.1503906250e-2,    4.8339843750e-2,    5.4687500000e-2,
        6.0546875000e-2,    6.6894531250e-2,    7.3242187500e-2,    8.0078125000e-2,
        8.6425781250e-2,    9.2285156250e-2,    9.6679687500e-2,    9.8632812500e-2,
        1.0000000150e-1,    2.0097656550e-1,    2.0341797170e-1,    2.0830078420e-1,
        2.1464844050e-1,    2.2099609670e-1,    2.2832031550e-1,    2.3515625300e-1,
        2.4150390920e-1,    2.4833984670e-1,    2.5517576930e-1,    2.6201170680e-1,
        2.6835936310e-1,    2.7519530060e-1,    2.8203123810e-1,    2.8886717560e-1,
        2.9570311310e-1,    3.0253905060e-1,    3.0937498810e-1,    3.1621092560e-1,
        3.2304686310e-1,    3.2988280060e-1,    3.3720701930e-1,    3.4404295680e-1,
        3.5087889430e-1,    3.5820311310e-1,    3.6503905060e-1,    3.7187498810e-1,
        3.7968748810e-1,    3.8652342560e-1,    3.9384764430e-1,    4.0117186310e-1,
        4.0849608180e-1,    4.1582030060e-1,    4.2314451930e-1,    4.3046873810e-1,
        4.3779295680e-1,    4.4511717560e-1,    4.5244139430e-1,    4.5976561310e-1,
        4.6708983180e-1,    4.7441405060e-1,    4.8222655060e-1,    4.9003905060e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n3_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        1.1342786250e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.9879518450e-1,    3.9759036900e-1,    1.9879518450e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n3_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        1.1342786250e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.9879518450e-1,    3.9759036900e-1,    1.9879518450e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n3_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        1.1342786250e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.9879518450e-1,    3.9759036900e-1,    1.9879518450e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n3_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        1.1342786250e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.9879518450e-1,    3.9759036900e-1,    1.9879518450e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n3_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        1.1342786250e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.9879518450e-1,    3.9759036900e-1,    1.9879518450e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,2,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n3_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        1.1342786250e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.9879518450e-1,    3.9759036900e-1,    1.9879518450e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n4_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        2.9650272800e-2,    1.2075543400e-1,
    ];
    let expected_deviations = vec![
        1.8221031130e-1,    3.6442062260e-1,    1.8221031130e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n4_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        3.4772951160e-2,    1.1779715120e-1,
    ];
    let expected_deviations = vec![
        1.8943883480e-1,    3.7887766960e-1,    1.8943883480e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    3.7999999520e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n4_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        3.4772951160e-2,    1.1779715120e-1,
    ];
    let expected_deviations = vec![
        1.8943883480e-1,    3.7887766960e-1,    1.8943883480e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    3.7999999520e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n4_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        3.4772951160e-2,    1.1779715120e-1,
    ];
    let expected_deviations = vec![
        1.8943883480e-1,    3.7887766960e-1,    1.8943883480e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    3.7999999520e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n4_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        3.4772951160e-2,    1.1779715120e-1,
    ];
    let expected_deviations = vec![
        1.8943883480e-1,    3.7887766960e-1,    1.8943883480e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    3.7999999520e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,2,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n4_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        3.4772951160e-2,    1.1779715120e-1,
    ];
    let expected_deviations = vec![
        1.8943883480e-1,    3.7887766960e-1,    1.8943883480e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    3.7999999520e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n11_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -2.8055142610e-3,    4.6177979560e-2,   -4.1350275280e-2,   -9.5180831850e-3,
        9.6230573950e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.8748285770e-2,    1.5749657150e-1,    7.8748285770e-2,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    2.1999999880e-1,    2.7000001070e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.3000000720e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n11_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        1.3242266140e-4,    4.6806637200e-2,   -3.9494015280e-2,   -1.0722141710e-2,
        9.6503615380e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.0774642530e-2,    1.6154928510e-1,    8.0774642530e-2,
    ];
    let expected_extremal_frequencies = vec![
        6.6666670140e-2,    2.1999999880e-1,    2.8666666150e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.4666665790e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n11_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -7.5297494190e-4,    4.6298243110e-2,   -3.8997214290e-2,   -9.9977776410e-3,
        9.6673347060e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.1495076420e-2,    1.6299015280e-1,    8.1495076420e-2,
    ];
    let expected_extremal_frequencies = vec![
        7.3333337900e-2,    2.1999999880e-1,    2.7999994160e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.4666659830e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n11_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -8.5504894380e-4,    4.6251747760e-2,   -3.8965813820e-2,   -9.9183507260e-3,
        9.6689492460e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.1538669760e-2,    1.6307733950e-1,    8.1538669760e-2,
    ];
    let expected_extremal_frequencies = vec![
        7.5000010430e-2,    2.1999999880e-1,    2.7624994520e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.4874992970e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n11_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        -8.2321849190e-4,    4.6256184580e-2,   -3.8961861280e-2,   -9.9449791010e-3,
        9.6683964130e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.1551820040e-2,    1.6310364010e-1,    8.1551820040e-2,
    ];
    let expected_extremal_frequencies = vec![
        7.5000002980e-2,    2.1999999880e-1,    2.7833348510e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.4666683670e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,2,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n11_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        -8.3988648840e-4,    4.6247135850e-2,   -3.8953162730e-2,   -9.9312104280e-3,
        9.6687205140e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.1564046440e-2,    1.6312809290e-1,    8.1564046440e-2,
    ];
    let expected_extremal_frequencies = vec![
        7.5675673780e-2,    2.1999999880e-1,    2.7675679330e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.4756782050e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n12_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -1.6991384330e-2,    1.1662342590e-2,   -5.8356644590e-3,   -8.4946379070e-2,
        4.6412225810e-2,    6.0056108980e-2,
    ];
    let expected_deviations = vec![
        7.3626302180e-2,    1.4725260440e-1,    7.3626302180e-2,
    ];
    let expected_extremal_frequencies = vec![
        4.1666667910e-2,    1.7000000180e-1,    2.1999999880e-1,    2.6166665550e-1,
        3.3000001310e-1,    3.7999999520e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n12_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -1.6582900660e-2,    1.0167076250e-2,   -4.7325314950e-3,   -8.2258455460e-2,
        4.6491362150e-2,    6.0617581010e-2,
    ];
    let expected_deviations = vec![
        7.8209385280e-2,    1.5641877060e-1,    7.8209385280e-2,
    ];
    let expected_extremal_frequencies = vec![
        5.5555555970e-2,    1.7000000180e-1,    2.1999999880e-1,    2.7555555110e-1,
        3.3000001310e-1,    3.7999999520e-1,    4.6333336830e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n12_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -1.6558373350e-2,    1.0061170910e-2,   -4.6699522060e-3,   -8.2170121370e-2,
        4.6515431260e-2,    6.0597851870e-2,
    ];
    let expected_deviations = vec![
        7.8438691790e-2,    1.5687738360e-1,    7.8438691790e-2,
    ];
    let expected_extremal_frequencies = vec![
        6.1111111190e-2,    1.7000000180e-1,    2.1999999880e-1,    2.8111118080e-1,
        3.3000001310e-1,    3.7999999520e-1,    4.5777797700e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n12_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -1.6565104950e-2,    1.0095829140e-2,   -4.6864929610e-3,   -8.2182742660e-2,
        4.6506188810e-2,    6.0612924400e-2,
    ];
    let expected_deviations = vec![
        7.8379035000e-2,    1.5675807000e-1,    7.8379035000e-2,
    ];
    let expected_extremal_frequencies = vec![
        5.7291660460e-2,    1.7000000180e-1,    2.1999999880e-1,    2.7729168530e-1,
        3.3000001310e-1,    3.7999999520e-1,    4.5812514420e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n12_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        -1.6555609180e-2,    1.0052833710e-2,   -4.6626250260e-3,   -8.2155376670e-2,
        4.6517245470e-2,    6.0600101950e-2,
    ];
    let expected_deviations = vec![
        7.8464157880e-2,    1.5692831580e-1,    7.8464157880e-2,
    ];
    let expected_extremal_frequencies = vec![
        6.0185171660e-2,    1.7000000180e-1,    2.1999999880e-1,    2.8018513320e-1,
        3.3000001310e-1,    3.7999999520e-1,    4.5870339870e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,2,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n12_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        -1.6556860880e-2,    1.0060496630e-2,   -4.6656797640e-3,   -8.2157731060e-2,
        4.6515706930e-2,    6.0603365300e-2,
    ];
    let expected_deviations = vec![
        7.8451439740e-2,    1.5690287950e-1,    7.8451439740e-2,
    ];
    let expected_extremal_frequencies = vec![
        5.8558546010e-2,    1.7000000180e-1,    2.1999999880e-1,    2.8081077340e-1,
        3.3000001310e-1,    3.7999999520e-1,    4.5882877710e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n16_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        2.6949010790e-3,    2.1713169290e-2,   -8.4248632190e-3,    3.5744622350e-2,
        -2.6815123860e-3,   -6.7240893840e-2,    3.8738157600e-2,    6.8313561380e-2,
    ];
    let expected_deviations = vec![
        5.6407529860e-2,    1.1281505970e-1,    5.6407529860e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.1250000000e-2,    1.2500000000e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.8249999880e-1,    3.3000001310e-1,    3.7999999520e-1,    4.1124999520e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n16_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        3.8004643280e-3,    2.1169977260e-2,   -7.8948438170e-3,    3.5158112650e-2,
        -3.3656992020e-3,   -6.6898569460e-2,    3.7345543500e-2,    6.9589875640e-2,
    ];
    let expected_deviations = vec![
        5.8267913760e-2,    1.1653582750e-1,    5.8267913760e-2,
    ];
    let expected_extremal_frequencies = vec![
        4.1666667910e-2,    1.2500000000e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.8249999880e-1,    3.3000001310e-1,    3.7999999520e-1,    4.2166668180e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n16_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        3.4357719120e-3,    2.0919337870e-2,   -8.1226062030e-3,    3.5105980930e-2,
        -2.9719397430e-3,   -6.6223859790e-2,    3.7813570350e-2,    6.9750383500e-2,
    ];
    let expected_deviations = vec![
        5.8793969450e-2,    1.1758793890e-1,    5.8793969450e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.7500001490e-2,    1.1666671190e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.8666660190e-1,    3.3000001310e-1,    3.7999999520e-1,    4.2583328490e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n16_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        3.4429607910e-3,    2.0923987030e-2,   -8.1201251600e-3,    3.5093143580e-2,
        -2.9843077060e-3,   -6.6229365770e-2,    3.7794880570e-2,    6.9750070570e-2,
    ];
    let expected_deviations = vec![
        5.8808684350e-2,    1.1761736870e-1,    5.8808684350e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-2,    1.1718750000e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.9031249880e-1,    3.3000001310e-1,    3.7999999520e-1,    4.2296874520e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n16_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        3.4399917350e-3,    2.0924989130e-2,   -8.1195421520e-3,    3.5086318850e-2,
        -2.9763355850e-3,   -6.6207624970e-2,    3.7796504800e-2,    6.9755211470e-2,
    ];
    let expected_deviations = vec![
        5.8836426590e-2,    1.1767285320e-1,    5.8836426590e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.8194455210e-2,    1.1631949250e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.8944417830e-1,    3.3000001310e-1,    3.7999999520e-1,    4.2340260740e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,2,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n16_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        3.4442855980e-3,    2.0925818010e-2,   -8.1168450420e-3,    3.5083442930e-2,
        -2.9792897400e-3,   -6.6207215190e-2,    3.7790112200e-2,    6.9757811730e-2,
    ];
    let expected_deviations = vec![
        5.8843079950e-2,    1.1768615990e-1,    5.8843079950e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.8851354270e-2,    1.1655401440e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.8925701980e-1,    3.3000001310e-1,    3.7999999520e-1,    4.2391908170e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n47_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -1.1725790100e-3,   -1.1696979640e-3,    2.0652972160e-3,   -1.6750779470e-4,
        8.8184897320e-5,    3.1074183060e-4,   -3.9497818800e-3,    2.4518612770e-3,
        4.4260043650e-3,   -4.3017938730e-3,    9.0193713550e-4,   -1.9723312000e-3,
        -4.4194180520e-3,    1.3148019090e-2,   -1.5462331940e-3,   -1.2854717670e-2,
        7.7716703530e-3,   -7.8701470050e-3,    6.7526618950e-3,    3.1828116630e-2,
        -4.5875854790e-2,   -3.0556550250e-2,    8.1854797900e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.8199387700e-3,    3.6398775410e-3,    1.8199387700e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0869565420e-2,    3.2608695330e-2,    5.4347828030e-2,    7.6086953280e-2,
        9.7826078530e-2,    1.1956520380e-1,    1.4130432900e-1,    1.5217389170e-1,
        1.7000000180e-1,    2.3086956140e-1,    2.4173912410e-1,    2.5260868670e-1,
        2.7434781190e-1,    2.9608693720e-1,    3.0695649980e-1,    3.1782606240e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.9086955790e-1,    4.0173912050e-1,
        4.2347824570e-1,    4.4521737100e-1,    4.6695649620e-1,    4.8869562150e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n47_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -1.2006278850e-3,   -1.4402503150e-3,    2.0281341860e-3,   -3.0295050240e-4,
        4.8276735470e-5,    4.6927866060e-4,   -3.9294622840e-3,    2.3031996100e-3,
        4.4274576940e-3,   -4.2733484880e-3,    9.4541884030e-4,   -1.8095877020e-3,
        -4.5814942570e-3,    1.2942502280e-2,   -1.3764188620e-3,   -1.2823479250e-2,
        7.8282365580e-3,   -7.7760899440e-3,    6.4764432610e-3,    3.1825438140e-2,
        -4.5700456950e-2,   -3.0648406590e-2,    8.1943154340e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.3112313360e-3,    4.6224626710e-3,    2.3112313360e-3,
    ];
    let expected_extremal_frequencies = vec![
        7.2463769470e-3,    2.8985507790e-2,    5.0724633040e-2,    7.2463758290e-2,
        9.4202883540e-2,    1.1594200880e-1,    1.3768114150e-1,    1.5942026670e-1,
        2.1999999880e-1,    2.2724637390e-1,    2.4173912410e-1,    2.5623187420e-1,
        2.7072462440e-1,    2.9246374960e-1,    3.0695649980e-1,    3.2144925000e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8724637030e-1,    4.0173912050e-1,
        4.2347824570e-1,    4.4521737100e-1,    4.6695649620e-1,    4.8869562150e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n47_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -1.1054897910e-3,   -1.5009421620e-3,    1.9959085620e-3,   -2.9849458950e-4,
        -1.0565272530e-5,    6.3062505800e-4,   -3.8046641280e-3,    2.1502883170e-3,
        4.3793660590e-3,   -4.3626651170e-3,    9.4925670420e-4,   -1.5242705120e-3,
        -4.5917504470e-3,    1.2734182180e-2,   -1.4052095360e-3,   -1.2904365550e-2,
        8.0039240420e-3,   -7.5320182370e-3,    6.2516028990e-3,    3.1691081820e-2,
        -4.5716267080e-2,   -3.0659049750e-2,    8.2275539640e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.5166969280e-3,    5.0333938560e-3,    2.5166969280e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0144928470e-2,    3.1884063040e-2,    5.3623199460e-2,    7.5362302360e-2,
        9.7101382910e-2,    1.1739119140e-1,    1.3768100740e-1,    1.5942008790e-1,
        2.1999999880e-1,    2.2579708700e-1,    2.3884053530e-1,    2.5623187420e-1,
        2.7507260440e-1,    2.9246404770e-1,    3.0985549090e-1,    3.2434836030e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8579714300e-1,    4.0173929930e-1,
        4.2202931640e-1,    4.4231933360e-1,    4.6550792460e-1,    4.8869651560e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n47_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -1.1104131120e-3,   -1.4985630990e-3,    1.9990752920e-3,   -2.9945618010e-4,
        -4.2990432120e-6,    6.2379019800e-4,   -3.8098893130e-3,    2.1585740610e-3,
        4.3792142530e-3,   -4.3583302760e-3,    9.5314532520e-4,   -1.5368682800e-3,
        -4.5882733540e-3,    1.2741536830e-2,   -1.4062961560e-3,   -1.2897429060e-2,
        7.9944320020e-3,   -7.5444635000e-3,    6.2627573500e-3,    3.1694855540e-2,
        -4.5716732740e-2,   -3.0654521660e-2,    8.2262225450e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.5185523550e-3,    5.0371047110e-3,    2.5185523550e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0869564490e-2,    3.2608691600e-2,    5.4347816850e-2,    7.4728250500e-2,
        9.6467375760e-2,    1.1684780570e-1,    1.3858699800e-1,    1.5896753970e-1,
        2.1999999880e-1,    2.2543480990e-1,    2.3902183770e-1,    2.5668489930e-1,
        2.7434784170e-1,    2.9336947200e-1,    3.0967372660e-1,    3.2461929320e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8679343460e-1,    4.0173900130e-1,
        4.2211931940e-1,    4.4249963760e-1,    4.6559733150e-1,    4.8869502540e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n47_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        -1.1098011160e-3,   -1.5008666090e-3,    1.9984454850e-3,   -2.9793567960e-4,
        -4.0077720770e-6,    6.2351475940e-4,   -3.8106432180e-3,    2.1556196730e-3,
        4.3814284730e-3,   -4.3564094230e-3,    9.5175695610e-4,   -1.5388457100e-3,
        -4.5909751210e-3,    1.2741960590e-2,   -1.4024425760e-3,   -1.2896972710e-2,
        7.9924501480e-3,   -7.5449869040e-3,    6.2616309150e-3,    3.1699080020e-2,
        -4.5713562520e-2,   -3.0659276990e-2,    8.2260191440e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.5241351690e-3,    5.0482703370e-3,    2.5241351690e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0869566350e-2,    3.2004814590e-2,    5.3743984550e-2,    7.4879214170e-2,
        9.6618250010e-2,    1.1775342380e-1,    1.3828489180e-1,    1.5942032640e-1,
        2.1999999880e-1,    2.2543482480e-1,    2.3871995510e-1,    2.5623202320e-1,
        2.7434766290e-1,    2.9306715730e-1,    3.1057894230e-1,    3.2446759940e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8664239650e-1,    4.0173876290e-1,
        4.2166596650e-1,    4.4280087950e-1,    4.6514350180e-1,    4.8869383340e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,2,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n47_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        -1.1104266160e-3,   -1.5011280780e-3,    1.9991216250e-3,   -2.9691844250e-4,
        -3.3383257690e-6,    6.2202569100e-4,   -3.8126725700e-3,    2.1558203730e-3,
        4.3833460660e-3,   -4.3541612100e-3,    9.5131178390e-4,   -1.5423577280e-3,
        -4.5921942220e-3,    1.2744236740e-2,   -1.3995831830e-3,   -1.2896578760e-2,
        7.9889101910e-3,   -7.5482497920e-3,    6.2633147460e-3,    3.1702794130e-2,
        -4.5712281020e-2,   -3.0661944300e-2,    8.2255862650e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.5244045540e-3,    5.0488091070e-3,    2.5244045540e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0575793680e-2,    3.2314911480e-2,    5.3466543560e-2,    7.5205720960e-2,
        9.6357353030e-2,    1.1750898510e-1,    1.3807290790e-1,    1.5922427180e-1,
        2.1999999880e-1,    2.2528783980e-1,    2.3880121110e-1,    2.5583994390e-1,
        2.7464163300e-1,    2.9285576940e-1,    3.1048235300e-1,    3.2399606700e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8646307590e-1,    4.0173944830e-1,
        4.2171624300e-1,    4.4286814330e-1,    4.6519514920e-1,    4.8810970780e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n48_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -9.7024225400e-5,   -1.5954010890e-3,    1.0002888740e-3,    8.8704976950e-4,
        -1.0648396560e-3,    1.1751780980e-3,   -2.0244559270e-3,   -1.5301260860e-3,
        5.6266449390e-3,   -1.5160316830e-3,   -3.0198264870e-3,    2.4126721550e-3,
        -6.0653812250e-3,    5.2569331600e-3,    1.0086135010e-2,   -1.3874810190e-2,
        -7.6041882860e-4,    4.2909840120e-3,   -9.8428437490e-3,    2.9454484580e-2,
        -1.2113302950e-3,   -6.3883103430e-2,    3.6421731110e-2,    6.7320436240e-2,
    ];
    let expected_deviations = vec![
        1.3008746320e-3,    2.6017492640e-3,    1.3008746320e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0416666980e-2,    3.1250000000e-2,    5.2083335820e-2,    7.2916664180e-2,
        9.3749992550e-2,    1.1458332090e-1,    1.3541665670e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.3041667040e-1,    2.4083334210e-1,    2.5124999880e-1,
        2.7208331230e-1,    2.9291662570e-1,    3.0333328250e-1,    3.1374993920e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.9041665200e-1,    4.0083330870e-1,
        4.2166662220e-1,    4.4249993560e-1,    4.6333324910e-1,    4.8416656260e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n48_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        4.8192861870e-5,   -1.6859207540e-3,    9.0922234810e-4,    7.7887979570e-4,
        -1.4247745280e-3,    1.1273158020e-3,   -1.7473084150e-3,   -1.4769039120e-3,
        5.7293265130e-3,   -1.5200241470e-3,   -3.3227244860e-3,    2.4404507130e-3,
        -5.8709657750e-3,    5.0485925750e-3,    1.0145047680e-2,   -1.3822433540e-2,
        -8.5280276830e-4,    4.6088127420e-3,   -9.9125616250e-3,    2.9077526180e-2,
        -1.1196937410e-3,   -6.3908986750e-2,    3.6597892640e-2,    6.7627213900e-2,
    ];
    let expected_deviations = vec![
        1.7686664360e-3,    3.5373328720e-3,    1.7686664360e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.3888888990e-2,    3.4722223880e-2,    5.5555555970e-2,    7.6388895510e-2,
        9.7222238780e-2,    1.2500002980e-1,    1.5972226860e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2694444660e-1,    2.4083334210e-1,    2.5472223760e-1,
        2.7555558090e-1,    2.8944447640e-1,    3.1027781960e-1,    3.1722226740e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8694444300e-1,    4.0083333850e-1,
        4.1472223400e-1,    4.3555557730e-1,    4.5638892050e-1,    4.7722226380e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n48_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        4.6175727040e-4,   -1.6160226660e-3,    7.9022371210e-4,    8.6944841310e-4,
        -1.8664593810e-3,    8.1852573200e-4,   -1.4973655340e-3,   -1.8426983150e-3,
        5.7345717210e-3,   -1.1269124220e-3,   -3.6358451470e-3,    2.7418776880e-3,
        -5.6038731710e-3,    4.4394815340e-3,    1.0369718070e-2,   -1.3761751350e-2,
        -1.2505766940e-3,    5.1937461830e-3,   -1.0014681150e-2,    2.8678398580e-2,
        -6.0292333360e-4,   -6.4266592260e-2,    3.6509200930e-2,    6.8169981240e-2,
    ];
    let expected_deviations = vec![
        2.1925566250e-3,    4.3851132500e-3,    2.1925566250e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.2500000190e-2,    3.6111112680e-2,    6.1111111190e-2,    8.6111173030e-2,
        1.1250013110e-1,    1.4027798180e-1,    1.6111136970e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2555556890e-1,    2.3944449420e-1,    2.5750002260e-1,
        2.7555543180e-1,    2.9361084100e-1,    3.1027737260e-1,    3.2416614890e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8555550580e-1,    4.0083315970e-1,
        4.1888856890e-1,    4.3833285570e-1,    4.5916602020e-1,    4.7861030700e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n48_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        4.6096683950e-4,   -1.6202525000e-3,    7.8457756900e-4,    8.6958974130e-4,
        -1.8633840370e-3,    8.2327326530e-4,   -1.4952830970e-3,   -1.8436489630e-3,
        5.7346895340e-3,   -1.1224445190e-3,   -3.6338530480e-3,    2.7380760290e-3,
        -5.6090601720e-3,    4.4376561420e-3,    1.0375713000e-2,   -1.3757103120e-2,
        -1.2517534200e-3,    5.1890588370e-3,   -1.0018331930e-2,    2.8680978340e-2,
        -5.9637799860e-4,   -6.4267925920e-2,    3.6503691230e-2,    6.8167641760e-2,
    ];
    let expected_deviations = vec![
        2.1934167020e-3,    4.3868334030e-3,    2.1934167020e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.1718749070e-2,    3.6458332090e-2,    6.1197891830e-2,    8.5937514900e-2,
        1.1197923120e-1,    1.4062502980e-1,    1.6145828370e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2520831230e-1,    2.3953117430e-1,    2.5645825270e-1,
        2.7598965170e-1,    2.9421895740e-1,    3.1114616990e-1,    3.2416710260e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8651046160e-1,    4.0083348750e-1,
        4.1906279330e-1,    4.3859419230e-1,    4.5812559130e-1,    4.7895908360e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n48_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 36);

    let expected_impulse_response = vec![
        4.6037687570e-4,   -1.6195403880e-3,    7.8693625980e-4,    8.6759979600e-4,
        -1.8644458610e-3,    8.2216807640e-4,   -1.4964063880e-3,   -1.8406726890e-3,
        5.7349121200e-3,   -1.1261373290e-3,   -3.6332209130e-3,    2.7383903510e-3,
        -5.6068450210e-3,    4.4402335770e-3,    1.0370913890e-2,   -1.3759500350e-2,
        -1.2486078780e-3,    5.1896059890e-3,   -1.0015799660e-2,    2.8680339460e-2,
        -6.0150958600e-4,   -6.4265489580e-2,    3.6507032810e-2,    6.8167470400e-2,
    ];
    let expected_deviations = vec![
        2.1977715660e-3,    4.3955431320e-3,    2.1977715660e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.2152776120e-2,    3.6458332090e-2,    6.0763951390e-2,    8.5648126900e-2,
        1.1226839570e-1,    1.4004607500e-1,    1.6145803030e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2578701380e-1,    2.3967584970e-1,    2.5703689460e-1,
        2.7555534240e-1,    2.9407379030e-1,    3.1085613370e-1,    3.2474496960e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8578701020e-1,    4.0025454760e-1,
        4.1877299550e-1,    4.3844884630e-1,    4.5870339870e-1,    4.7895795110e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,2,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n48_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::Differentiator, &bands, 37);

    let expected_impulse_response = vec![
        4.5984305320e-4,   -1.6206983710e-3,    7.8636209950e-4,    8.6714664940e-4,
        -1.8639718180e-3,    8.2325248510e-4,   -1.4963366560e-3,   -1.8400338010e-3,
        5.7351817380e-3,   -1.1262237090e-3,   -3.6329643340e-3,    2.7378103700e-3,
        -5.6076995100e-3,    4.4404864310e-3,    1.0371175590e-2,   -1.3758959250e-2,
        -1.2483657340e-3,    5.1885913130e-3,   -1.0015827600e-2,    2.8680846100e-2,
        -6.0147792100e-4,   -6.4265191560e-2,    3.6506697540e-2,    6.8166524170e-2,
    ];
    let expected_deviations = vec![
        2.1981701720e-3,    4.3963403440e-3,    2.1981701720e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.1824322860e-2,    3.6036022010e-2,    6.0810782020e-2,    8.5585542020e-2,
        1.1261255290e-1,    1.4020282030e-1,    1.6159948710e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2563070060e-1,    2.3970745500e-1,    2.5659939650e-1,
        2.7518022060e-1,    2.9376104470e-1,    3.1121575830e-1,    3.2472908500e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8619360330e-1,    4.0026998520e-1,
        4.1885080930e-1,    4.3799468870e-1,    4.5882773400e-1,    4.7909772400e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,2,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n81_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        -2.5954053850e-5,   -6.3178456910e-5,    5.9250847700e-5,   -2.3018037610e-5,
        6.1138212910e-5,    1.3605927230e-4,   -3.3796278880e-4,   -1.4747049140e-5,
        2.9487055140e-4,   -1.8483155870e-4,    4.2616203430e-4,   -1.7956076770e-4,
        -1.0487979740e-3,    9.3550700690e-4,    3.7745822920e-4,   -4.7240295680e-4,
        9.4491872010e-4,   -1.8697801280e-3,   -6.4640177880e-4,    3.2934863120e-3,
        -1.1870980960e-3,   -4.7394598370e-4,    6.5074523450e-4,   -4.3425420300e-3,
        3.9003873240e-3,    4.4045397080e-3,   -6.1187194660e-3,    1.3679967960e-3,
        -1.5252598100e-3,   -3.9095366370e-3,    1.3784363870e-2,   -3.2901633530e-3,
        -1.3379056010e-2,    9.1043077410e-3,   -7.3734885080e-3,    7.0040063000e-3,
        3.0612513420e-2,   -4.6964682640e-2,   -2.9224805530e-2,    8.2755997780e-2,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.4797655800e-5,    1.2959531160e-4,    6.4797655800e-5,
    ];
    let expected_extremal_frequencies = vec![
        6.2500000930e-3,    1.8750000750e-2,    3.1250000000e-2,    4.3750002980e-2,
        5.6250005960e-2,    6.8750008940e-2,    8.1250011920e-2,    9.3750014900e-2,
        1.0625001790e-1,    1.1875002090e-1,    1.2500001490e-1,    1.3750000300e-1,
        1.4999999110e-1,    1.5624998510e-1,    1.6249997910e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2624999280e-1,    2.3249998690e-1,    2.3874998090e-1,
        2.5124996900e-1,    2.6374995710e-1,    2.7624994520e-1,    2.8874993320e-1,
        2.9499992730e-1,    3.0749991540e-1,    3.1374990940e-1,    3.1999990340e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8624998930e-1,    3.9249998330e-1,
        3.9874997740e-1,    4.1124996540e-1,    4.1749995950e-1,    4.2999994750e-1,
        4.4249993560e-1,    4.5499992370e-1,    4.6749991180e-1,    4.7999989990e-1,
        4.9249988790e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,2,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n81_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        -5.2403527660e-5,   -1.0170966560e-4,    1.2161411000e-4,   -5.6596763900e-6,
        1.7421552910e-5,    1.4226723580e-4,   -4.0837062990e-4,   -5.9987360150e-7,
        4.6876841220e-4,   -2.6070052990e-4,    3.1372689410e-4,   -1.4852517050e-4,
        -1.1231221720e-3,    1.1197620770e-3,    5.3801271130e-4,   -8.1489072180e-4,
        8.8142341700e-4,   -1.7182958550e-3,   -6.4550776730e-4,    3.5863565280e-3,
        -1.3406838990e-3,   -9.8915444690e-4,    9.4148941570e-4,   -4.0887002830e-3,
        3.8781191690e-3,    4.5517929830e-3,   -6.6594807430e-3,    1.1625442420e-3,
        -7.8897154890e-4,   -3.8943393160e-3,    1.3585054320e-2,   -3.3772196620e-3,
        -1.3919513670e-2,    9.5747746530e-3,   -6.6851810550e-3,    6.4637269820e-3,
        3.0381852760e-2,   -4.7053027900e-2,   -2.9344130310e-2,    8.3597242830e-2,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.0523905080e-4,    2.1047810150e-4,    1.0523905080e-4,
    ];
    let expected_extremal_frequencies = vec![
        4.1666668840e-3,    1.6666667540e-2,    2.9166666790e-2,    4.1666667910e-2,
        5.4166667160e-2,    6.6666670140e-2,    7.9166680570e-2,    9.1666691010e-2,
        1.0416670140e-1,    1.1666671190e-1,    1.2916670740e-1,    1.3750003280e-1,
        1.5000002090e-1,    1.5833334620e-1,    1.6250000890e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2416666150e-1,    2.3249998690e-1,    2.4499997500e-1,
        2.5749996300e-1,    2.6999995110e-1,    2.7833327650e-1,    2.9083326460e-1,
        3.0333325270e-1,    3.1166657810e-1,    3.1999990340e-1,    3.2416656610e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8416665790e-1,    3.8833332060e-1,
        3.9666664600e-1,    4.0916663410e-1,    4.2166662220e-1,    4.2999994750e-1,
        4.4249993560e-1,    4.5499992370e-1,    4.6749991180e-1,    4.7999989990e-1,
        4.9249988790e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,2,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n81_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        -4.4423530200e-5,   -1.2492266250e-4,    1.1277553860e-4,   -9.5496798170e-6,
        1.6660182150e-5,    1.9889586840e-4,   -4.0610443100e-4,   -7.5942094550e-5,
        4.6772856150e-4,   -2.5436689610e-4,    3.6037329120e-4,   -5.1691342380e-5,
        -1.2388413310e-3,    1.0018511680e-3,    6.2100368090e-4,   -7.7029818200e-4,
        9.9143665280e-4,   -1.7125642630e-3,   -9.3841459600e-4,    3.6000902760e-3,
        -1.1193525280e-3,   -9.5457339190e-4,    1.0146470740e-3,   -4.3173017910e-3,
        3.6039382680e-3,    4.8959446140e-3,   -6.4528072250e-3,    1.0586875720e-3,
        -8.5923098960e-4,   -4.2491657660e-3,    1.3673345560e-2,   -2.8343335730e-3,
        -1.4033846560e-2,    9.3285944310e-3,   -6.8116788750e-3,    6.2935985620e-3,
        3.0873913320e-2,   -4.6779654920e-2,   -2.9838647690e-2,    8.3453327420e-2,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.4130126510e-4,    2.8260253020e-4,    1.4130126510e-4,
    ];
    let expected_extremal_frequencies = vec![
        5.8333338240e-3,    1.8333332610e-2,    3.0833320690e-2,    4.3333310630e-2,
        5.5833298710e-2,    6.8333290520e-2,    8.0833278600e-2,    9.2499934140e-2,
        1.0499992220e-1,    1.1666657780e-1,    1.2833324070e-1,    1.3999989630e-1,
        1.5083321930e-1,    1.5999987720e-1,    1.6749987010e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2333332900e-1,    2.3166665430e-1,    2.4333330990e-1,
        2.5583329800e-1,    2.6833328600e-1,    2.7999994160e-1,    2.9166659710e-1,
        3.0249992010e-1,    3.1249991060e-1,    3.2166656850e-1,    3.2749989630e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8249999280e-1,    3.8916665320e-1,
        3.9833331110e-1,    4.0916663410e-1,    4.1999995710e-1,    4.3166661260e-1,
        4.4416660070e-1,    4.5666658880e-1,    4.6833324430e-1,    4.8083323240e-1,
        4.9333322050e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,2,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n81_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(81, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        -4.4184114810e-5,   -1.2486455670e-4,    1.1236281720e-4,   -9.5649302240e-6,
        1.6719313860e-5,    1.9940944910e-4,   -4.0580617500e-4,   -7.6841541160e-5,
        4.6745777950e-4,   -2.5386325430e-4,    3.6067239130e-4,   -5.0843402280e-5,
        -1.2395215450e-3,    1.0004590730e-3,    6.2200537650e-4,   -7.6963478930e-4,
        9.9178426900e-4,   -1.7123207220e-3,   -9.4067386820e-4,    3.5999566320e-3,
        -1.1171291120e-3,   -9.5466926000e-4,    1.0147530120e-3,   -4.3184561650e-3,
        3.6017380190e-3,    4.8986622130e-3,   -6.4511112870e-3,    1.0571709140e-3,
        -8.5933646190e-4,   -4.2513813820e-3,    1.3673672450e-2,   -2.8301160780e-3,
        -1.4035182070e-2,    9.3266889450e-3,   -6.8118921480e-3,    6.2920269560e-3,
        3.0877392740e-2,   -4.6777594830e-2,   -2.9842503370e-2,    8.3452910180e-2,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.4118697440e-4,    2.8237394870e-4,    1.4118697440e-4,
    ];
    let expected_extremal_frequencies = vec![
        6.2500005590e-3,    1.8749998880e-2,    3.1249986960e-2,    4.3749976900e-2,
        5.5468715730e-2,    6.7968726160e-2,    8.0468773840e-2,    9.2968821530e-2,
        1.0468761620e-1,    1.1718766390e-1,    1.2890645860e-1,    1.3984400030e-1,
        1.5078154210e-1,    1.6015657780e-1,    1.6718785460e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2312501070e-1,    2.3093754050e-1,    2.4343758820e-1,
        2.5593751670e-1,    2.6843732600e-1,    2.8015589710e-1,    2.9109323020e-1,
        3.0281180140e-1,    3.1296789650e-1,    3.2156151530e-1,    3.2781142000e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8234370950e-1,    3.8937485220e-1,
        3.9796847110e-1,    4.0890580420e-1,    4.1984313730e-1,    4.3234294650e-1,
        4.4406151770e-1,    4.5656132700e-1,    4.6906113620e-1,    4.8156094550e-1,
        4.9406075480e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,2,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n82_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        1.3303425480e-5,   -1.2472487290e-4,    2.2941530910e-5,    1.3873854190e-4,
        -6.3480423710e-5,    1.2433614760e-4,   -1.5876890390e-4,   -3.9458036190e-4,
        5.2713567860e-4,    1.5264510880e-4,   -2.3859215430e-4,    3.1038557060e-4,
        -9.8910334050e-4,   -9.7215524870e-6,    1.6626195280e-3,   -7.7468209200e-4,
        -1.4996661050e-4,    2.1442356230e-5,   -2.1410591430e-3,    2.5802394380e-3,
        1.7843304670e-3,   -3.1220491510e-3,    1.0281952560e-3,   -1.6839209710e-3,
        -1.2765183350e-3,    7.3377462100e-3,   -2.7757394130e-3,   -4.9898335710e-3,
        3.6527039480e-3,   -5.0939060750e-3,    5.6199743410e-3,    9.8572270940e-3,
        -1.5977956350e-2,   -6.7091779780e-4,    6.6406405530e-3,   -9.5037296410e-3,
        2.8255227950e-2,   -2.4215728040e-3,   -6.4260870220e-2,    3.8243528460e-2,
        6.8831153210e-2,
    ];
    let expected_deviations = vec![
        5.9048023100e-5,    1.1809604620e-4,    5.9048023100e-5,
    ];
    let expected_extremal_frequencies = vec![
        6.0975607480e-3,    1.8292682250e-2,    3.0487803740e-2,    4.2682923380e-2,
        5.4878048600e-2,    6.7073173820e-2,    7.9268299040e-2,    9.1463424270e-2,
        1.0365854950e-1,    1.1585367470e-1,    1.2804879250e-1,    1.3414634760e-1,
        1.4634145800e-1,    1.5243901310e-1,    1.5853656830e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.3219510910e-1,    2.3829266430e-1,    2.5048777460e-1,
        2.6268288490e-1,    2.6878044010e-1,    2.8097555040e-1,    2.9317066070e-1,
        3.0536577110e-1,    3.1146332620e-1,    3.1756088140e-1,    3.2365843650e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8609755040e-1,    3.9219510560e-1,
        3.9829266070e-1,    4.1048777100e-1,    4.1658532620e-1,    4.2878043650e-1,
        4.4097554680e-1,    4.5317065720e-1,    4.6536576750e-1,    4.7756087780e-1,
        4.8975598810e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,2,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n82_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        3.2057607310e-5,   -1.4966183520e-4,   -2.2350359360e-6,    1.2872007210e-4,
        -7.4555027820e-5,    1.8955431010e-4,   -1.1327698300e-4,   -4.7586168510e-4,
        4.6627031410e-4,    1.3730447970e-4,   -1.7024778940e-4,    4.7686847390e-4,
        -1.0441582420e-3,   -2.3129134210e-4,    1.6150820300e-3,   -6.6696206340e-4,
        8.4502353270e-5,    1.1300033660e-4,   -2.4842128620e-3,    2.3279555610e-3,
        1.9749745260e-3,   -2.8000213210e-3,    1.2141068000e-3,   -1.9688562020e-3,
        -1.7903175900e-3,    7.4345283210e-3,   -2.2426550280e-3,   -4.7432333230e-3,
        3.4173759630e-3,   -5.6781023740e-3,    5.4488447490e-3,    1.0496429170e-2,
        -1.5502987430e-2,   -9.4664283100e-4,    6.0543958100e-3,   -9.8146395760e-3,
        2.8735354540e-2,   -1.6952101140e-3,   -6.4409285780e-2,    3.7533074620e-2,
        6.8503513930e-2,
    ];
    let expected_deviations = vec![
        9.1203844930e-5,    1.8240768990e-4,    9.1203844930e-5,
    ];
    let expected_extremal_frequencies = vec![
        8.1300809980e-3,    2.0325202490e-2,    3.2520323990e-2,    4.4715445490e-2,
        5.6910566990e-2,    6.9105684760e-2,    8.1300795080e-2,    9.3495905400e-2,
        1.0569101570e-1,    1.1788612600e-1,    1.3008123640e-1,    1.3821130990e-1,
        1.5040642020e-1,    1.5853649380e-1,    1.6260153060e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2813007240e-1,    2.4032518270e-1,    2.4845525620e-1,
        2.6065036650e-1,    2.7284547690e-1,    2.8097555040e-1,    2.9317066070e-1,
        3.0536577110e-1,    3.1349584460e-1,    3.2162591810e-1,    3.2569095490e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8406503200e-1,    3.8813006880e-1,
        3.9626014230e-1,    4.0845525260e-1,    4.1658532620e-1,    4.2878043650e-1,
        4.4097554680e-1,    4.5317065720e-1,    4.6536576750e-1,    4.7756087780e-1,
        4.8975598810e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,2,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n82_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        4.6632056180e-5,   -1.4977518000e-4,   -1.8315564380e-5,    1.0414874120e-4,
        -9.1151421660e-5,    2.3373108710e-4,   -6.0002035750e-5,   -5.0760334120e-4,
        3.9535079850e-4,    1.0076972830e-4,   -1.1347369580e-4,    6.0733727880e-4,
        -1.0494338350e-3,   -3.9377470970e-4,    1.5284222320e-3,   -5.8773608180e-4,
        2.7661354400e-4,    1.9337321280e-4,   -2.7081668380e-3,    2.1060318690e-3,
        2.0763594660e-3,   -2.5292944630e-3,    1.3658544050e-3,   -2.1725385450e-3,
        -2.1659932100e-3,    7.4601150120e-3,   -1.8446333710e-3,   -4.5177442950e-3,
        3.2323747870e-3,   -6.1169029210e-3,    5.3163720290e-3,    1.0948066600e-2,
        -1.5118547720e-2,   -1.1342475190e-3,    5.5853733790e-3,   -1.0044839230e-2,
        2.9099974780e-2,   -1.1582951990e-3,   -6.4500510690e-2,    3.6992587150e-2,
        6.8235129120e-2,
    ];
    let expected_deviations = vec![
        1.0966344420e-4,    2.1932688830e-4,    1.0966344420e-4,
    ];
    let expected_extremal_frequencies = vec![
        6.5040653570e-3,    1.8699185920e-2,    3.1707305460e-2,    4.3902415780e-2,
        5.6910533460e-2,    6.9105647500e-2,    8.1300757830e-2,    9.3495868150e-2,
        1.0569097850e-1,    1.1788608880e-1,    1.2926819920e-1,    1.4065030220e-1,
        1.5121939780e-1,    1.6097548600e-1,    1.6747954490e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2325202820e-1,    2.3300811650e-1,    2.4764224890e-1,
        2.5983735920e-1,    2.7121946220e-1,    2.8260156510e-1,    2.9317066070e-1,
        3.0373975630e-1,    3.1349584460e-1,    3.2162591810e-1,    3.2812997700e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8243901730e-1,    3.8894307610e-1,
        3.9707314970e-1,    4.0682923790e-1,    4.1821134090e-1,    4.2959344390e-1,
        4.4097554680e-1,    4.5235764980e-1,    4.6455276010e-1,    4.7593486310e-1,
        4.8812997340e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,2,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n82_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(82, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        4.6503242630e-5,   -1.4996327810e-4,   -1.8060251020e-5,    1.0430192920e-4,
        -9.1126130430e-5,    2.3356157180e-4,   -6.0612488600e-5,   -5.0739606380e-4,
        3.9621035100e-4,    1.0082890370e-4,   -1.1382771480e-4,    6.0646689960e-4,
        -1.0498503690e-3,   -3.9217475570e-4,    1.5292691530e-3,   -5.8873789380e-4,
        2.7548446090e-4,    1.9244747820e-4,   -2.7066483160e-3,    2.1084530740e-3,
        2.0751797130e-3,   -2.5315438400e-3,    1.3650578910e-3,   -2.1714754400e-3,
        -2.1625827070e-3,    7.4602290990e-3,   -1.8485968470e-3,   -4.5191412790e-3,
        3.2338690940e-3,   -6.1138221060e-3,    5.3182216360e-3,    1.0943971570e-2,
        -1.5122057870e-2,   -1.1319434270e-3,    5.5885682810e-3,   -1.0042864830e-2,
        2.9097549620e-2,   -1.1635832490e-3,   -6.4499259000e-2,    3.6997318270e-2,
        6.8236559630e-2,
    ];
    let expected_deviations = vec![
        1.0985875270e-4,    2.1971750540e-4,    1.0985875270e-4,
    ];
    let expected_extremal_frequencies = vec![
        6.0975607480e-3,    1.9054876640e-2,    3.1249986960e-2,    4.4207293540e-2,
        5.6402403860e-2,    6.9359712300e-2,    8.1554822620e-2,    9.3749932940e-2,
        1.0594504330e-1,    1.1814015360e-1,    1.2957307700e-1,    1.4100599290e-1,
        1.5167671440e-1,    1.6082304720e-1,    1.6768279670e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2304877640e-1,    2.3295730350e-1,    2.4743899700e-1,
        2.5963410740e-1,    2.7106702330e-1,    2.8249993920e-1,    2.9317066070e-1,
        3.0384138230e-1,    3.1374990940e-1,    3.2213404770e-1,    3.2746940850e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8228657840e-1,    3.8838413360e-1,
        3.9753046630e-1,    4.0743899350e-1,    4.1810971500e-1,    4.2954263090e-1,
        4.4097554680e-1,    4.5240846280e-1,    4.6460357310e-1,    4.7603648900e-1,
        4.8823159930e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,2,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n128_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 2);

    let expected_impulse_response = vec![
        1.2790744680e-6,   -2.7626231260e-7,   -1.5970967980e-6,    2.8764920900e-6,
        -4.6946074690e-6,   -6.0764523370e-7,    1.1156159420e-5,   -9.6306903290e-6,
        2.5336012190e-6,    5.9041708480e-6,   -3.1893687260e-5,    3.1590847360e-5,
        2.4721282900e-5,   -5.0136230130e-5,    3.5236407710e-5,   -3.6996534620e-5,
        -4.6699729860e-5,    1.6689451880e-4,   -6.3508079620e-5,   -1.0941911020e-4,
        1.1170500510e-4,   -1.7324095820e-4,    1.4648560320e-4,    2.9751437250e-4,
        -4.6311944610e-4,    1.6086341930e-5,    1.8067515340e-4,   -2.8511189160e-4,
        6.9719948810e-4,   -1.5494442780e-4,   -1.0214069630e-3,    7.5896817730e-4,
        9.3522176030e-6,   -7.1712100180e-5,    1.1098806280e-3,   -1.7662234600e-3,
        -5.2644469540e-4,    2.1135741840e-3,   -9.0616638770e-4,    6.5258314130e-4,
        1.7211207890e-4,   -3.5940487870e-3,    2.7391794140e-3,    2.3538270030e-3,
        -2.8779399580e-3,    1.8083960750e-3,   -2.8110526040e-3,   -2.3966911250e-3,
        8.3677656950e-3,   -2.1709108260e-3,   -4.5717754400e-3,    3.3506951290e-3,
        -6.7748078150e-3,    5.8480123990e-3,    1.1287626810e-2,   -1.5571372580e-2,
        -9.0933591130e-4,    5.2844318560e-3,   -1.0151108730e-2,    2.9636332770e-2,
        -1.3872068380e-3,   -6.4475327730e-2,    3.6988418550e-2,    6.7875936630e-2,
    ];
    let expected_deviations = vec![
        8.1063478770e-7,    1.6212695750e-6,    8.1063478770e-7,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    1.1718750000e-2,    1.9531250000e-2,    2.7343750000e-2,
        3.5156250000e-2,    4.2968750000e-2,    5.0781250000e-2,    5.8593750000e-2,
        6.6406250000e-2,    7.4218750000e-2,    8.2031250000e-2,    8.9843750000e-2,
        9.7656250000e-2,    1.0546875000e-1,    1.0937500000e-1,    1.1718750000e-1,
        1.2500000000e-1,    1.3281250000e-1,    1.4062500000e-1,    1.4453125000e-1,
        1.5234375000e-1,    1.5625000000e-1,    1.6015625000e-1,    1.6406250000e-1,
        1.7000000180e-1,    2.1999999880e-1,    2.2390624880e-1,    2.2781249880e-1,
        2.3171874880e-1,    2.3562499880e-1,    2.4343749880e-1,    2.4734374880e-1,
        2.5515624880e-1,    2.6296874880e-1,    2.6687499880e-1,    2.7468749880e-1,
        2.8249999880e-1,    2.9031249880e-1,    2.9421874880e-1,    3.0203124880e-1,
        3.0984374880e-1,    3.1374999880e-1,    3.1765624880e-1,    3.2156249880e-1,
        3.2546874880e-1,    3.3000001310e-1,    3.7999999520e-1,    3.8390624520e-1,
        3.8781249520e-1,    3.9171874520e-1,    3.9562499520e-1,    3.9953124520e-1,
        4.0734374520e-1,    4.1515624520e-1,    4.2296874520e-1,    4.3078124520e-1,
        4.3468749520e-1,    4.4249999520e-1,    4.5031249520e-1,    4.5812499520e-1,
        4.6593749520e-1,    4.7374999520e-1,    4.8546874520e-1,    4.9328124520e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,2,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n128_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 3);

    let expected_impulse_response = vec![
        1.5916236860e-6,    1.1810902830e-7,   -1.6055358860e-6,    3.9779852160e-6,
        -6.4404953260e-6,   -2.4638084140e-6,    1.4044923770e-5,   -9.2815462270e-6,
        4.1542089090e-6,    6.5361605270e-6,   -4.2029514590e-5,    3.3787498980e-5,
        3.5166976890e-5,   -5.2007424530e-5,    3.8617894460e-5,   -4.8385627450e-5,
        -6.3160179710e-5,    1.9284289740e-4,   -5.0995247880e-5,   -1.2169260300e-4,
        1.0885830850e-4,   -2.0489099550e-4,    1.5573963170e-4,    3.5657006080e-4,
        -4.8030866310e-4,   -1.3951910660e-5,    1.6676864470e-4,   -3.1097547620e-4,
        7.7142188090e-4,   -1.0838118030e-4,   -1.1055431100e-3,    7.2935578650e-4,
        5.7719075810e-6,   -4.1966071880e-5,    1.2247086270e-3,   -1.8193833530e-3,
        -6.5717811230e-4,    2.1390861370e-3,   -8.6541735800e-4,    7.3931313820e-4,
        2.2226630240e-4,   -3.7750380580e-3,    2.6667623320e-3,    2.4791196920e-3,
        -2.8073906430e-3,    1.8674110760e-3,   -2.9064663690e-3,   -2.5988861450e-3,
        8.4620937710e-3,   -1.9876370210e-3,   -4.5591816310e-3,    3.2971028700e-3,
        -6.9545768200e-3,    5.7883011180e-3,    1.1535069910e-2,   -1.5463488180e-2,
        -1.0320274160e-3,    5.1569473000e-3,   -1.0254636410e-2,    2.9770743100e-2,
        -1.1405348780e-3,   -6.4558170740e-2,    3.6778662350e-2,    6.7820474510e-2,
    ];
    let expected_deviations = vec![
        1.7634916960e-6,    3.5269833920e-6,    1.7634916960e-6,
    ];
    let expected_extremal_frequencies = vec![
        2.6041667440e-3,    1.0416666980e-2,    1.8229166050e-2,    2.6041664180e-2,
        3.3854164180e-2,    4.1666667910e-2,    4.9479171630e-2,    5.7291675360e-2,
        6.5104171630e-2,    7.2916664180e-2,    8.0729156730e-2,    8.8541649280e-2,
        9.6354141830e-2,    1.0416663440e-1,    1.1197912690e-1,    1.1979161950e-1,
        1.2760411200e-1,    1.3541662690e-1,    1.4062497020e-1,    1.4843748510e-1,
        1.5364582840e-1,    1.5885417160e-1,    1.6406251490e-1,    1.6666668650e-1,
        1.7000000180e-1,    2.1999999880e-1,    2.2260417040e-1,    2.2520834210e-1,
        2.2781251370e-1,    2.3562502860e-1,    2.4083337190e-1,    2.4604171510e-1,
        2.5385421510e-1,    2.6166668530e-1,    2.6947915550e-1,    2.7468746900e-1,
        2.8249993920e-1,    2.9031240940e-1,    2.9552072290e-1,    3.0333319310e-1,
        3.0854150650e-1,    3.1635397670e-1,    3.2156229020e-1,    3.2416644690e-1,
        3.2677060370e-1,    3.3000001310e-1,    3.7999999520e-1,    3.8260415200e-1,
        3.8520830870e-1,    3.8781246540e-1,    3.9562493560e-1,    4.0083324910e-1,
        4.0864571930e-1,    4.1385403280e-1,    4.2166650300e-1,    4.2947897320e-1,
        4.3729144330e-1,    4.4510391350e-1,    4.5291638370e-1,    4.6072885390e-1,
        4.6854132410e-1,    4.7635379430e-1,    4.8416626450e-1,    4.9197873470e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
128,2,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n128_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 15);

    let expected_impulse_response = vec![
        2.0029322060e-6,    4.8407650870e-7,   -2.2919389270e-6,    4.8072006390e-6,
        -7.2402049230e-6,   -3.7581162360e-6,    1.6861789850e-5,   -1.0148207370e-5,
        3.1633526310e-6,    9.3240441860e-6,   -4.7520210500e-5,    3.4853135730e-5,
        4.2764269890e-5,   -5.8811445340e-5,    4.0603521480e-5,   -4.6858440330e-5,
        -7.5648058560e-5,    2.0747746750e-4,   -4.5215041610e-5,   -1.3933934680e-4,
        1.2036458070e-4,   -2.1329119040e-4,    1.4667742650e-4,    3.9211407420e-4,
        -4.9815489910e-4,   -3.4624797990e-5,    1.9245062140e-4,   -3.3626199000e-4,
        7.8772113190e-4,   -7.0341746320e-5,   -1.1660129530e-3,    7.3520484150e-4,
        3.7197198250e-5,   -8.0678873930e-5,    1.2787606100e-3,   -1.8222669610e-3,
        -7.4178387880e-4,    2.2034291180e-3,   -8.5726426910e-4,    7.0426391900e-4,
        2.9760837790e-4,   -3.8498435170e-3,    2.6200544090e-3,    2.5945033410e-3,
        -2.8565020770e-3,    1.8624196530e-3,   -2.8494875880e-3,   -2.7262454390e-3,
        8.5125230250e-3,   -1.8841519490e-3,   -4.6683279800e-3,    3.3434673680e-3,
        -6.9534145300e-3,    5.6736199190e-3,    1.1681765320e-2,   -1.5449571420e-2,
        -1.1536492970e-3,    5.2501242610e-3,   -1.0319096970e-2,    2.9735172170e-2,
        -9.6872076390e-4,   -6.4661748710e-2,    3.6715187130e-2,    6.7923419180e-2,
    ];
    let expected_deviations = vec![
        2.7288795080e-6,    5.4577590160e-6,    2.7288795080e-6,
    ];
    let expected_extremal_frequencies = vec![
        3.6458333490e-3,    1.1458336380e-2,    1.9270835440e-2,    2.7083327990e-2,
        3.4895822410e-2,    4.2708314960e-2,    5.0520807500e-2,    5.8333300050e-2,
        6.6145792600e-2,    7.3958285150e-2,    8.1770777700e-2,    8.9583270250e-2,
        9.6874929960e-2,    1.0468742250e-1,    1.1249991510e-1,    1.1979157480e-1,
        1.2708325680e-1,    1.3437502090e-1,    1.4166678490e-1,    1.4843770860e-1,
        1.5520863230e-1,    1.6093787550e-1,    1.6562543810e-1,    1.6875047980e-1,
        1.7000000180e-1,    2.1999999880e-1,    2.2104167940e-1,    2.2416672110e-1,
        2.2885428370e-1,    2.3406268660e-1,    2.4031277000e-1,    2.4708369370e-1,
        2.5385451320e-1,    2.6166689400e-1,    2.6843762400e-1,    2.7572917940e-1,
        2.8302073480e-1,    2.9031229020e-1,    2.9708302020e-1,    3.0385375020e-1,
        3.1010365490e-1,    3.1635355950e-1,    3.2156181340e-1,    3.2624924180e-1,
        3.2885336880e-1,    3.3000001310e-1,    3.7999999520e-1,    3.8104164600e-1,
        3.8416659830e-1,    3.8833320140e-1,    3.9406228070e-1,    4.0031218530e-1,
        4.0708291530e-1,    4.1437447070e-1,    4.2166602610e-1,    4.2895758150e-1,
        4.3624913690e-1,    4.4406151770e-1,    4.5187389850e-1,    4.5968627930e-1,
        4.6801948550e-1,    4.7583186630e-1,    4.8364424710e-1,    4.9197745320e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
128,2,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn differentiator_bandpass_n128_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(128, FilterType::Differentiator, &bands, 16);

    let expected_impulse_response = vec![
        2.0054289960e-6,    4.8159063230e-7,   -2.3060733840e-6,    4.8059728210e-6,
        -7.2161756180e-6,   -3.7568677270e-6,    1.6848258380e-5,   -1.0169380400e-5,
        3.1379563550e-6,    9.3889193520e-6,   -4.7475423340e-5,    3.4770386260e-5,
        4.2730531280e-5,   -5.8819914560e-5,    4.0659761000e-5,   -4.6700464740e-5,
        -7.5783835200e-5,    2.0728536770e-4,   -4.5097272960e-5,   -1.3926661630e-4,
        1.2052253440e-4,   -2.1329720040e-4,    1.4618836580e-4,    3.9221046610e-4,
        -4.9772870260e-4,   -3.4652373870e-5,    1.9249452450e-4,   -3.3672971770e-4,
        7.8737037260e-4,   -6.9431553130e-5,   -1.1658119040e-3,    7.3470018110e-4,
        3.7044941560e-5,   -8.1259531730e-5,    1.2794454810e-3,   -1.8211232960e-3,
        -7.4283598220e-4,    2.2027981470e-3,   -8.5712515280e-4,    7.0431351200e-4,
        2.9915216150e-4,   -3.8501853120e-3,    2.6180113200e-3,    2.5951161510e-3,
        -2.8558296620e-3,    1.8629792610e-3,   -2.8486964290e-3,   -2.7286205440e-3,
        8.5118087010e-3,   -1.8818264360e-3,   -4.6681966630e-3,    3.3435761000e-3,
        -6.9543262940e-3,    5.6713628580e-3,    1.1684050780e-2,   -1.5447672460e-2,
        -1.1552399960e-3,    5.2497382280e-3,   -1.0320504200e-2,    2.9735449700e-2,
        -9.6526369450e-4,   -6.4662881200e-2,    3.6712951960e-2,    6.7923761900e-2,
    ];
    let expected_deviations = vec![
        2.7248847800e-6,    5.4497695600e-6,    2.7248847800e-6,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    1.1718750000e-2,    1.9531250000e-2,    2.7343750000e-2,
        3.5156250000e-2,    4.2968750000e-2,    5.0781250000e-2,    5.8593750000e-2,
        6.5917968750e-2,    7.3730468750e-2,    8.1542968750e-2,    8.9355468750e-2,
        9.7167968750e-2,    1.0449218750e-1,    1.1230468750e-1,    1.1962890620e-1,
        1.2744140620e-1,    1.3476562500e-1,    1.4160156250e-1,    1.4843750000e-1,
        1.5527343750e-1,    1.6064453120e-1,    1.6552734380e-1,    1.6894531250e-1,
        1.7000000180e-1,    2.1999999880e-1,    2.2097656130e-1,    2.2390624880e-1,
        2.2878906130e-1,    2.3416015510e-1,    2.4050781130e-1,    2.4734374880e-1,
        2.5417968630e-1,    2.6150390510e-1,    2.6833984260e-1,    2.7566406130e-1,
        2.8298828010e-1,    2.9031249880e-1,    2.9714843630e-1,    3.0398437380e-1,
        3.1033203010e-1,    3.1619140510e-1,    3.2156249880e-1,    3.2595703010e-1,
        3.2888671760e-1,    3.3000001310e-1,    3.7999999520e-1,    3.8097655770e-1,
        3.8390624520e-1,    3.8878905770e-1,    3.9416015150e-1,    4.0050780770e-1,
        4.0734374520e-1,    4.1417968270e-1,    4.2150390150e-1,    4.2882812020e-1,
        4.3664062020e-1,    4.4396483900e-1,    4.5177733900e-1,    4.6007812020e-1,
        4.6789062020e-1,    4.7570312020e-1,    4.8400390150e-1,    4.9181640150e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n3_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        5.7735025880e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        0.0000000000e0,
    ];
    let expected_extremal_frequencies = vec![
        1.6666667160e-1,    3.3333334330e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n3_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        8.3164674040e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.5418183800e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.3333333430e-1,    4.6666666870e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n3_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        8.3675682540e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.7351365090e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.5000000000e-1,    4.6875000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n3_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        9.1983139510e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.3966279030e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3888888990e-2,    2.5000005960e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n3_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        9.2259019610e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.4351778030e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3513513840e-2,    2.5675681230e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n4_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        1.9891238210e-1,    6.1991435290e-1,
    ];
    let expected_deviations = vec![
        1.5799595420e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.2500000000e-1,    2.5000000000e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n4_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        2.7480390670e-1,    6.3293528560e-1,
    ];
    let expected_deviations = vec![
        2.8373715280e-1,
    ];
    let expected_extremal_frequencies = vec![
        8.3333335820e-2,    2.5000000000e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n4_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        5.2948409320e-1,    6.4612942930e-1,
    ];
    let expected_deviations = vec![
        7.6670920850e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.6666667540e-2,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n4_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        5.3596317770e-1,    6.4631867410e-1,
    ];
    let expected_deviations = vec![
        7.7928906680e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.5625000000e-2,    2.0312500000e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n4_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        5.9519350530e-1,    6.4826273920e-1,
    ];
    let expected_deviations = vec![
        8.9386141300e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.9444444960e-3,    2.0138895510e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n4_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        5.9652704000e-1,    6.4824801680e-1,
    ];
    let expected_deviations = vec![
        8.9655786750e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.7567569200e-3,    1.9594593350e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n11_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        1.0916339610e-1,   -7.6799002710e-9,    1.8790662290e-1,    7.4643180530e-9,
        6.2456035610e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        9.1634325680e-2,
    ];
    let expected_extremal_frequencies = vec![
        5.0000000750e-2,    1.0000000150e-1,    1.5000000600e-1,    2.5000000000e-1,
        3.5000002380e-1,    4.0000003580e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n11_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.7559368910e-1,   -6.0996363520e-8,    1.9724638760e-1,    1.5754384460e-8,
        6.3830280300e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.9856448470e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.3333335070e-2,    6.6666670140e-2,    1.6666667160e-1,    2.3333333430e-1,
        3.3333334330e-1,    4.3333333730e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n11_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        4.3815743920e-1,    1.1067613310e-7,    2.2273236510e-1,    6.4361096010e-8,
        6.4082723860e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.0830219980e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.6666668280e-3,    6.0000006110e-2,    1.5333332120e-1,    2.5333324070e-1,
        3.4666648510e-1,    4.3999972940e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n11_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        4.4534718990e-1,   -8.3883463730e-8,    2.2272896770e-1,    5.5205951810e-8,
        6.3922393320e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.2368454930e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.2500000930e-3,    6.2500007450e-2,    1.5624998510e-1,    2.4999989570e-1,
        3.4374982120e-1,    4.3749973180e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n11_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        5.1703411340e-1,   -4.9542880020e-8,    2.2564560170e-1,    1.2712546040e-7,
        6.4056038860e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.6389756200e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.7777778450e-3,    6.1111111190e-2,    1.5277785060e-1,    2.5000032780e-1,
        3.4722280500e-1,    4.3888971210e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n11_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        5.1893544200e-1,    1.2617556420e-7,    2.2582715750e-1,    7.8670254310e-8,
        6.4091724160e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.6721634860e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.7027027680e-3,    6.2162149700e-2,    1.5405406060e-1,    2.5135120750e-1,
        3.4594616290e-1,    4.3783840540e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n12_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        8.9716330170e-2,    5.9545964000e-2,    7.7454447750e-2,    1.1655704680e-1,
        2.0547857880e-1,    6.3434439900e-1,
    ];
    let expected_deviations = vec![
        1.2440403550e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.1666667910e-2,    8.3333335820e-2,    1.6666667160e-1,    2.5000000000e-1,
        3.3333331350e-1,    4.1666662690e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n12_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.4899484810e-1,    7.2529882190e-2,    8.6670055990e-2,    1.2721651790e-1,
        2.1358826760e-1,    6.3712835310e-1,
    ];
    let expected_deviations = vec![
        2.2475649420e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.7777777980e-2,    8.3333335820e-2,    1.3888889550e-1,    2.2222222390e-1,
        3.3333337310e-1,    4.1666674610e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n12_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        4.0253674980e-1,    7.9136461020e-2,    9.5899254080e-2,    1.3081890340e-1,
        2.1406924720e-1,    6.3714891670e-1,
    ];
    let expected_deviations = vec![
        7.3080116510e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.5555556900e-3,    5.5555555970e-2,    1.4444445070e-1,    2.2777777910e-1,
        3.2222241160e-1,    4.1111153360e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n12_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        4.0918114780e-1,    7.8522384170e-2,    9.6296519040e-2,    1.3132214550e-1,
        2.1414369340e-1,    6.3721895220e-1,
    ];
    let expected_deviations = vec![
        7.4511456490e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.2083334890e-3,    5.7291660460e-2,    1.4062500000e-1,    2.2916658220e-1,
        3.1770834330e-1,    4.1145852210e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n12_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        4.7491100430e-1,    8.0442994830e-2,    9.7735881810e-2,    1.3186383250e-1,
        2.1485960480e-1,    6.3750004770e-1,
    ];
    let expected_deviations = vec![
        8.7539839740e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.3148148320e-3,    5.5555544790e-2,    1.4120368660e-1,    2.2916688020e-1,
        3.1944444780e-1,    4.0972188120e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n12_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        4.7650849820e-1,    8.0546319480e-2,    9.7820520400e-2,    1.3185328250e-1,
        2.1476602550e-1,    6.3744664190e-1,
    ];
    let expected_deviations = vec![
        8.7849658730e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.2522523070e-3,    5.6306295100e-2,    1.4189183710e-1,    2.2972962260e-1,
        3.1981965900e-1,    4.0990969540e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n16_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        8.0567307770e-2,    4.3119259180e-2,    4.7118321060e-2,    6.0222670440e-2,
        8.2171127200e-2,    1.2084665890e-1,    2.0825800300e-1,    6.3529503350e-1,
    ];
    let expected_deviations = vec![
        1.1726237090e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.1250000000e-2,    6.2500000000e-2,    1.2500000000e-1,    1.8750000000e-1,
        2.5000000000e-1,    3.1250000000e-1,    3.7500000000e-1,    4.3750000000e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n16_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.3639324900e-1,    4.8712655900e-2,    5.3531497720e-2,    7.2740837930e-2,
        9.2477023600e-2,    1.2655347590e-1,    2.1254098420e-1,    6.3740909100e-1,
    ];
    let expected_deviations = vec![
        2.1905307470e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.0833333950e-2,    4.1666667910e-2,    1.0416667160e-1,    1.6666665670e-1,
        2.2916664180e-1,    3.1250000000e-1,    3.7500002980e-1,    4.3750005960e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n16_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        3.8976261020e-1,    5.5416017770e-2,    6.2209367750e-2,    7.4133276940e-2,
        9.3731880190e-2,    1.2940317390e-1,    2.1345859770e-1,    6.3703131680e-1,
    ];
    let expected_deviations = vec![
        7.2635513540e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.1666668840e-3,    4.1666667910e-2,    1.0416670140e-1,    1.7083333430e-1,
        2.3333327470e-1,    2.9999989270e-1,    3.6666649580e-1,    4.3333309890e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n16_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        3.9673832060e-1,    5.5804848670e-2,    6.2963694330e-2,    7.4438214300e-2,
        9.3600928780e-2,    1.2914431100e-1,    2.1328538660e-1,    6.3697898390e-1,
    ];
    let expected_deviations = vec![
        7.4044090510e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    4.2968750000e-2,    1.0546875000e-1,    1.6796875000e-1,
        2.3437500000e-1,    3.0078125000e-1,    3.6718750000e-1,    4.3359375000e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n16_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        4.6358865500e-1,    5.6776702400e-2,    6.3965380190e-2,    7.5501620770e-2,
        9.4319581990e-2,    1.2970834970e-1,    2.1360200640e-1,    6.3704109190e-1,
    ];
    let expected_deviations = vec![
        8.7289261820e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7361111240e-3,    4.1666679080e-2,    1.0416670890e-1,    1.6840265690e-1,
        2.3437462750e-1,    3.0034661290e-1,    3.6805468800e-1,    4.3402665850e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n16_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        4.6512380240e-1,    5.6754082440e-2,    6.4000308510e-2,    7.5455129150e-2,
        9.4340205190e-2,    1.2973469500e-1,    2.1358305220e-1,    6.3706231120e-1,
    ];
    let expected_deviations = vec![
        8.7607944010e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.6891892300e-3,    4.0540542450e-2,    1.0304050890e-1,    1.6891904180e-1,
        2.3479767140e-1,    3.0067628620e-1,    3.6655491590e-1,    4.3412274120e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n47_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        6.3395671550e-2,    1.0012001890e-6,    2.6511177420e-2,   -3.1720219340e-8,
        2.8263919060e-2,   -1.4735292100e-7,    3.1315639620e-2,   -1.7915118630e-7,
        3.5900637510e-2,   -2.6115006340e-7,    4.2444318530e-2,   -2.2964796640e-7,
        5.1728099580e-2,   -2.4531945540e-7,    6.5279603000e-2,   -2.1183812750e-7,
        8.6435467000e-2,   -2.0964844790e-7,    1.2395906450e-1,   -1.9965477800e-7,
        2.1013009550e-1,   -1.4179624940e-7,    6.3591814040e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.0085048530e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0869565420e-2,    2.1739130840e-2,    3.2608695330e-2,    5.4347828030e-2,
        7.6086953280e-2,    9.7826078530e-2,    1.1956520380e-1,    1.4130432900e-1,
        1.6304345430e-1,    1.8478257950e-1,    2.0652170480e-1,    2.2826083000e-1,
        2.4999995530e-1,    2.7173909540e-1,    2.9347822070e-1,    3.1521734600e-1,
        3.3695647120e-1,    3.5869559650e-1,    3.8043472170e-1,    4.0217384700e-1,
        4.2391297220e-1,    4.4565209750e-1,    4.6739122270e-1,    4.7826078530e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n47_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.2683261930e-1,   -2.1486280840e-7,    2.1344587210e-2,    1.8793056710e-7,
        2.9186218980e-2,    8.4905181550e-8,    3.4716427330e-2,    1.7494478750e-7,
        4.0901601310e-2,    3.2687367480e-8,    4.8015147450e-2,    9.7649518690e-8,
        5.7374626400e-2,    8.2363158070e-8,    7.0449233060e-2,    2.3935882610e-8,
        9.0833574530e-2,    1.7126694730e-8,    1.2726527450e-1,    5.1532424550e-8,
        2.1220427750e-1,   -2.8767601630e-8,    6.3660824300e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.0593385400e-1,
    ];
    let expected_extremal_frequencies = vec![
        7.2463769470e-3,    1.4492753890e-2,    3.6231882870e-2,    5.7971008120e-2,
        7.9710133370e-2,    1.0144925860e-1,    1.2318838390e-1,    1.4492751660e-1,
        1.6666664180e-1,    1.8840576710e-1,    2.1014489230e-1,    2.3188401760e-1,
        2.5362315770e-1,    2.6811590790e-1,    2.8985503320e-1,    3.1159415840e-1,
        3.3333328370e-1,    3.5507240890e-1,    3.7681153420e-1,    3.9855065940e-1,
        4.2028978470e-1,    4.4202890990e-1,    4.6376803520e-1,    4.8550716040e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n47_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        3.7344983220e-1,    1.8318286270e-7,    3.4382343290e-2,    3.2869650110e-7,
        3.7910282610e-2,   -3.6888536670e-7,    4.1053175930e-2,    1.0726218990e-6,
        4.5457571740e-2,    6.8444671800e-7,    5.1314771180e-2,    1.3223193490e-6,
        5.9715807440e-2,    9.4130109570e-7,    7.2090566160e-2,    1.0749895410e-6,
        9.1938734050e-2,    6.7846031020e-7,    1.2797224520e-1,    5.9519197750e-7,
        2.1259176730e-1,    2.4814517020e-7,    6.3673412800e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.1110880370e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.4492754130e-3,    1.4492755760e-2,    3.3333338800e-2,    5.5072475220e-2,
        7.6811574400e-2,    9.8550654950e-2,    1.2028973550e-1,    1.4202882350e-1,
        1.6376790400e-1,    1.8550698460e-1,    2.0724606510e-1,    2.2898514570e-1,
        2.5072422620e-1,    2.7101424340e-1,    2.9275354740e-1,    3.1449285150e-1,
        3.3623215560e-1,    3.5797145960e-1,    3.7971076370e-1,    4.0145006780e-1,
        4.2318937180e-1,    4.4492867590e-1,    4.6666798000e-1,    4.8550871010e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n47_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        3.8051763180e-1,    2.2099106540e-8,    3.5133391620e-2,    1.5371351480e-7,
        3.7069410090e-2,   -2.6320003600e-7,    4.0293723340e-2,    1.0248743370e-7,
        4.4890344140e-2,    2.3413537060e-7,    5.1153302190e-2,    3.5672218250e-7,
        5.9789359570e-2,    2.3209378240e-7,    7.2326362130e-2,    3.7594702460e-8,
        9.2163324360e-2,   -1.7809628620e-7,    1.2815856930e-1,   -2.4539315290e-7,
        2.1268343930e-1,   -1.7862632260e-7,    6.3677358630e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.2654944660e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3586956770e-3,    1.3586955150e-2,    3.3967386930e-2,    5.5706512180e-2,
        7.6086945830e-2,    9.7826071080e-2,    1.1956519630e-1,    1.4130440350e-1,
        1.6304364800e-1,    1.8478289250e-1,    2.0652213690e-1,    2.2826138140e-1,
        2.5000062580e-1,    2.7173963190e-1,    2.9347863790e-1,    3.1521764400e-1,
        3.3695665000e-1,    3.5869565610e-1,    3.8043466210e-1,    4.0217366810e-1,
        4.2391267420e-1,    4.4429299240e-1,    4.6603199840e-1,    4.8641231660e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n47_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        4.5014631750e-1,    4.1960109340e-7,    3.5819441080e-2,    1.9342706990e-8,
        3.8045495750e-2,    2.8948971930e-7,    4.1410863400e-2,    9.2671734820e-7,
        4.5854270460e-2,   -1.0496023610e-8,    5.1784634590e-2,   -5.3153803490e-7,
        6.0134351250e-2,   -1.1229417400e-7,    7.2544813160e-2,    7.3854937450e-7,
        9.2359840870e-2,    8.3591544350e-7,    1.2833321090e-1,    2.5385133990e-7,
        2.1279823780e-1,   -1.7630077310e-7,    6.3681173320e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.6526739600e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.0386472610e-4,    1.3285025950e-2,    3.3816412090e-2,    5.4951716210e-2,
        7.6690800490e-2,    9.8429836330e-2,    1.1956501010e-1,    1.4130423960e-1,
        1.6304354370e-1,    1.8478284780e-1,    2.0652215180e-1,    2.2826145590e-1,
        2.5000074510e-1,    2.7173951270e-1,    2.9347828030e-1,    3.1521704790e-1,
        3.3695581560e-1,    3.5869458320e-1,    3.8043335080e-1,    4.0156826380e-1,
        4.2330703140e-1,    4.4504579900e-1,    4.6618071200e-1,    4.8671177030e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n47_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        4.5187079910e-1,    4.3266794590e-7,    3.5801410680e-2,    1.1915340110e-6,
        3.8331210610e-2,    2.7907310600e-7,    4.1392147540e-2,    1.0647632960e-6,
        4.5778751370e-2,   -4.3986733540e-7,    5.1741063590e-2,   -1.4517445380e-7,
        6.0166835780e-2,   -6.9301870550e-7,    7.2538793090e-2,   -1.7454465250e-7,
        9.2300832270e-2,   -2.0834033880e-7,    1.2824076410e-1,    1.6965236680e-7,
        2.1274077890e-1,    7.6241690290e-8,    6.3678836820e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.6859315630e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.8754405470e-4,    1.3513515700e-2,    3.4077547490e-2,    5.5229179560e-2,
        7.6380811630e-2,    9.8119989040e-2,    1.1985916640e-1,    1.4159813520e-1,
        1.6333703700e-1,    1.8507593870e-1,    2.0681484040e-1,    2.2855374220e-1,
        2.5029265880e-1,    2.7144455910e-1,    2.9318401220e-1,    3.1492346530e-1,
        3.3666291830e-1,    3.5840237140e-1,    3.8014182450e-1,    4.0188127760e-1,
        4.2362073060e-1,    4.4477263090e-1,    4.6592453120e-1,    4.8648887870e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n48_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        6.3660763200e-2,    1.7099395390e-2,    9.6304416660e-3,    1.0411277410e-2,
        1.3397492470e-2,    1.5391752120e-2,    1.5774950390e-2,    1.5737831590e-2,
        1.6552314160e-2,    1.8401071430e-2,    2.0603299140e-2,    2.2610694170e-2,
        2.4585306640e-2,    2.7108550070e-2,    3.0567765240e-2,    3.4972906110e-2,
        4.0311336520e-2,    4.6992987390e-2,    5.6025147440e-2,    6.9153636690e-2,
        8.9758992200e-2,    1.2654566760e-1,    2.1177530290e-1,    6.3648188110e-1,
    ];
    let expected_deviations = vec![
        1.0347357390e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0416666980e-2,    2.0833333950e-2,    3.1250000000e-2,    5.2083335820e-2,
        7.2916664180e-2,    9.3749992550e-2,    1.1458332090e-1,    1.3541665670e-1,
        1.5625000000e-1,    1.7708334330e-1,    2.0833335820e-1,    2.2916670140e-1,
        2.5000002980e-1,    2.7083334330e-1,    2.9166665670e-1,    3.1249997020e-1,
        3.3333328370e-1,    3.5416659710e-1,    3.7499991060e-1,    3.9583322410e-1,
        4.1666653750e-1,    4.3749985100e-1,    4.5833316450e-1,    4.7916647790e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n48_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.1539263280e-1,    1.0872013870e-2,    1.0981559750e-2,    1.8044695260e-2,
        1.4981463550e-2,    1.5369579200e-2,    1.8722608690e-2,    1.8765777350e-2,
        1.9583195450e-2,    2.2230342030e-2,    2.3435503240e-2,    2.4964213370e-2,
        2.7904659510e-2,    3.0341953040e-2,    3.3242106440e-2,    3.7643373010e-2,
        4.2521864180e-2,    4.8822075130e-2,    5.8035790920e-2,    7.0805609230e-2,
        9.0840697290e-2,    1.2745165820e-1,    2.1223324540e-1,    6.3651561740e-1,
    ];
    let expected_deviations = vec![
        2.1209950750e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.9444444960e-3,    1.3888888990e-2,    3.4722223880e-2,    5.5555555970e-2,
        7.6388895510e-2,    9.7222238780e-2,    1.1805558200e-1,    1.3888892530e-1,
        1.5972226860e-1,    1.8055561180e-1,    2.0138895510e-1,    2.2222229840e-1,
        2.4305564170e-1,    2.6388898490e-1,    2.8472232820e-1,    3.0555567150e-1,
        3.3333346250e-1,    3.5416680570e-1,    3.7500014900e-1,    3.9583349230e-1,
        4.1666683550e-1,    4.3750017880e-1,    4.5833352210e-1,    4.7916686530e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n48_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        3.6702317000e-1,    1.6401022670e-2,    1.6866207120e-2,    1.7492979760e-2,
        1.8055886030e-2,    1.8845379350e-2,    1.9779771570e-2,    2.0787656310e-2,
        2.1959006790e-2,    2.3326456550e-2,    2.4836897850e-2,    2.6572465900e-2,
        2.8659403320e-2,    3.0997335910e-2,    3.3750772480e-2,    3.8532197480e-2,
        4.3061017990e-2,    4.9455940720e-2,    5.8473765850e-2,    7.1057140830e-2,
        9.1144502160e-2,    1.2749135490e-1,    2.1228992940e-1,    6.3663446900e-1,
    ];
    let expected_deviations = vec![
        7.1661967040e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3888889230e-3,    1.3888888990e-2,    3.3333335070e-2,    5.4166667160e-2,
        7.5000032780e-2,    9.5833420750e-2,    1.1805570130e-1,    1.3888908920e-1,
        1.5972247720e-1,    1.8055586520e-1,    2.0277814570e-1,    2.2361153360e-1,
        2.4444492160e-1,    2.6666700840e-1,    2.8750017290e-1,    3.0833333730e-1,
        3.3055537940e-1,    3.5138854380e-1,    3.7222170830e-1,    3.9305487280e-1,
        4.1527691480e-1,    4.3611007930e-1,    4.5694324370e-1,    4.7916528580e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n48_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        3.7432870270e-1,    1.6783803700e-2,    1.7212271690e-2,    1.7723500730e-2,
        1.8333524470e-2,    1.9046813250e-2,    1.9864261150e-2,    2.0794093610e-2,
        2.1856963630e-2,    2.3088395600e-2,    2.4541616440e-2,    2.6285767560e-2,
        2.8413355350e-2,    3.1039178370e-2,    3.4328281880e-2,    3.8488805290e-2,
        4.2781770230e-2,    4.9408555030e-2,    5.8236956600e-2,    7.1067273620e-2,
        9.1241836550e-2,    1.2756562230e-1,    2.1236467360e-1,    6.3667595390e-1,
    ];
    let expected_deviations = vec![
        7.3108166460e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3020833720e-3,    1.3020832090e-2,    3.3854167910e-2,    5.4687481370e-2,
        7.5520828370e-2,    9.6354201440e-2,    1.1718757450e-1,    1.3802087310e-1,
        1.6015620530e-1,    1.8098945920e-1,    2.0182271300e-1,    2.2395804520e-1,
        2.4479129910e-1,    2.6562473180e-1,    2.8776031730e-1,    3.0859380960e-1,
        3.2942730190e-1,    3.5156288740e-1,    3.7239637970e-1,    3.9322987200e-1,
        4.1536545750e-1,    4.3619894980e-1,    4.5703244210e-1,    4.7916802760e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,1,36
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n48_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        4.4251361490e-1,    1.7065823080e-2,    1.7511695620e-2,    1.8029093740e-2,
        1.8631219860e-2,    1.9340217110e-2,    2.0164608960e-2,    2.1116912360e-2,
        2.2214710710e-2,    2.3485183720e-2,    2.4980485440e-2,    2.6803672310e-2,
        2.8722882270e-2,    3.1428158280e-2,    3.4512996670e-2,    3.8315415380e-2,
        4.3171703820e-2,    4.9576878550e-2,    5.8370709420e-2,    7.1137547490e-2,
        9.1264724730e-2,    1.2756872180e-1,    2.1236312390e-1,    6.3662433620e-1,
    ];
    let expected_deviations = vec![
        8.6787027120e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.7870370800e-4,    1.3310182840e-2,    3.2986100760e-2,    5.3819488730e-2,
        7.5231499970e-2,    9.6064753830e-2,    1.1747670920e-1,    1.3830997050e-1,
        1.5972192590e-1,    1.8113388120e-1,    2.0254583660e-1,    2.2337909040e-1,
        2.4479104580e-1,    2.6620301600e-1,    2.8761497140e-1,    3.0844822530e-1,
        3.2986018060e-1,    3.5127213600e-1,    3.7268409130e-1,    3.9351734520e-1,
        4.1492930050e-1,    4.3634125590e-1,    4.5775321130e-1,    4.7858646510e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,1,37
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n48_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        4.4414520260e-1,    1.7079293730e-2,    1.7513632770e-2,    1.8029332160e-2,
        1.8644869330e-2,    1.9350290300e-2,    2.0172417160e-2,    2.1124780180e-2,
        2.2220849990e-2,    2.3497462270e-2,    2.4953663350e-2,    2.6853799820e-2,
        2.8800547120e-2,    3.1311571600e-2,    3.4411787990e-2,    3.8290798660e-2,
        4.3186843400e-2,    4.9615740780e-2,    5.8409810070e-2,    7.1158111100e-2,
        9.1261029240e-2,    1.2752819060e-1,    2.1235883240e-1,    6.3666498660e-1,
    ];
    let expected_deviations = vec![
        8.7115955350e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.6306307670e-4,    1.2950448320e-2,    3.3220708370e-2,    5.4054029290e-2,
        7.4887350200e-2,    9.6283733840e-2,    1.1711705480e-1,    1.3851360980e-1,
        1.5991027650e-1,    1.8130694330e-1,    2.0214053990e-1,    2.2353720660e-1,
        2.4493387340e-1,    2.6633009310e-1,    2.8716313840e-1,    3.0855923890e-1,
        3.2995533940e-1,    3.5135144000e-1,    3.7218448520e-1,    3.9358058570e-1,
        4.1497668620e-1,    4.3637278680e-1,    4.5720583200e-1,    4.7860193250e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n81_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        4.5916454840e-5,    6.5278694030e-2,   -2.5428720620e-5,    2.0330533390e-2,
        -2.7854501240e-5,    1.3196416200e-2,   -8.2034384830e-6,    1.3408966360e-2,
        1.3679566110e-5,    1.6111508010e-2,    1.7748940080e-5,    1.8799722190e-2,
        3.4095010050e-6,    2.0521953700e-2,   -1.1912391530e-5,    2.1684303880e-2,
        -1.2042675730e-5,    2.3259416220e-2,    1.3683384170e-6,    2.5907874110e-2,
        1.2525581040e-5,    2.9652431610e-2,    9.3084845500e-6,    3.4205645320e-2,
        -3.6982264650e-6,    3.9527773860e-2,   -1.1507145250e-5,    4.6148657800e-2,
        -5.7578117780e-6,    5.5182009940e-2,    6.5912477110e-6,    6.8414509300e-2,
        1.1166735930e-5,    8.9216828350e-2,    3.0059345590e-6,    1.2623560430e-1,
        -8.3599288700e-6,    2.1165138480e-1,   -9.9255012170e-6,    6.3645535710e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.0401435200e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.2500000930e-3,    1.2500000190e-2,    1.8750000750e-2,    3.1250000000e-2,
        4.3750002980e-2,    5.6250005960e-2,    6.8750008940e-2,    8.1250011920e-2,
        1.0000001640e-1,    1.1250001940e-1,    1.2500001490e-1,    1.3750000300e-1,
        1.4999999110e-1,    1.6249997910e-1,    1.7499996720e-1,    1.8749995530e-1,
        1.9999994340e-1,    2.1249993150e-1,    2.2499991950e-1,    2.3749990760e-1,
        2.4999989570e-1,    2.6249989870e-1,    2.7499988680e-1,    2.8749987480e-1,
        2.9999986290e-1,    3.1249985100e-1,    3.2499983910e-1,    3.3749982710e-1,
        3.4999981520e-1,    3.6249980330e-1,    3.7499979140e-1,    3.8749977950e-1,
        3.9999976750e-1,    4.1874974970e-1,    4.3124973770e-1,    4.4374972580e-1,
        4.5624971390e-1,    4.6874970200e-1,    4.8124969010e-1,    4.8749968410e-1,
        4.9374967810e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n81_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.1651518430e-4,    1.1729557070e-1,   -1.1390725560e-4,    1.3864241540e-2,
        -3.6578321670e-6,    1.4039471750e-2,    3.9979695430e-5,    2.1817624570e-2,
        -3.4771048380e-5,    1.9424453380e-2,   -4.4856792560e-6,    2.0071119070e-2,
        2.6764169890e-5,    2.4275571110e-2,   -2.1014409870e-5,    2.5223881010e-2,
        -4.5601450440e-6,    2.6712685820e-2,    2.0169904020e-5,    3.0749857430e-2,
        -1.4822886440e-5,    3.3638864760e-2,   -4.0566446840e-6,    3.6891996860e-2,
        1.5914813050e-5,    4.2760342360e-2,   -1.1573758460e-5,    4.9251139160e-2,
        -3.0833696200e-6,    5.7492613790e-2,    1.2908222740e-5,    7.0936918260e-2,
        -1.0191954060e-5,    9.1266930100e-2,   -1.6205276550e-6,    1.2699872260e-1,
        1.0912397560e-5,    2.1227568390e-1,   -1.0024234140e-5,    6.3694047930e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.1287263930e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.1666668840e-3,    8.3333337680e-3,    2.0833333950e-2,    3.3333335070e-2,
        4.5833334330e-2,    5.8333333580e-2,    7.0833340290e-2,    8.3333350720e-2,
        9.5833361150e-2,    1.0833337160e-1,    1.2083338200e-1,    1.3333337010e-1,
        1.4583335820e-1,    1.5833334620e-1,    1.7499999700e-1,    1.8749998510e-1,
        1.9999997320e-1,    2.1249996130e-1,    2.2499994930e-1,    2.3749993740e-1,
        2.4999992550e-1,    2.6249992850e-1,    2.7499991660e-1,    2.8749990460e-1,
        2.9999989270e-1,    3.1249988080e-1,    3.2499986890e-1,    3.3749985690e-1,
        3.5416650770e-1,    3.6666649580e-1,    3.7916648390e-1,    3.9166647200e-1,
        4.0416646000e-1,    4.1666644810e-1,    4.2916643620e-1,    4.4166642430e-1,
        4.5416641240e-1,    4.6666640040e-1,    4.7916638850e-1,    4.9166637660e-1,
        4.9583303930e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n81_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        1.8447291800e-3,    3.6917933820e-1,   -1.7381771470e-3,    1.9923329350e-2,
        2.0874722400e-8,    2.0643860100e-2,   -1.4512792400e-5,    2.1501362320e-2,
        -2.5970221030e-5,    2.2513478990e-2,    1.8798498790e-5,    2.3805081840e-2,
        -3.8626676540e-6,    2.5251060720e-2,    1.7921411200e-6,    2.7011632920e-2,
        6.9688758230e-6,    2.9117226600e-2,    1.8330676540e-5,    3.1630694870e-2,
        -2.4178239980e-5,    3.4854710100e-2,    5.8544108470e-5,    3.8519382480e-2,
        -1.3998505890e-5,    4.3339073660e-2,   -4.9445130570e-5,    4.9596667290e-2,
        -8.3552135040e-5,    5.7985365390e-2,    1.2627318210e-4,    7.1573078630e-2,
        -3.7350910130e-5,    9.1360628600e-2,    1.5197969330e-5,    1.2762224670e-1,
        1.0652860510e-5,    2.1238636970e-1,   -7.7329023040e-6,    6.3666319850e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.1756911280e-1,
    ];
    let expected_extremal_frequencies = vec![
        8.3333335350e-4,    8.3333337680e-3,    1.9999997690e-2,    3.2499987630e-2,
        4.4999975710e-2,    5.8333296330e-2,    7.0833288130e-2,    8.3333276210e-2,
        9.6666596830e-2,    1.0916658490e-1,    1.2166657300e-1,    1.3499990110e-1,
        1.4749988910e-1,    1.5999987720e-1,    1.7333319780e-1,    1.8583318590e-1,
        1.9833317400e-1,    2.1166649460e-1,    2.2416648270e-1,    2.3749980330e-1,
        2.4999979140e-1,    2.6249977950e-1,    2.7583310010e-1,    2.8833308820e-1,
        3.0083307620e-1,    3.1416639690e-1,    3.2666638490e-1,    3.3999970560e-1,
        3.5249969360e-1,    3.6499968170e-1,    3.7833300230e-1,    3.9083299040e-1,
        4.0333297850e-1,    4.1666629910e-1,    4.2916628720e-1,    4.4166627530e-1,
        4.5499959590e-1,    4.6749958400e-1,    4.7999957200e-1,    4.9166622760e-1,
        4.9916622040e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n81_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        1.5095686540e-2,    3.7563890220e-1,   -1.4260130000e-2,    2.0715206860e-2,
        -8.5543142630e-6,    2.0846247670e-2,    3.4908764060e-5,    2.1742075680e-2,
        4.3003063180e-5,    2.2778719660e-2,   -1.1777156030e-5,    2.3969620470e-2,
        -6.5409054510e-5,    2.5445938110e-2,   -3.4007825890e-5,    2.7091681960e-2,
        5.7925935830e-6,    2.9082179070e-2,   -1.2540788160e-4,    3.1565845010e-2,
        -8.4825733210e-5,    3.4662663940e-2,   -9.9439814220e-5,    3.8358449940e-2,
        -1.0395812570e-4,    4.3185830120e-2,    1.8009683120e-5,    4.9696385860e-2,
        2.2103718950e-4,    5.8632493020e-2,    1.4844053660e-4,    7.1078240870e-2,
        -3.7032898400e-4,    9.1595947740e-2,    2.1849118640e-6,    1.2754184010e-1,
        7.3762959800e-5,    2.1233177190e-1,   -3.4584081730e-6,    6.3665068150e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.3157209160e-1,
    ];
    let expected_extremal_frequencies = vec![
        7.8125001160e-4,    7.8125009310e-3,    2.0312497390e-2,    3.2812487330e-2,
        4.5312475410e-2,    5.7812463490e-2,    7.0312485100e-2,    8.3593785760e-2,
        9.6093833450e-2,    1.0859388110e-1,    1.2187518180e-1,    1.3437522950e-1,
        1.4687527720e-1,    1.6015657780e-1,    1.7265662550e-1,    1.8593792620e-1,
        1.9843797390e-1,    2.1093802150e-1,    2.2421932220e-1,    2.3671936990e-1,
        2.5000065570e-1,    2.6250046490e-1,    2.7500027420e-1,    2.8828132150e-1,
        3.0078113080e-1,    3.1406217810e-1,    3.2656198740e-1,    3.3906179670e-1,
        3.5234284400e-1,    3.6484265330e-1,    3.7812370060e-1,    3.9062350990e-1,
        4.0312331910e-1,    4.1640436650e-1,    4.2890417580e-1,    4.4218522310e-1,
        4.5468503240e-1,    4.6718484160e-1,    4.7968465090e-1,    4.9218446020e-1,
        4.9921560290e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n82_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        6.0107924040e-2,    1.0152831670e-2,    2.8129592540e-3,    5.4209083320e-3,
        8.7474584580e-3,    8.9562386270e-3,    7.4334070090e-3,    6.9568157200e-3,
        8.0109983680e-3,    9.0153291820e-3,    8.9699327950e-3,    8.5939019920e-3,
        8.9252293110e-3,    9.8366141320e-3,    1.0468825700e-2,    1.0582804680e-2,
        1.0785847900e-2,    1.1524021630e-2,    1.2475326660e-2,    1.3140723110e-2,
        1.3604670760e-2,    1.4344587920e-2,    1.5495479110e-2,    1.6709953550e-2,
        1.7748475070e-2,    1.8851637840e-2,    2.0370304580e-2,    2.2278189660e-2,
        2.4304687980e-2,    2.6442080740e-2,    2.9048383240e-2,    3.2438784840e-2,
        3.6637991670e-2,    4.1701734070e-2,    4.8171579840e-2,    5.7142615320e-2,
        7.0254385470e-2,    9.0674936770e-2,    1.2706995010e-1,    2.1192616220e-1,
        6.3648653030e-1,
    ];
    let expected_deviations = vec![
        1.0247828070e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.0975607480e-3,    1.2195121500e-2,    1.8292682250e-2,    3.0487803740e-2,
        4.2682923380e-2,    5.4878048600e-2,    6.7073173820e-2,    7.9268299040e-2,
        9.1463424270e-2,    1.0365854950e-1,    1.1585367470e-1,    1.2804879250e-1,
        1.4024390280e-1,    1.5243901310e-1,    1.6463412340e-1,    1.7682923380e-1,
        1.8902434410e-1,    2.0121945440e-1,    2.1341456470e-1,    2.3170723020e-1,
        2.4390234050e-1,    2.5609746580e-1,    2.6829257610e-1,    2.8048768640e-1,
        2.9268279670e-1,    3.0487790700e-1,    3.1707301740e-1,    3.2926812770e-1,
        3.4146323800e-1,    3.5365834830e-1,    3.6585345860e-1,    3.7804856900e-1,
        3.9024367930e-1,    4.0243878960e-1,    4.1463389990e-1,    4.2682901020e-1,
        4.3902412060e-1,    4.5121923090e-1,    4.6341434120e-1,    4.7560945150e-1,
        4.8780456190e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n82_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.1174333840e-1,    5.3022503850e-3,    5.3161382680e-3,    1.1402361090e-2,
        7.1522444490e-3,    7.2031170130e-3,    9.7044110300e-3,    8.5439682010e-3,
        8.6466372010e-3,    1.0198339820e-2,    9.8475068810e-3,    1.0016545650e-2,
        1.1200353500e-2,    1.1223495010e-2,    1.1480897670e-2,    1.2521654370e-2,
        1.2798339130e-2,    1.3176560400e-2,    1.4204114680e-2,    1.4716356990e-2,
        1.5268474820e-2,    1.6393929720e-2,    1.7186880110e-2,    1.8003553150e-2,
        1.9366979600e-2,    2.0562052730e-2,    2.1814852950e-2,    2.3651212450e-2,
        2.5518149140e-2,    2.7571678160e-2,    3.0378401280e-2,    3.3574998380e-2,
        3.7358105180e-2,    4.2489945890e-2,    4.9046754840e-2,    5.7792842390e-2,
        7.0766210560e-2,    9.1023862360e-2,    1.2724542620e-1,    2.1221697330e-1,
        6.3669800760e-1,
    ];
    let expected_deviations = vec![
        2.1031156180e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.0650404990e-3,    8.1300809980e-3,    2.0325202490e-2,    3.2520323990e-2,
        4.4715445490e-2,    5.6910566990e-2,    6.9105684760e-2,    8.1300795080e-2,
        9.3495905400e-2,    1.0569101570e-1,    1.1788612600e-1,    1.3008123640e-1,
        1.4227634670e-1,    1.5447145700e-1,    1.6666656730e-1,    1.7886167760e-1,
        1.9105678800e-1,    2.0325189830e-1,    2.1544700860e-1,    2.2764211890e-1,
        2.3983722930e-1,    2.5203233960e-1,    2.6422744990e-1,    2.7642256020e-1,
        2.8861767050e-1,    3.0081278090e-1,    3.1300789120e-1,    3.2520300150e-1,
        3.4146314860e-1,    3.5365825890e-1,    3.6585336920e-1,    3.7804847960e-1,
        3.9024358990e-1,    4.0243870020e-1,    4.1463381050e-1,    4.2682892080e-1,
        4.3902403120e-1,    4.5121914150e-1,    4.6341425180e-1,    4.7560936210e-1,
        4.8780447240e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n82_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        3.6254966260e-1,    9.5179677010e-3,    9.6667110920e-3,    9.8333954810e-3,
        1.0017365220e-2,    1.0215133430e-2,    1.0417819020e-2,    1.0649174450e-2,
        1.0873138900e-2,    1.1118590830e-2,    1.1380553250e-2,    1.1636704210e-2,
        1.1886119840e-2,    1.2128025290e-2,    1.2356162070e-2,    1.3700664040e-2,
        1.3668954370e-2,    1.4184713360e-2,    1.4754116540e-2,    1.5383243560e-2,
        1.6082525250e-2,    1.6845226290e-2,    1.7738163470e-2,    1.8669664860e-2,
        1.9757449630e-2,    2.1000385280e-2,    2.2403001790e-2,    2.3989737030e-2,
        2.5796473030e-2,    2.7872204780e-2,    3.0706048010e-2,    3.3781468870e-2,
        3.7694752220e-2,    4.2657136920e-2,    4.9156486990e-2,    5.8031201360e-2,
        7.0836424830e-2,    9.1097474100e-2,    1.2737953660e-1,    2.1222043040e-1,
        6.3662171360e-1,
    ];
    let expected_deviations = vec![
        7.1441710000e-1,
    ];
    let expected_extremal_frequencies = vec![
        8.1300811140e-4,    8.1300819290e-3,    1.9512193280e-2,    3.1707305460e-2,
        4.3902415780e-2,    5.6097526100e-2,    6.8292640150e-2,    8.0487750470e-2,
        9.2682860790e-2,    1.0487797110e-1,    1.1707308140e-1,    1.3008120660e-1,
        1.4227631690e-1,    1.5447142720e-1,    1.6666653750e-1,    1.7886164780e-1,
        1.9105675820e-1,    2.0406487580e-1,    2.1625998620e-1,    2.2845509650e-1,
        2.4065020680e-1,    2.5284531710e-1,    2.6585343480e-1,    2.7804854510e-1,
        2.9024365540e-1,    3.0243876580e-1,    3.1463387610e-1,    3.2682898640e-1,
        3.3983710410e-1,    3.5203221440e-1,    3.6422732470e-1,    3.7642243500e-1,
        3.8861754540e-1,    4.0162566300e-1,    4.1382077340e-1,    4.2601588370e-1,
        4.3821099400e-1,    4.5040610430e-1,    4.6260121460e-1,    4.7560933230e-1,
        4.8780444260e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n82_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        3.6993753910e-1,    9.5055103300e-3,    9.6253156660e-3,    9.8092854020e-3,
        9.9668502810e-3,    1.0126769540e-2,    1.0324776170e-2,    1.0551601650e-2,
        1.0788530110e-2,    1.1050641540e-2,    1.1359095570e-2,    1.1703580620e-2,
        1.2072771790e-2,    1.2500375510e-2,    1.3027966020e-2,    1.3607561590e-2,
        1.3238668440e-2,    1.4192581180e-2,    1.4794349670e-2,    1.5362024310e-2,
        1.6044318680e-2,    1.6829252240e-2,    1.7686367030e-2,    1.8632709980e-2,
        1.9712388520e-2,    2.0943164830e-2,    2.2334337230e-2,    2.3929774760e-2,
        2.5799572470e-2,    2.7986347680e-2,    3.0531406400e-2,    3.3590495590e-2,
        3.7723302840e-2,    4.2634963990e-2,    4.9053907390e-2,    5.8007657530e-2,
        7.0870757100e-2,    9.1038107870e-2,    1.2738108630e-1,    2.1224689480e-1,
        6.3663625720e-1,
    ];
    let expected_deviations = vec![
        7.2931760550e-1,
    ];
    let expected_extremal_frequencies = vec![
        7.6219509360e-4,    7.6219514010e-3,    1.9054876640e-2,    3.1249986960e-2,
        4.3445099150e-2,    5.5640209470e-2,    6.7835323510e-2,    8.0792628230e-2,
        9.2987738550e-2,    1.0518284890e-1,    1.1737795920e-1,    1.2957307700e-1,
        1.4176818730e-1,    1.5472549200e-1,    1.6692060230e-1,    1.7911571260e-1,
        1.9131082300e-1,    2.0350593330e-1,    2.1646323800e-1,    2.2865834830e-1,
        2.4085345860e-1,    2.5304856900e-1,    2.6524367930e-1,    2.7820098400e-1,
        2.9039609430e-1,    3.0259120460e-1,    3.1478631500e-1,    3.2698142530e-1,
        3.3917653560e-1,    3.5213384030e-1,    3.6432895060e-1,    3.7652406100e-1,
        3.8871917130e-1,    4.0091428160e-1,    4.1387158630e-1,    4.2606669660e-1,
        4.3826180700e-1,    4.5045691730e-1,    4.6265202760e-1,    4.7560933230e-1,
        4.8780444260e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,3,1,2
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n128_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        5.8194641020e-2,    6.5329335630e-3,   -3.1796842810e-4,    3.3583119510e-3,
        6.5126717090e-3,    5.7308003310e-3,    3.9188191290e-3,    4.0093883870e-3,
        5.2616298200e-3,    5.5730268360e-3,    4.8499479890e-3,    4.5580044390e-3,
        5.1535815000e-3,    5.6483075020e-3,    5.4382309320e-3,    5.1505863670e-3,
        5.4254829880e-3,    5.9118568900e-3,    5.9982985260e-3,    5.8270990850e-3,
        5.9468001130e-3,    6.3702464100e-3,    6.6339373590e-3,    6.6161155700e-3,
        6.6886991260e-3,    7.0471763610e-3,    7.4146687980e-3,    7.5558573010e-3,
        7.6598674060e-3,    7.9806447030e-3,    8.4120631220e-3,    8.7057501080e-3,
        8.9007765050e-3,    9.2308223250e-3,    9.7174048420e-3,    1.0159641500e-2,
        1.0496020320e-2,    1.0898292060e-2,    1.1465668680e-2,    1.2071460490e-2,
        1.2604653840e-2,    1.3164222240e-2,    1.3882517810e-2,    1.4709770680e-2,
        1.5529483560e-2,    1.6379505400e-2,    1.7395287750e-2,    1.8596112730e-2,
        1.9894927740e-2,    2.1298021080e-2,    2.2940337660e-2,    2.4911463260e-2,
        2.7198016640e-2,    2.9841721060e-2,    3.3042132850e-2,    3.7054777150e-2,
        4.2128026490e-2,    4.8682868480e-2,    5.7595074180e-2,    7.0514798160e-2,
        9.0807497500e-2,    1.2722438570e-1,    2.1211647990e-1,    6.3657808300e-1,
    ];
    let expected_deviations = vec![
        1.0202699150e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    7.8125000000e-3,    1.1718750000e-2,    1.9531250000e-2,
        2.7343750000e-2,    3.5156250000e-2,    4.2968750000e-2,    5.0781250000e-2,
        5.8593750000e-2,    6.6406250000e-2,    7.4218750000e-2,    8.2031250000e-2,
        8.9843750000e-2,    9.7656250000e-2,    1.0546875000e-1,    1.1328125000e-1,
        1.2109375000e-1,    1.2890625000e-1,    1.3671875000e-1,    1.4453125000e-1,
        1.5234375000e-1,    1.6015625000e-1,    1.6796875000e-1,    1.7578125000e-1,
        1.8359375000e-1,    1.9140625000e-1,    1.9921875000e-1,    2.0703125000e-1,
        2.1484375000e-1,    2.2265625000e-1,    2.3437500000e-1,    2.4218750000e-1,
        2.5000000000e-1,    2.5781250000e-1,    2.6562500000e-1,    2.7343750000e-1,
        2.8125000000e-1,    2.8906250000e-1,    2.9687500000e-1,    3.0468750000e-1,
        3.1250000000e-1,    3.2031250000e-1,    3.2812500000e-1,    3.3593750000e-1,
        3.4375000000e-1,    3.5156250000e-1,    3.5937500000e-1,    3.6718750000e-1,
        3.7500000000e-1,    3.8281250000e-1,    3.9062500000e-1,    3.9843750000e-1,
        4.0625000000e-1,    4.1406250000e-1,    4.2187500000e-1,    4.2968750000e-1,
        4.3750000000e-1,    4.4531250000e-1,    4.5312500000e-1,    4.6093750000e-1,
        4.6875000000e-1,    4.7656250000e-1,    4.8437500000e-1,    4.9218750000e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,3,1,3
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n128_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.1001493780e-1,    2.9294937850e-3,    2.9311478140e-3,    8.5781812670e-3,
        3.9225891230e-3,    3.9339512590e-3,    6.1115324500e-3,    4.6191066500e-3,
        4.6436637640e-3,    5.8883875610e-3,    5.2075237040e-3,    5.2462518220e-3,
        6.0981214050e-3,    5.7533681390e-3,    5.8095455170e-3,    6.4587742090e-3,
        6.2933266160e-3,    6.3692182300e-3,    6.9042146210e-3,    6.8513602020e-3,
        6.9515556100e-3,    7.4200332160e-3,    7.4497163300e-3,    7.5784176590e-3,
        8.0115348100e-3,    8.1095993520e-3,    8.2730948930e-3,    8.6939632890e-3,
        8.8560581210e-3,    9.0626776220e-3,    9.4913244250e-3,    9.7202658650e-3,
        9.9817514420e-3,    1.0438472030e-2,    1.0743737220e-2,    1.1076807980e-2,
        1.1585921050e-2,    1.1985749010e-2,    1.2415230270e-2,    1.3009011750e-2,
        1.3533890250e-2,    1.4098346230e-2,    1.4825105670e-2,    1.5526026490e-2,
        1.6289174560e-2,    1.7226904630e-2,    1.8194168810e-2,    1.9267737870e-2,
        2.0556211470e-2,    2.1960616110e-2,    2.3562729360e-2,    2.5482535360e-2,
        2.7688324450e-2,    3.0304908750e-2,    3.3520877360e-2,    3.7456393240e-2,
        4.2435050010e-2,    4.8981249330e-2,    5.7880342010e-2,    7.0732235910e-2,
        9.0951681140e-2,    1.2732648850e-1,    2.1220552920e-1,    6.3662159440e-1,
    ];
    let expected_deviations = vec![
        2.0926253500e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.6041667440e-3,    5.2083334890e-3,    1.3020833950e-2,    2.0833332090e-2,
        2.8645830230e-2,    3.6458332090e-2,    4.4270835820e-2,    5.2083339540e-2,
        5.9895843270e-2,    6.7708335820e-2,    7.5520828370e-2,    8.3333320920e-2,
        9.1145813470e-2,    9.8958306010e-2,    1.0677079860e-1,    1.1458329110e-1,
        1.2239578370e-1,    1.3020828370e-1,    1.3802079860e-1,    1.4583331350e-1,
        1.5364582840e-1,    1.6145834330e-1,    1.6927085820e-1,    1.7708337310e-1,
        1.8489588800e-1,    1.9270840290e-1,    2.0052091780e-1,    2.0833343270e-1,
        2.1614594760e-1,    2.2395846250e-1,    2.3177097740e-1,    2.3958349230e-1,
        2.4739600720e-1,    2.5520849230e-1,    2.6302096250e-1,    2.7083343270e-1,
        2.7864590290e-1,    2.8645837310e-1,    2.9427084330e-1,    3.0208331350e-1,
        3.0989578370e-1,    3.1770825390e-1,    3.2552072410e-1,    3.3593735100e-1,
        3.4374982120e-1,    3.5156229140e-1,    3.5937476160e-1,    3.6718723180e-1,
        3.7499970200e-1,    3.8281217220e-1,    3.9062464240e-1,    3.9843711260e-1,
        4.0624958280e-1,    4.1406205300e-1,    4.2187452320e-1,    4.2968699340e-1,
        4.3749946360e-1,    4.4531193380e-1,    4.5312440400e-1,    4.6093687420e-1,
        4.6874934430e-1,    4.7656181450e-1,    4.8437428470e-1,    4.9218675490e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,3,1,15
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n128_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        3.6031413080e-1,    6.1427950860e-3,    6.1936378480e-3,    6.2527954580e-3,
        6.3050091270e-3,    6.3574314120e-3,    6.4031183720e-3,    6.4651668070e-3,
        6.5172612670e-3,    6.5571665760e-3,    6.6070258620e-3,    6.6571831700e-3,
        6.6901147370e-3,    6.7234933380e-3,    6.7480504510e-3,    7.7559053900e-3,
        7.2883963580e-3,    7.4130296710e-3,    7.5614750390e-3,    7.7025592330e-3,
        7.8565776350e-3,    8.0086886880e-3,    8.1854462620e-3,    8.3678960800e-3,
        8.5351467130e-3,    8.7242126460e-3,    8.9225173000e-3,    9.1118216510e-3,
        9.2985630040e-3,    9.4763636590e-3,    9.9833011630e-3,    1.0137915610e-2,
        1.0431051250e-2,    1.0757982730e-2,    1.1101365090e-2,    1.1473000050e-2,
        1.1865258220e-2,    1.2299835680e-2,    1.2783944610e-2,    1.3268768790e-2,
        1.3816773890e-2,    1.4416694640e-2,    1.5059709550e-2,    1.5746891500e-2,
        1.6481339930e-2,    1.7458379270e-2,    1.8385827540e-2,    1.9468128680e-2,
        2.0704448220e-2,    2.2108376030e-2,    2.3724913600e-2,    2.5592386720e-2,
        2.7792930600e-2,    3.0451953410e-2,    3.3608973030e-2,    3.7540912630e-2,
        4.2535305020e-2,    4.9061417580e-2,    5.7945251460e-2,    7.0758223530e-2,
        9.1012954710e-2,    1.2736213210e-1,    2.1222400670e-1,    6.3662636280e-1,
    ];
    let expected_deviations = vec![
        7.1327370410e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.2083336050e-4,    5.2083334890e-3,    1.2500003910e-2,    1.9791668280e-2,
        2.7604160830e-2,    3.5416655240e-2,    4.3229147790e-2,    5.1562473180e-2,
        5.9374965730e-2,    6.7187458280e-2,    7.4999950830e-2,    8.2812443380e-2,
        9.0624935930e-2,    9.8437428470e-2,    1.0624992100e-1,    1.1406241360e-1,
        1.2187490610e-1,    1.3020829860e-1,    1.3802090290e-1,    1.4583350720e-1,
        1.5364611150e-1,    1.6145871580e-1,    1.6927132010e-1,    1.7708392440e-1,
        1.8489652870e-1,    1.9270913300e-1,    2.0104257760e-1,    2.0885518190e-1,
        2.1666778620e-1,    2.2448039050e-1,    2.3229299490e-1,    2.4010559920e-1,
        2.4791820350e-1,    2.5573062900e-1,    2.6354300980e-1,    2.7187621590e-1,
        2.7968859670e-1,    2.8750097750e-1,    2.9531335830e-1,    3.0312573910e-1,
        3.1093811990e-1,    3.1875050070e-1,    3.2656288150e-1,    3.3489608760e-1,
        3.4270846840e-1,    3.5052084920e-1,    3.5833323000e-1,    3.6614561080e-1,
        3.7395799160e-1,    3.8177037240e-1,    3.8958275320e-1,    3.9739513400e-1,
        4.0572834010e-1,    4.1354072090e-1,    4.2135310170e-1,    4.2916548250e-1,
        4.3697786330e-1,    4.4479024410e-1,    4.5260262490e-1,    4.6041500570e-1,
        4.6874821190e-1,    4.7656059270e-1,    4.8437297340e-1,    4.9218535420e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,3,1,16
0.000000000000,0.500000000000
1.000000000000
1.000000000000
*/
#[test]
fn hilbert_allpass_n128_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.5f32,
        desired_value: 1f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        3.6773476000e-1,    5.9408545490e-3,    6.0067772870e-3,    6.0697197910e-3,
        6.1337947850e-3,    6.2035322190e-3,    6.2801539900e-3,    6.3653290270e-3,
        6.4599514010e-3,    6.5634548660e-3,    6.6774785520e-3,    6.8022906780e-3,
        6.9389641280e-3,    7.0901513100e-3,    7.2611570360e-3,    7.4638724330e-3,
        6.7461729050e-3,    7.4307322500e-3,    7.5566172600e-3,    7.6947212220e-3,
        7.8459680080e-3,    8.0076456070e-3,    8.1791877750e-3,    8.3590149880e-3,
        8.5469484330e-3,    8.7435245510e-3,    8.9475512500e-3,    9.1595649720e-3,
        9.3785524370e-3,    9.6027255060e-3,    9.8259449010e-3,    1.0028660300e-2,
        1.0586798190e-2,    1.0762512680e-2,    1.1110842230e-2,    1.1484384540e-2,
        1.1881232260e-2,    1.2305259700e-2,    1.2760639190e-2,    1.3253331180e-2,
        1.3788640500e-2,    1.4372467990e-2,    1.5013575550e-2,    1.5720605850e-2,
        1.6505360600e-2,    1.7382740970e-2,    1.8376350400e-2,    1.9536793230e-2,
        2.0625591280e-2,    2.2110402580e-2,    2.3724079130e-2,    2.5591611860e-2,
        2.7788817880e-2,    3.0409872530e-2,    3.3589124680e-2,    3.7520766260e-2,
        4.2506217960e-2,    4.9029946330e-2,    5.7928919790e-2,    7.0786595340e-2,
        9.0992689130e-2,    1.2736451630e-1,    2.1222901340e-1,    6.3656520840e-1,
    ];
    let expected_deviations = vec![
        7.2840905190e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.8828125000e-4,    4.8828125000e-3,    1.2207031250e-2,    2.0019531250e-2,
        2.7832031250e-2,    3.5644531250e-2,    4.3457031250e-2,    5.1269531250e-2,
        5.9082031250e-2,    6.6894531250e-2,    7.4707031250e-2,    8.2519531250e-2,
        9.0820312500e-2,    9.8632812500e-2,    1.0644531250e-1,    1.1425781250e-1,
        1.2207031250e-1,    1.2988281250e-1,    1.3769531250e-1,    1.4550781250e-1,
        1.5380859380e-1,    1.6162109380e-1,    1.6943359380e-1,    1.7724609380e-1,
        1.8505859380e-1,    1.9287109380e-1,    2.0068359380e-1,    2.0849609380e-1,
        2.1679687500e-1,    2.2460937500e-1,    2.3242187500e-1,    2.4023437500e-1,
        2.4804687500e-1,    2.5585937500e-1,    2.6367187500e-1,    2.7148437500e-1,
        2.7978515620e-1,    2.8759765620e-1,    2.9541015620e-1,    3.0322265620e-1,
        3.1103515620e-1,    3.1884765620e-1,    3.2666015620e-1,    3.3447265620e-1,
        3.4277343750e-1,    3.5058593750e-1,    3.5839843750e-1,    3.6621093750e-1,
        3.7402343750e-1,    3.8183593750e-1,    3.8964843750e-1,    3.9746093750e-1,
        4.0576171880e-1,    4.1357421880e-1,    4.2138671880e-1,    4.2919921880e-1,
        4.3701171880e-1,    4.4482421880e-1,    4.5263671880e-1,    4.6044921880e-1,
        4.6875000000e-1,    4.7656250000e-1,    4.8437500000e-1,    4.9218750000e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n3_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        3.2491970060e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.1803400520e-1,    6.1803400520e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n3_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        3.2491970060e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.1803400520e-1,    6.1803400520e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n3_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        4.1582340000e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.2709091900e-1,    8.2709091900e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.3333335070e-2,    2.6666668060e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n3_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        4.1946035620e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.3633464570e-1,    8.3633464570e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.1250000000e-2,    2.6249998810e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n3_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        4.6017354730e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        9.1978645320e-1,    9.1978645320e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3888888990e-2,    2.5555557010e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n3_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        4.6104982500e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        9.2180049420e-1,    9.2180049420e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3513513840e-2,    2.5405406950e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n4_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        2.7968627210e-1,    8.6213052270e-3,
    ];
    let expected_deviations = vec![
        5.4212987420e-1,    5.4212987420e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n4_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        2.7968627210e-1,    8.6213052270e-3,
    ];
    let expected_deviations = vec![
        5.4212987420e-1,    5.4212987420e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.0000000150e-1,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n4_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        4.4356808070e-1,    1.3672947880e-2,
    ];
    let expected_deviations = vec![
        8.5979008670e-1,    8.5979008670e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.6666667540e-2,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n4_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        4.4746017460e-1,    1.3792902230e-2,
    ];
    let expected_deviations = vec![
        8.6733436580e-1,    8.6733436580e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.5625000000e-2,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n4_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        4.8297488690e-1,    1.4887660740e-2,
    ];
    let expected_deviations = vec![
        9.3617427350e-1,    9.3617427350e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.9444444960e-3,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n4_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        4.8380830880e-1,    1.4913380150e-2,
    ];
    let expected_deviations = vec![
        9.3778979780e-1,    9.3778979780e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.7567569200e-3,    2.0000000300e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n11_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        6.5579704940e-2,    1.2447761740e-1,    1.9722500440e-1,    1.8593418600e-1,
        1.2222883110e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.8832907080e-2,    1.8832907080e-2,
    ];
    let expected_extremal_frequencies = vec![
        5.0000000750e-2,    1.0000000150e-1,    2.5000000000e-1,    3.0000001190e-1,
        3.5000002380e-1,    4.0000003580e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n11_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.1826926470e-1,    1.6338232160e-1,    2.0857614280e-1,    1.9411462550e-1,
        1.2013638020e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        9.9259436130e-2,    9.9259436130e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.3333335070e-2,    1.0000000150e-1,    2.0000000300e-1,    2.6666668060e-1,
        3.6666667460e-1,    4.6666666870e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n11_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        4.7961896660e-1,    2.4701356890e-1,    1.4495068790e-1,    1.9099980590e-2,
        -3.3464193340e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.8144679070e-1,    6.8144679070e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.6666668280e-3,    5.3333338350e-2,    2.0000000300e-1,    2.6666662100e-1,
        3.5333320500e-1,    4.5333310960e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n11_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        4.8893797400e-1,    2.4941068890e-1,    1.4478868250e-1,    1.4573514460e-2,
        -3.6952078340e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.9777143000e-1,    6.9777143000e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.2500000930e-3,    5.6250005960e-2,    2.0000000300e-1,    2.6249995830e-1,
        3.5624986890e-1,    4.4999977950e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n11_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        5.8456736800e-1,    2.7536359430e-1,    1.3572871690e-1,   -2.4208098650e-2,
        -7.2503626350e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.4969979520e-1,    8.4969979520e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.7777778450e-3,    5.2777778360e-2,    2.0000000300e-1,    2.6388904450e-1,
        3.5555595160e-1,    4.5277842880e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n11_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        5.8696877960e-1,    2.7601659300e-1,    1.3514864440e-1,   -2.5365620850e-2,
        -7.3453366760e-2,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.5342735050e-1,    8.5342735050e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.7027027680e-3,    5.4054044190e-2,    2.0000000300e-1,    2.6216211920e-1,
        3.5675707460e-1,    4.5135203000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n12_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        5.2461303770e-2,    1.1366332320e-1,    1.7478089030e-1,    2.0129993560e-1,
        1.6176086660e-1,    6.2305808070e-2,
    ];
    let expected_deviations = vec![
        2.3467978460e-2,    2.3467978460e-2,
    ];
    let expected_extremal_frequencies = vec![
        4.1666667910e-2,    1.0000000150e-1,    2.0000000300e-1,    2.8333333130e-1,
        3.2499998810e-1,    4.0833330150e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n12_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.5224774180e-1,    1.1281500760e-1,    1.4860552550e-1,    1.6940003630e-1,
        1.3204771280e-1,    5.0585329530e-2,
    ];
    let expected_deviations = vec![
        2.0020069180e-1,    2.0020069180e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.7777777980e-2,    5.5555555970e-2,    1.0000000150e-1,    2.2777777910e-1,
        3.1111115220e-1,    4.2222231630e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n12_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        4.5847454670e-1,    1.9827559590e-1,    1.5342706440e-1,    8.1408679490e-2,
        2.4448752400e-2,    2.2763013840e-3,
    ];
    let expected_deviations = vec![
        7.0877844100e-1,    7.0877844100e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.5555556900e-3,    5.0000000750e-2,    1.0000000150e-1,    2.3888888960e-1,
        3.2222241160e-1,    4.1111153360e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n12_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        4.6760317680e-1,    1.9982251520e-1,    1.5326964860e-1,    7.9199135300e-2,
        2.1320164200e-2,    1.2555122380e-3,
    ];
    let expected_deviations = vec![
        7.2383093830e-1,    7.2383093830e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.2083334890e-3,    5.2083328370e-2,    1.0000000150e-1,    2.3645830150e-1,
        3.2500010730e-1,    4.1354194280e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n12_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        5.5126535890e-1,    2.2003805640e-1,    1.5223664050e-1,    5.6257247920e-2,
        -6.7815184590e-3,   -1.1573612690e-2,
    ];
    let expected_deviations = vec![
        8.6399680380e-1,    8.6399680380e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.3148148320e-3,    5.0925917920e-2,    1.0000000150e-1,    2.3703713720e-1,
        3.2268503310e-1,    4.1064766050e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n12_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        5.5322211980e-1,    2.2040301560e-1,    1.5205842260e-1,    5.5696845050e-2,
        -7.3081851010e-3,   -1.1827468870e-2,
    ];
    let expected_deviations = vec![
        8.6739879850e-1,    8.6739879850e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.2522523070e-3,    4.9549542370e-2,    1.0000000150e-1,    2.3828826840e-1,
        3.2387381790e-1,    4.1171160340e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n16_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        8.0360263590e-2,    1.1984303590e-2,    2.9905244710e-2,    7.6090022920e-2,
        1.5044592320e-1,    2.1291613580e-1,    1.8628555540e-1,    8.0407440660e-2,
    ];
    let expected_deviations = vec![
        1.3119782510e-1,    1.3119782510e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.1250000000e-2,    6.2500000000e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.3125000300e-1,    2.9374998810e-1,    3.5624998810e-1,    4.1874998810e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n16_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.3390058280e-1,    3.4993618730e-2,    3.6687225100e-2,    7.3500871660e-2,
        1.5241077540e-1,    2.2151944040e-1,    2.0584303140e-1,    8.5885643960e-2,
    ];
    let expected_deviations = vec![
        2.2588364780e-1,    2.2588364780e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.0833333950e-2,    6.2500000000e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.4166665970e-1,    3.0416667460e-1,    3.6666670440e-1,    4.2916673420e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n16_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        4.2736727000e-1,    8.8921904560e-2,    4.9904584880e-3,   -9.9155902860e-3,
        1.0636991260e-1,    2.6884573700e-1,    3.0464118720e-1,    1.3453757760e-1,
    ];
    let expected_deviations = vec![
        7.2195541860e-1,    7.2195541860e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.1666668840e-3,    4.1666667910e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.4166662990e-1,    3.0416658520e-1,    3.7083318830e-1,    4.3333312870e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n16_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        4.3590971830e-1,    9.1254681350e-2,    3.8120746610e-3,   -1.2384057040e-2,
        1.0608369110e-1,    2.7017492060e-1,    3.0743426080e-1,    1.3628351690e-1,
    ];
    let expected_deviations = vec![
        7.3581862450e-1,    7.3581862450e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    4.2968750000e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.4296875300e-1,    3.0546873810e-1,    3.6796873810e-1,    4.3437498810e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n16_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        5.1450759170e-1,    1.0462325810e-1,   -5.0584673880e-3,   -3.6084353920e-2,
        9.2575132850e-2,    2.8246521950e-1,    3.3376389740e-1,    1.4967966080e-1,
    ];
    let expected_deviations = vec![
        8.7020558120e-1,    8.7020558120e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7361111240e-3,    3.9930567150e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.4166651070e-1,    3.0416628720e-1,    3.6840215330e-1,    4.3437412380e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n16_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        5.1650327440e-1,    1.0498332980e-1,   -5.3398013110e-3,   -3.6635398860e-2,
        9.2344224450e-2,    2.8276026250e-1,    3.3441799880e-1,    1.5011477470e-1,
    ];
    let expected_deviations = vec![
        8.7340235710e-1,    8.7340235710e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.6891892300e-3,    4.0540542450e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.4054069820e-1,    3.0304092170e-1,    3.6891955140e-1,    4.3479818110e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n47_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        1.4807622130e-1,    2.3760373890e-1,    1.3555328550e-1,   -2.5060144070e-1,
        -6.5300357340e-1,   -6.0876226430e-1,    1.6149166230e-1,     1.2267510890e0,
        1.6050263640e0,    5.7407104970e-1,    -1.4116206170e0,    -2.7600374220e0,
        -1.9513106350e0,    8.9317005870e-1,     3.6750602720e0,     3.8387289050e0,
        7.1689677240e-1,    -3.4870767590e0,    -5.2419157030e0,    -2.7179477210e0,
        2.3740015030e0,     5.9481649400e0,     4.9396743770e0,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.0427782830e-2,    7.0427782830e-2,
    ];
    let expected_extremal_frequencies = vec![
        1.0869565420e-2,    2.1739130840e-2,    3.2608695330e-2,    5.4347828030e-2,
        6.5217390660e-2,    7.6086953280e-2,    8.6956515910e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.1086956560e-1,    2.2173912820e-1,    2.4347825350e-1,
        2.6521739360e-1,    2.8695651890e-1,    2.9782608150e-1,    3.1956520680e-1,
        3.4130433200e-1,    3.6304345730e-1,    3.8478258250e-1,    4.0652170780e-1,
        4.2826083300e-1,    4.4999995830e-1,    4.7173908350e-1,    4.8260864620e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n47_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        3.0028298500e-1,    4.3983936310e-1,    1.9304981830e-1,   -5.4296791550e-1,
        -1.2523126600e0,    -1.0367308860e0,    4.5611816640e-1,     2.3368856910e0,
        2.8230276110e0,    7.7214622500e-1,    -2.8029618260e0,    -4.9845170970e0,
        -3.2542247770e0,     1.9011566640e0,     6.6421670910e0,     6.6193447110e0,
        9.3396306040e-1,    -6.3899607660e0,    -9.2623653410e0,    -4.7004313470e0,
        4.1809134480e0,     1.0325968740e1,     8.5387115480e0,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.7577099800e-1,    1.7577099800e-1,
    ];
    let expected_extremal_frequencies = vec![
        7.2463769470e-3,    1.4492753890e-2,    2.8985507790e-2,    5.0724633040e-2,
        6.5217383210e-2,    7.9710133370e-2,    8.6956508460e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0724637810e-1,    2.2898550330e-1,    2.4347825350e-1,
        2.6521739360e-1,    2.7971014380e-1,    3.0144926910e-1,    3.2318839430e-1,
        3.4492751960e-1,    3.6666664480e-1,    3.8115939500e-1,    4.0289852020e-1,
        4.2463764550e-1,    4.4637677070e-1,    4.6811589600e-1,    4.8985502120e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n47_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        8.8380569220e-1,     1.0206860300e0,    2.2673147920e-1,    -1.5110795500e0,
        -2.7674541470e0,    -1.7893837690e0,     1.6013460160e0,     5.0320663450e0,
        5.0494284630e0,    3.3048701290e-1,    -6.2306628230e0,    -9.0806312560e0,
        -4.6348953250e0,     4.8602442740e0,     1.2140161510e1,     1.0544713970e1,
        3.9982795720e-2,    -1.1797063830e1,    -1.5358884810e1,    -6.9586367610e0,
        7.5695443150e0,     1.7036701200e1,     1.3779428480e1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.8616795540e-1,    6.8616795540e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.4492754130e-3,    1.3043479990e-2,    3.0434789140e-2,    4.9275372180e-2,
        6.6666670140e-2,    8.2608662550e-2,    9.5652110870e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0724636320e-1,    2.2318835560e-1,    2.4057962000e-1,
        2.5942024590e-1,    2.7971026300e-1,    3.0000028010e-1,    3.2029029730e-1,
        3.4202960130e-1,    3.6231961850e-1,    3.8405892250e-1,    4.0434893970e-1,
        4.2608824370e-1,    4.4637826090e-1,    4.6811756490e-1,    4.8985686900e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n47_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        9.0260970590e-1,     1.0403857230e0,    2.2885894780e-1,    -1.5435528760e0,
        -2.8209612370e0,    -1.8177230360e0,     1.6398085360e0,     5.1276326180e0,
        5.1318216320e0,    3.1839513780e-1,    -6.3531951900e0,    -9.2325696950e0,
        -4.6919651030e0,     4.9635047910e0,     1.2343834880e1,     1.0697360040e1,
        1.3561248780e-2,    -1.1997344970e1,    -1.5593029020e1,    -7.0507369040e0,
        7.6945219040e0,     1.7293813710e1,     1.3981916430e1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.0230239630e-1,    7.0230239630e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3586956770e-3,    1.2228259820e-2,    2.9891299080e-2,    4.8913035540e-2,
        6.6576078530e-2,    8.2880422470e-2,    9.5108680430e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0815221970e-1,    2.2309795020e-1,    2.4076108630e-1,
        2.5978282090e-1,    2.8016313910e-1,    3.0054345730e-1,    3.2092377540e-1,
        3.4130409360e-1,    3.6304309960e-1,    3.8342341780e-1,    4.0516242380e-1,
        4.2554274200e-1,    4.4728174810e-1,    4.6766206620e-1,    4.8940107230e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n47_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        1.0888264180e0,     1.2465780970e0,    2.6062440870e-1,    -1.8738813400e0,
        -3.3978459840e0,    -2.1703271870e0,     1.9979270700e0,     6.1757764820e0,
        6.1454038620e0,    3.3323073390e-1,    -7.6744980810e0,    -1.1090588570e1,
        -5.5932030680e0,     6.0068893430e0,     1.4827991490e1,     1.2799268720e1,
        -4.3273925780e-2,    -1.4432884220e1,    -1.8708709720e1,    -8.4486389160e0,
        9.2160100940e0,     2.0701820370e1,     1.6732717510e1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.5187298060e-1,    8.5187298060e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.0386472610e-4,    1.2077296150e-2,    3.0193220820e-2,    4.8913057890e-2,
        6.7029006780e-2,    8.2729421560e-2,    9.5410525800e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0724643770e-1,    2.2294704620e-1,    2.4106313290e-1,
        2.5978285070e-1,    2.7971005440e-1,    3.0024111270e-1,    3.2077217100e-1,
        3.4190708400e-1,    3.6243814230e-1,    3.8357305530e-1,    4.0470796820e-1,
        4.2584288120e-1,    4.4697779420e-1,    4.6811270710e-1,    4.8924762010e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n47_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        1.0930318830e0,     1.2505500320e0,    2.6117587090e-1,    -1.8803453450e0,
        -3.4084970950e0,    -2.1760344510e0,     2.0053458210e0,     6.1947383880e0,
        6.1619954110e0,    3.3134579660e-1,    -7.6983470920e0,    -1.1120964050e1,
        -5.6052808760e0,     6.0268702510e0,     1.4868568420e1,     1.2830347060e1,
        -4.7669410710e-2,    -1.4472691540e1,    -1.8755901340e1,    -8.4676513670e0,
        9.2409238820e0,     2.0753593440e1,     1.6773622510e1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        8.5544675590e-1,    8.5544675590e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.8754405470e-4,    1.2338426900e-2,    3.0552279200e-2,    4.8766180870e-2,
        6.6980086270e-2,    8.2843810320e-2,    9.5182262360e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0763799550e-1,    2.2291398050e-1,    2.4054011700e-1,
        2.5992912050e-1,    2.7990591530e-1,    3.0047026280e-1,    3.2103461030e-1,
        3.4159895780e-1,    3.6275085810e-1,    3.8390275840e-1,    4.0446710590e-1,
        4.2561900620e-1,    4.4677090640e-1,    4.6792280670e-1,    4.8966225980e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n48_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        1.5183085200e-1,    2.6549273730e-1,    1.7645138500e-1,   -2.2566068170e-1,
        -7.1244859700e-1,   -7.7215504650e-1,   -3.2857418060e-2,     1.1993149520e0,
        1.9051271680e0,     1.1192537550e0,    -1.0228705410e0,    -2.9873127940e0,
        -2.8794193270e0,   -1.7520213130e-1,     3.3858325480e0,     4.8587946890e0,
        2.5373926160e0,    -2.2487230300e0,    -5.8033051490e0,    -4.9764895440e0,
        1.9808769230e-2,     5.4504685400e0,     6.9433751110e0,     3.1242008210e0,
    ];
    let expected_deviations = vec![
        7.3871612550e-2,    7.3871612550e-2,
    ];
    let expected_extremal_frequencies = vec![
        1.0416666980e-2,    2.0833333950e-2,    3.1250000000e-2,    5.2083335820e-2,
        6.2500000000e-2,    7.2916664180e-2,    8.3333328370e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.1041667460e-1,    2.2083334620e-1,    2.4166668950e-1,
        2.6250001790e-1,    2.7291667460e-1,    2.9374998810e-1,    3.1458330150e-1,
        3.3541661500e-1,    3.5624992850e-1,    3.7708324190e-1,    3.9791655540e-1,
        4.1874986890e-1,    4.3958318230e-1,    4.6041649580e-1,    4.8124980930e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n48_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        2.6712989810e-1,    3.3368164300e-1,    1.1523860690e-1,   -4.0328830480e-1,
        -8.2784640790e-1,   -6.1193829770e-1,    3.5387796160e-1,     1.4190604690e0,
        1.5726280210e0,    3.5964655880e-1,    -1.5077997450e0,    -2.5193133350e0,
        -1.5971138480e0,    8.3945131300e-1,     2.9988274570e0,     3.0857691760e0,
        8.0811691280e-1,    -2.2540843490e0,    -3.7970490460e0,    -2.5616269110e0,
        6.6340780260e-1,     3.5919954780e0,     4.1022906300e0,     1.7799739840e0,
    ];
    let expected_deviations = vec![
        1.8476507070e-1,    1.8476507070e-1,
    ];
    let expected_extremal_frequencies = vec![
        6.9444444960e-3,    1.3888888990e-2,    3.4722223880e-2,    4.8611111940e-2,
        6.9444447760e-2,    8.3333343270e-2,    9.0277791020e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0694445070e-1,    2.2083334620e-1,    2.3472224180e-1,
        2.5555557010e-1,    2.7638891340e-1,    2.9722225670e-1,    3.1805559990e-1,
        3.3194449540e-1,    3.5277783870e-1,    3.7361118200e-1,    3.9444452520e-1,
        4.1527786850e-1,    4.3611121180e-1,    4.5694455500e-1,    4.7777789830e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n48_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        8.5338371990e-1,    9.1915446520e-1,    1.2674307820e-1,    -1.4230029580e0,
        -2.3536796570e0,    -1.2498666050e0,     1.6695464850e0,     4.1183514600e0,
        3.4873654840e0,   -5.1690626140e-1,    -4.9727697370e0,    -5.9310855870e0,
        -2.0033054350e0,     4.0245637890e0,     7.3207063670e0,     5.1222352980e0,
        -8.9891815190e-1,    -6.1208920480e0,    -6.8139367100e0,    -2.8730335240e0,
        2.6354298590e0,     6.1362800600e0,     5.8651180270e0,     2.3629736900e0,
    ];
    let expected_deviations = vec![
        6.9382727150e-1,    6.9382727150e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3888889230e-3,    1.2500000190e-2,    3.0555555600e-2,    4.8611111940e-2,
        6.6666677590e-2,    8.1944495440e-2,    9.4444528220e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0555557310e-1,    2.1944449840e-1,    2.3611120880e-1,
        2.5555562970e-1,    2.7361103890e-1,    2.9444420340e-1,    3.1388849020e-1,
        3.3472165470e-1,    3.5555481910e-1,    3.7499910590e-1,    3.9583227040e-1,
        4.1666543480e-1,    4.3749859930e-1,    4.5833176370e-1,    4.7916492820e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n48_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        8.7121266130e-1,    9.3619972470e-1,    1.2678802010e-1,    -1.4548695090e0,
        -2.3992280960e0,    -1.2660019400e0,     1.7143391370e0,     4.2013311390e0,
        3.5370473860e0,   -5.5574083330e-1,    -5.0827922820e0,    -6.0231080060e0,
        -1.9927144050e0,     4.1362061500e0,     7.4410476680e0,     5.1511282920e0,
        -9.7975826260e-1,    -6.2344284060e0,    -6.8667812350e0,    -2.8394398690e0,
        2.7046289440e0,     6.1783428190e0,     5.8666057590e0,     2.3562917710e0,
    ];
    let expected_deviations = vec![
        7.0897036790e-1,    7.0897036790e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.3020833720e-3,    1.1718749070e-2,    2.9947921630e-2,    4.8177070920e-2,
        6.6406227650e-2,    8.2031257450e-2,    9.5052115620e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0520831640e-1,    2.1953117850e-1,    2.3645819720e-1,
        2.5468733910e-1,    2.7421873810e-1,    2.9375013710e-1,    3.1458362940e-1,
        3.3411502840e-1,    3.5494852070e-1,    3.7578201290e-1,    3.9661550520e-1,
        4.1744899750e-1,    4.3828248980e-1,    4.5781388880e-1,    4.7864738110e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,2,36
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n48_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        1.0459802150e0,     1.1168603900e0,    1.3796401020e-1,    -1.7583754060e0,
        -2.8772239690e0,    -1.4989068510e0,     2.0798683170e0,     5.0384650230e0,
        4.2036771770e0,   -7.2335004810e-1,    -6.1233682630e0,    -7.1894860270e0,
        -2.3152375220e0,     5.0076427460e0,     8.8859815600e0,     6.0670390130e0,
        -1.2728881840e0,    -7.4760875700e0,    -8.1327371600e0,    -3.2970495220e0,
        3.2414307590e0,     7.2684688570e0,     6.8515753750e0,     2.7419338230e0,
    ];
    let expected_deviations = vec![
        8.5573506360e-1,    8.5573506360e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.7870370800e-4,    1.2152776120e-2,    3.0092580240e-2,    4.8611141740e-2,
        6.5972276030e-2,    8.2175917920e-2,    9.4907350840e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0578701790e-1,    2.1909715240e-1,    2.3645819720e-1,
        2.5497666000e-1,    2.7407380940e-1,    2.9432836170e-1,    3.1400421260e-1,
        3.3483746650e-1,    3.5509201880e-1,    3.7592527270e-1,    3.9617982510e-1,
        4.1701307890e-1,    4.3784633280e-1,    4.5867958660e-1,    4.7951284050e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,2,37
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n48_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        1.0503340960e0,     1.1217209100e0,    1.3851046560e-1,    -1.7663915160e0,
        -2.8906395440e0,    -1.5062100890e0,     2.0887842180e0,     5.0617198940e0,
        4.2243595120e0,   -7.2496461870e-1,    -6.1514239310e0,    -7.2252731320e0,
        -2.3298320770e0,     5.0286674500e0,     8.9295091630e0,     6.1009459500e0,
        -1.2742776870e0,    -7.5120358470e0,    -8.1778059010e0,    -3.3205285070e0,
        3.2543196680e0,     7.3077206610e0,     6.8919553760e0,     2.7587585450e0,
    ];
    let expected_deviations = vec![
        8.5932642220e-1,    8.5932642220e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.6306307670e-4,    1.1824322860e-2,    2.9842330140e-2,    4.8423402010e-2,
        6.6441409290e-2,    8.2207165660e-2,    9.5157608390e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0563070480e-1,    2.1914438900e-1,    2.3603649440e-1,
        2.5461769100e-1,    2.7432462570e-1,    2.9403156040e-1,    3.1430155040e-1,
        3.3457154040e-1,    3.5484153030e-1,    3.7567457560e-1,    3.9650762080e-1,
        4.1677761080e-1,    4.3761065600e-1,    4.5844370130e-1,    4.7927674650e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n81_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        5.3812038900e-1,     1.7131571770e0,     2.0364384650e0,    -1.4967352150e0,
        -9.4397840500e0,    -1.4615024570e1,    -3.9776778220e0,     2.6788675310e1,
        5.5906501770e1,     4.0532848360e1,    -4.2046283720e1,    -1.4485314940e2,
        -1.5604710390e2,     6.7773132320e0,     2.7303875730e2,     4.0435076900e2,
        1.7375463870e2,    -3.6618215940e2,    -7.9534765620e2,    -6.0994940190e2,
        2.6963580320e2,     1.2339326170e3,     1.3503405760e3,     2.0579449460e2,
        -1.4968813480e3,    -2.2908527830e3,    -1.1719951170e3,     1.2970466310e3,
        3.1423632810e3,     2.5408474120e3,    -4.3037658690e2,    -3.5111557620e3,
        -3.9781865230e3,    -1.0654326170e3,     3.0827365720e3,     5.0001103520e3,
        2.8518117680e3,    -1.8118100590e3,    -5.1845683590e3,    -4.3916928710e3,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.3114970920e-2,    6.3114970920e-2,
    ];
    let expected_extremal_frequencies = vec![
        6.2500000930e-3,    1.2500000190e-2,    1.8750000750e-2,    2.5000000370e-2,
        3.7500001490e-2,    5.0000004470e-2,    5.6250005960e-2,    6.8750008940e-2,
        7.5000010430e-2,    8.1250011920e-2,    8.7500013410e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0624999700e-1,    2.1249999110e-1,    2.1874998510e-1,
        2.2499997910e-1,    2.3124997320e-1,    2.4374996130e-1,    2.4999995530e-1,
        2.6249995830e-1,    2.6874995230e-1,    2.8124994040e-1,    2.9374992850e-1,
        3.0624991660e-1,    3.1874990460e-1,    3.2499989870e-1,    3.3749988680e-1,
        3.4999987480e-1,    3.6249986290e-1,    3.7499985100e-1,    3.8749983910e-1,
        3.9999982710e-1,    4.1249981520e-1,    4.2499980330e-1,    4.3124979730e-1,
        4.4374978540e-1,    4.5624977350e-1,    4.6874976160e-1,    4.8124974970e-1,
        4.9374973770e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n81_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        8.4406143430e-1,     2.2407298090e0,     2.0627751350e0,    -2.7689542770e0,
        -1.1324498180e1,    -1.4543071750e1,   -1.9024276730e-1,     3.1283872600e1,
        5.4064376830e1,     2.9452615740e1,    -5.2558189390e1,    -1.3698240660e2,
        -1.2494512940e2,     3.5224952700e1,     2.5601290890e2,     3.2948535160e2,
        9.6599266050e1,    -3.5074563600e2,    -6.4873321530e2,    -4.2827340700e2,
        2.9941625980e2,     1.0057081300e3,     9.9429132080e2,     4.3067749020e1,
        -1.2266922610e3,    -1.7111972660e3,    -7.5929003910e2,     1.0933176270e3,
        2.3571972660e3,     1.7792647710e3,    -4.5375152590e2,    -2.6356484380e3,
        -2.8495029300e3,    -6.6053039550e2,     2.3131240230e3,     3.6077622070e3,
        1.9940588380e3,    -1.3587404790e3,    -3.7403984380e3,    -3.1450397950e3,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.7313304540e-1,    1.7313304540e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.1666668840e-3,    8.3333337680e-3,    1.6666667540e-2,    2.9166666790e-2,
        3.7500001490e-2,    5.0000000750e-2,    5.8333333580e-2,    6.6666670140e-2,
        7.5000010430e-2,    8.3333350720e-2,    8.7500020860e-2,    9.1666691010e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0416666570e-1,    2.1249999110e-1,
        2.1666665380e-1,    2.2916664180e-1,    2.3749996720e-1,    2.4583329260e-1,
        2.5833329560e-1,    2.7083328370e-1,    2.7916660900e-1,    2.9166659710e-1,
        3.0416658520e-1,    3.1249991060e-1,    3.2499989870e-1,    3.3749988680e-1,
        3.4999987480e-1,    3.6249986290e-1,    3.7499985100e-1,    3.8749983910e-1,
        3.9583316450e-1,    4.0833315250e-1,    4.2083314060e-1,    4.3333312870e-1,
        4.4583311680e-1,    4.5833310480e-1,    4.7083309290e-1,    4.8333308100e-1,
        4.9583306910e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n81_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        2.0521576400e0,     4.3678870200e0,     2.7348687650e0,    -6.8547601700e0,
        -1.9628036500e1,    -1.9537321090e1,     7.7271423340e0,     5.1846374510e1,
        7.0377655030e1,     1.9368749620e1,    -9.1393142700e1,    -1.7364622500e2,
        -1.1666616820e2,     9.7160560610e1,     3.2186682130e2,     3.2677355960e2,
        6.0978240970e0,    -4.5464678960e2,    -6.5323376460e2,    -2.9935235600e2,
        4.5477093510e2,     1.0202459720e3,     8.1229766850e2,    -1.8712844850e2,
        -1.2642575680e3,    -1.4646165770e3,    -4.2205230710e2,     1.1855599370e3,
        2.0527031250e3,     1.3070603030e3,    -6.4892193600e2,    -2.3104409180e3,
        -2.2397602540e3,    -3.1513012700e2,     2.0325345460e3,     2.8993398440e3,
        1.4785828860e3,    -1.1947462160e3,    -3.0109350590e3,    -2.4865605470e3,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        6.8394327160e-1,    6.8394327160e-1,
    ];
    let expected_extremal_frequencies = vec![
        8.3333335350e-4,    7.5000007640e-3,    1.7500000070e-2,    2.8333323080e-2,
        3.9166647940e-2,    4.9999970940e-2,    5.9999961410e-2,    6.9999955590e-2,
        7.9166613520e-2,    8.7499938910e-2,    9.4166599210e-2,    9.8333261910e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0250000060e-1,    2.0749999580e-1,
        2.1499998870e-1,    2.2416664660e-1,    2.3416663710e-1,    2.4499996010e-1,
        2.5583329800e-1,    2.6666662100e-1,    2.7833327650e-1,    2.8999993210e-1,
        3.0166658760e-1,    3.1333324310e-1,    3.2499989870e-1,    3.3666655420e-1,
        3.4833320980e-1,    3.6083319780e-1,    3.7249985340e-1,    3.8499984150e-1,
        3.9666649700e-1,    4.0916648510e-1,    4.2083314060e-1,    4.3333312870e-1,
        4.4499978420e-1,    4.5749977230e-1,    4.6999976040e-1,    4.8166641590e-1,
        4.9416640400e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n81_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        2.1005620960e0,     4.4707083700e0,     2.8001027110e0,    -7.0152263640e0,
        -2.0094150540e1,    -2.0011032100e1,     7.8907318120e0,     5.3068450930e1,
        7.2082633970e1,     1.9906538010e1,    -9.3507621770e1,    -1.7783694460e2,
        -1.1964769740e2,     9.9266372680e1,     3.2957455440e2,     3.3496505740e2,
        6.7632904050e0,    -4.6537847900e2,    -6.6944824220e2,    -3.0753173830e2,
        4.6511685180e2,     1.0453686520e3,     8.3341082760e2,    -1.9036306760e2,
        -1.2951035160e3,    -1.5020334470e3,    -4.3442211910e2,     1.2140084230e3,
        2.1046967770e3,     1.3418140870e3,    -6.6356115720e2,    -2.3686835940e3,
        -2.2979624020e3,    -3.2489013670e2,     2.0836269530e3,     2.9740561520e3,
        1.5176052250e3,    -1.2247395020e3,    -3.0884223630e3,    -2.5508835450e3,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        7.0035487410e-1,    7.0035487410e-1,
    ];
    let expected_extremal_frequencies = vec![
        7.8125001160e-4,    7.0312507450e-3,    1.7187500370e-2,    2.8124989940e-2,
        3.9062481370e-2,    4.9999970940e-2,    6.0156211260e-2,    7.0312485100e-2,
        7.8906267880e-2,    8.7500050660e-2,    9.3750074510e-2,    9.8437592390e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0234376190e-1,    2.0781253280e-1,
        2.1562506260e-1,    2.2421884540e-1,    2.3437513410e-1,    2.4453142290e-1,
        2.5546884540e-1,    2.6640617850e-1,    2.7812474970e-1,    2.8984332080e-1,
        3.0156189200e-1,    3.1328046320e-1,    3.2499903440e-1,    3.3671760560e-1,
        3.4843617680e-1,    3.6093598600e-1,    3.7265455720e-1,    3.8515436650e-1,
        3.9687293770e-1,    4.0859150890e-1,    4.2109131810e-1,    4.3280988930e-1,
        4.4530969860e-1,    4.5780950780e-1,    4.6952807900e-1,    4.8202788830e-1,
        4.9374645950e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n82_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        4.1855442520e-1,     1.1507024770e0,     1.1840397120e0,    -1.1179261210e0,
        -5.6077165600e0,    -8.0013313290e0,    -1.7690629960e0,     1.4269761090e1,
        2.8474048610e1,     2.0469167710e1,    -1.9141765590e1,    -6.7751853940e1,
        -7.5045135500e1,    -4.6167907710e0,     1.1437222290e2,     1.8146156310e2,
        9.7499389650e1,    -1.2638919070e2,    -3.2759130860e2,    -2.9471521000e2,
        3.4713500980e1,     4.5265585330e2,     5.9002606200e2,     2.2575671390e2,
        -4.5167807010e2,    -9.0523962400e2,    -6.6381616210e2,     2.1989208980e2,
        1.0971896970e3,     1.1916436770e3,     2.7893457030e2,    -1.0140682370e3,
        -1.6316229250e3,    -9.5933959960e2,     5.8209765620e2,     1.7860870360e3,
        1.6224272460e3,     1.3153930660e2,    -1.5371789550e3,    -2.0327092290e3,
        -9.2287182620e2,
    ];
    let expected_deviations = vec![
        7.1910955010e-2,    7.1910955010e-2,
    ];
    let expected_extremal_frequencies = vec![
        6.0975607480e-3,    1.2195121500e-2,    1.8292682250e-2,    3.0487803740e-2,
        3.6585364490e-2,    4.8780485990e-2,    5.4878048600e-2,    6.7073173820e-2,
        7.3170736430e-2,    7.9268299040e-2,    8.5365861650e-2,    9.1463424270e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0609755810e-1,    2.1219511330e-1,
        2.1829266850e-1,    2.2439022360e-1,    2.3658533390e-1,    2.4268288910e-1,
        2.5487801430e-1,    2.6707312460e-1,    2.7926823500e-1,    2.8536579010e-1,
        2.9756090040e-1,    3.0975601080e-1,    3.2195112110e-1,    3.3414623140e-1,
        3.4634134170e-1,    3.5853645210e-1,    3.7073156240e-1,    3.8292667270e-1,
        3.9512178300e-1,    4.0731689330e-1,    4.1951200370e-1,    4.2560955880e-1,
        4.3780466910e-1,    4.4999977950e-1,    4.6219488980e-1,    4.7439000010e-1,
        4.8658511040e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n82_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        6.7489093540e-1,     1.6062247750e0,     1.2842133050e0,    -2.0696079730e0,
        -7.3308591840e0,    -8.6291418080e0,    6.6230869290e-1,     1.8783056260e1,
        3.0227376940e1,     1.4898323060e1,    -2.9489715580e1,    -7.1821357730e1,
        -6.2914150240e1,     1.8242431640e1,     1.2480062870e2,     1.5826181030e2,
        4.9591918950e1,    -1.5468859860e2,    -2.9220202640e2,    -2.0432676700e2,
        1.0438119510e2,     4.1587768550e2,     4.4250079350e2,     8.0332641600e1,
        -4.4427642820e2,    -7.0251123050e2,    -4.0862695310e2,     2.9224328610e2,
        8.7076422120e2,     8.1294342040e2,     7.1320434570e1,    -8.2649328610e2,
        -1.1545986330e3,    -5.8154821780e2,     5.1023333740e2,     1.2789974370e3,
        1.0844984130e3,     2.7166503910e1,    -1.0936152340e3,    -1.3971123050e3,
        -6.2764294430e2,
    ];
    let expected_deviations = vec![
        1.7889048160e-1,    1.7889048160e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.0650404990e-3,    8.1300809980e-3,    1.6260162000e-2,    2.8455283490e-2,
        3.6585364490e-2,    4.8780485990e-2,    5.6910566990e-2,    6.9105684760e-2,
        7.7235758300e-2,    8.5365831850e-2,    8.9430868630e-2,    9.3495905400e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0406503980e-1,    2.0813007650e-1,
        2.1626015010e-1,    2.2439022360e-1,    2.3252029720e-1,    2.4471540750e-1,
        2.5284549590e-1,    2.6504060630e-1,    2.7723571660e-1,    2.8536579010e-1,
        2.9756090040e-1,    3.0975601080e-1,    3.2195112110e-1,    3.3414623140e-1,
        3.4634134170e-1,    3.5853645210e-1,    3.7073156240e-1,    3.7886163590e-1,
        3.9105674620e-1,    4.0325185660e-1,    4.1544696690e-1,    4.2764207720e-1,
        4.3983718750e-1,    4.5203229780e-1,    4.6422740820e-1,    4.7642251850e-1,
        4.8861762880e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n82_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        1.9782936570e0,     4.0621500020e0,     2.2846693990e0,    -6.7269554140e0,
        -1.7851175310e1,    -1.6274482730e1,     9.3561611180e0,     4.6915527340e1,
        5.7833206180e1,     8.1697235110e0,    -8.4737709050e1,    -1.4084736630e2,
        -7.6152084350e1,     1.0222802730e2,     2.5922463990e2,     2.2239753720e2,
        -4.8553527830e1,    -3.6949105830e2,    -4.4357431030e2,    -1.2503509520e2,
        3.9348999020e2,     6.8104602050e2,     4.2491094970e2,    -2.5155340580e2,
        -8.2414166260e2,    -7.8109436040e2,    -8.1755676270e1,     7.5577368160e2,
        1.0526763920e3,     5.3440301510e2,    -4.2280413820e2,    -1.0838450930e3,
        -9.4145648190e2,    -1.0858459470e2,     7.9078906250e2,     1.1128549800e3,
        6.5868872070e2,    -2.2532189940e2,    -9.3132897950e2,    -1.0079125370e3,
        -4.2736737060e2,
    ];
    let expected_deviations = vec![
        6.8845963480e-1,    6.8845963480e-1,
    ];
    let expected_extremal_frequencies = vec![
        8.1300811140e-4,    7.3170736430e-3,    1.7073171210e-2,    2.7642266820e-2,
        3.9024371650e-2,    4.9593467270e-2,    5.9349555520e-2,    6.9918654860e-2,
        7.8861735760e-2,    8.6991809310e-2,    9.3495868150e-2,    9.8373912270e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0162601770e-1,    2.0650406180e-1,
        2.1382112800e-1,    2.2195120160e-1,    2.3170728980e-1,    2.4227638540e-1,
        2.5284549590e-1,    2.6341459160e-1,    2.7479669450e-1,    2.8617879750e-1,
        2.9756090040e-1,    3.0894300340e-1,    3.2113811370e-1,    3.3252021670e-1,
        3.4471532700e-1,    3.5609743000e-1,    3.6829254030e-1,    3.8048765060e-1,
        3.9186975360e-1,    4.0406486390e-1,    4.1625997420e-1,    4.2764207720e-1,
        4.3983718750e-1,    4.5203229780e-1,    4.6422740820e-1,    4.7560951110e-1,
        4.8780462150e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n82_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        2.0226039890e0,     4.1507129670e0,     2.3266849520e0,    -6.8909912110e0,
        -1.8250465390e1,    -1.6599645610e1,     9.6382274630e0,     4.8003936770e1,
        5.9033210750e1,     8.1142501830e0,    -8.6837654110e1,    -1.4387457280e2,
        -7.7293746950e1,     1.0515885930e2,     2.6504568480e2,     2.2635598750e2,
        -5.1156341550e1,    -3.7837225340e2,    -4.5209628300e2,    -1.2496905520e2,
        4.0425677490e2,     6.9496777340e2,     4.2994635010e2,    -2.6142022710e2,
        -8.4233172610e2,    -7.9276782230e2,    -7.6325561520e1,     7.7484576420e2,
        1.0701071780e3,     5.3596057130e2,    -4.3813592530e2,    -1.1033132320e3,
        -9.4982904050e2,    -1.0076855470e2,     8.0687768550e2,     1.1245501710e3,
        6.5922912600e2,    -2.3328759770e2,    -9.4072070310e2,    -1.0139204100e3,
        -4.2922985840e2,
    ];
    let expected_deviations = vec![
        7.0434761050e-1,    7.0434761050e-1,
    ];
    let expected_extremal_frequencies = vec![
        7.6219509360e-4,    6.8597560750e-3,    1.7530487850e-2,    2.8201209380e-2,
        3.8871932770e-2,    4.9542654310e-2,    5.9451181440e-2,    6.9359712300e-2,
        7.8506045040e-2,    8.6890183390e-2,    9.3749932940e-2,    9.8323099320e-2,
        1.0000000150e-1,    2.0000000300e-1,    2.0152439180e-1,    2.0685975250e-1,
        2.1371950210e-1,    2.2210364040e-1,    2.3201216760e-1,    2.4192069470e-1,
        2.5259143110e-1,    2.6402434710e-1,    2.7469506860e-1,    2.8612798450e-1,
        2.9756090040e-1,    3.0899381640e-1,    3.2118892670e-1,    3.3262184260e-1,
        3.4481695290e-1,    3.5624986890e-1,    3.6844497920e-1,    3.7987789510e-1,
        3.9207300540e-1,    4.0426811580e-1,    4.1570103170e-1,    4.2789614200e-1,
        4.4009125230e-1,    4.5228636260e-1,    4.6371927860e-1,    4.7591438890e-1,
        4.8810949920e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,3,2,2
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n128_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        4.6315330270e-1,     1.2053865430e1,     2.1203479770e1,     1.0022430420e0,
        -8.3167747500e1,    -1.9506823730e2,    -1.4990097050e2,     2.8398812870e2,
        9.9761596680e2,     1.1777539060e3,    -3.0475146480e2,    -3.3989362790e3,
        -5.4018242190e3,    -1.9761840820e3,     8.1122797850e3,     1.7827246090e4,
        1.3538777340e4,    -1.2318437500e4,    -4.5497359380e4,    -5.0120125000e4,
        2.4312578120e3,     9.1624335940e4,     1.3728940620e5,     5.5921906250e4,
        -1.4132965620e5,    -3.0196093750e5,    -2.2414432810e5,     1.4047831250e5,
        5.4800793750e5,     5.8113412500e5,     2.1685750000e4,    -8.1731562500e5,
        -1.1863566250e6,    -5.0629637500e5,     9.5241900000e5,     2.0157540000e6,
        1.4753283750e6,    -6.9124487500e5,    -2.8931340000e6,    -2.9962430000e6,
        -2.7486500000e5,     3.4612567500e6,     4.9304540000e6,     2.1721525000e6,
        -3.2390437500e6,    -6.8660775000e6,    -4.9806690000e6,     1.7808970000e6,
        8.1602370000e6,     8.3139680000e6,     1.1030720000e6,    -8.1186810000e6,
        -1.1429173000e7,    -5.1619420000e6,     6.2692330000e6,     1.3410302000e7,
        9.6639260000e6,    -2.6170120000e6,    -1.3482724000e7,    -1.3562029000e7,
        -2.2509390000e6,     1.1334334000e7,     1.5828818000e7,     7.2914480000e6,
    ];
    let expected_deviations = vec![
        6.8185679610e-2,    6.8185679610e-2,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    7.8125000000e-3,    1.1718750000e-2,    1.9531250000e-2,
        2.3437500000e-2,    3.1250000000e-2,    3.9062500000e-2,    4.2968750000e-2,
        5.0781250000e-2,    5.4687500000e-2,    6.2500000000e-2,    6.6406250000e-2,
        7.4218750000e-2,    7.8125000000e-2,    8.2031250000e-2,    8.5937500000e-2,
        8.9843750000e-2,    9.3750000000e-2,    1.0000000150e-1,    2.0000000300e-1,
        2.0390625300e-1,    2.0781250300e-1,    2.1171875300e-1,    2.1562500300e-1,
        2.1953125300e-1,    2.2343750300e-1,    2.2734375300e-1,    2.3125000300e-1,
        2.3906250300e-1,    2.4296875300e-1,    2.5078123810e-1,    2.5859373810e-1,
        2.6249998810e-1,    2.7031248810e-1,    2.7812498810e-1,    2.8593748810e-1,
        2.8984373810e-1,    2.9765623810e-1,    3.0546873810e-1,    3.1328123810e-1,
        3.2109373810e-1,    3.2890623810e-1,    3.3671873810e-1,    3.4062498810e-1,
        3.4843748810e-1,    3.5624998810e-1,    3.6406248810e-1,    3.7187498810e-1,
        3.7968748810e-1,    3.8749998810e-1,    3.9531248810e-1,    4.0312498810e-1,
        4.1093748810e-1,    4.1874998810e-1,    4.2265623810e-1,    4.3046873810e-1,
        4.3828123810e-1,    4.4609373810e-1,    4.5390623810e-1,    4.6171873810e-1,
        4.6953123810e-1,    4.7734373810e-1,    4.8515623810e-1,    4.9296873810e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,3,2,3
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n128_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.8777949810e0,     1.3369090080e1,     2.0522821430e1,    -6.2504768370e0,
        -9.2218307500e1,    -1.8215560910e2,    -9.7442871090e1,     3.2414706420e2,
        8.9059631350e2,     8.6358593750e2,    -5.4334057620e2,    -2.9453564450e3,
        -3.9799438480e3,    -6.4249853520e2,     7.0180908200e3,     1.2947566410e4,
        7.7878154300e3,    -1.1530922850e4,    -3.2542949220e4,    -3.0868488280e4,
        8.2657265620e3,     6.5069562500e4,     8.5546875000e4,     2.3268421880e4,
        -1.0191095310e5,    -1.8767675000e5,    -1.1867978120e5,     1.1182187500e5,
        3.3905653120e5,     3.2218481250e5,    -3.2400375000e4,    -5.0611803120e5,
        -6.6522456250e5,    -2.2459931250e5,     6.0111400000e5,     1.1317203750e6,
        7.4469006250e5,    -4.8202037500e5,    -1.6232057500e6,    -1.5593505000e6,
        -1.3827500000e4,     1.9470805000e6,     2.5891385000e6,     1.0007085000e6,
        -1.8504665000e6,    -3.6127142500e6,    -2.4608482500e6,     1.1051317500e6,
        4.2935085000e6,     4.1865437500e6,     3.7996350000e5,    -4.2751600000e6,
        -5.7904820000e6,    -2.4677780000e6,     3.3211190000e6,     6.8028620000e6,
        4.7777270000e6,    -1.4413070000e6,    -6.8315970000e6,    -6.7732630000e6,
        -1.0590425000e6,     5.7226400000e6,     7.9319390000e6,     3.6458305000e6,
    ];
    let expected_deviations = vec![
        1.7543707790e-1,    1.7543707790e-1,
    ];
    let expected_extremal_frequencies = vec![
        2.6041667440e-3,    5.2083334890e-3,    1.0416666980e-2,    1.8229166050e-2,
        2.3437498140e-2,    3.1249996270e-2,    3.6458332090e-2,    4.4270835820e-2,
        4.9479171630e-2,    5.7291675360e-2,    6.2500007450e-2,    6.7708335820e-2,
        7.5520828370e-2,    8.0729156730e-2,    8.5937485100e-2,    8.8541649280e-2,
        9.1145813470e-2,    9.3749977650e-2,    9.6354141830e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0260417460e-1,    2.0520834620e-1,    2.1041668950e-1,
        2.1302086110e-1,    2.1822920440e-1,    2.2343754770e-1,    2.2864589100e-1,
        2.3645840590e-1,    2.4166674910e-1,    2.4947926400e-1,    2.5468757750e-1,
        2.6250004770e-1,    2.6770836110e-1,    2.7552083130e-1,    2.8333330150e-1,
        2.9114577170e-1,    2.9635408520e-1,    3.0416655540e-1,    3.1197902560e-1,
        3.1979149580e-1,    3.2760396600e-1,    3.3281227950e-1,    3.4062474970e-1,
        3.4843721990e-1,    3.5624969010e-1,    3.6406216030e-1,    3.7187463050e-1,
        3.7968710060e-1,    3.8489541410e-1,    3.9270788430e-1,    4.0052035450e-1,
        4.0833282470e-1,    4.1614529490e-1,    4.2395776510e-1,    4.3177023530e-1,
        4.3958270550e-1,    4.4739517570e-1,    4.5520764590e-1,    4.6302011610e-1,
        4.7083258630e-1,    4.7604089980e-1,    4.8385337000e-1,    4.9166584010e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,3,2,15
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n128_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        6.5043425560e0,     2.6569074630e1,     3.4865982060e1,    -2.5386058810e1,
        -1.8161166380e2,    -3.0467391970e2,    -8.9887939450e1,     6.4393041990e2,
        1.4397846680e3,     1.1031840820e3,    -1.2942502440e3,    -4.6574921880e3,
        -5.2810419920e3,     4.9045751950e2,     1.1109449220e4,     1.7201378910e4,
        7.2402421880e3,    -1.9342621090e4,    -4.3037390620e4,    -3.3708472660e4,
        2.0300636720e4,     8.6172328120e4,     9.7087812500e4,     8.7129062500e3,
        -1.3784220310e5,    -2.1550851560e5,    -1.0656559380e5,     1.6483706250e5,
        3.9155575000e5,     3.2057028120e5,    -1.0138200000e5,    -5.9020812500e5,
        -6.8351981250e5,    -1.4400062500e5,     7.2062268750e5,     1.1781505000e6,
        6.5805643750e5,    -6.3800850000e5,    -1.7032557500e6,    -1.4708312500e6,
        1.7939175000e5,     2.0645120000e6,     2.5003937500e6,     7.6695675000e5,
        -2.0107893750e6,    -3.5252350000e6,    -2.1803050000e6,     1.3191210000e6,
        4.2138005000e6,     3.8539600000e6,     9.6223500000e4,    -4.2195610000e6,
        -5.4089590000e6,    -2.0984090000e6,     3.3159952500e6,     6.3904575000e6,
        4.3166310000e6,    -1.5163660000e6,    -6.4221995000e6,    -6.2327910000e6,
        -8.8236200000e5,     5.3586835000e6,     7.3450060000e6,     3.3651895000e6,
    ];
    let expected_deviations = vec![
        6.8421661850e-1,    6.8421661850e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.2083336050e-4,    4.6875001860e-3,    1.0937502610e-2,    1.7708336930e-2,
        2.4479163810e-2,    3.1249990690e-2,    3.8020819430e-2,    4.4791646300e-2,
        5.1562473180e-2,    5.7812467220e-2,    6.4062461260e-2,    7.0312455300e-2,
        7.6041616500e-2,    8.1770777700e-2,    8.6458273230e-2,    9.1145768760e-2,
        9.4791598620e-2,    9.7395762800e-2,    9.9479094150e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0104168360e-1,    2.0312504470e-1,    2.0677092670e-1,
        2.1093764900e-1,    2.1614605190e-1,    2.2135445480e-1,    2.2708369790e-1,
        2.3333378140e-1,    2.4010470510e-1,    2.4635478850e-1,    2.5312560800e-1,
        2.6041716340e-1,    2.6718789340e-1,    2.7447944880e-1,    2.8125017880e-1,
        2.8854173420e-1,    2.9583328960e-1,    3.0312484500e-1,    3.1041640040e-1,
        3.1822878120e-1,    3.2552033660e-1,    3.3281189200e-1,    3.4062427280e-1,
        3.4791582820e-1,    3.5520738360e-1,    3.6301976440e-1,    3.7031131980e-1,
        3.7812370060e-1,    3.8541525600e-1,    3.9322763680e-1,    4.0104001760e-1,
        4.0833157300e-1,    4.1614395380e-1,    4.2343550920e-1,    4.3124789000e-1,
        4.3906027080e-1,    4.4635182620e-1,    4.5416420700e-1,    4.6197658780e-1,
        4.6926814320e-1,    4.7708052400e-1,    4.8489290480e-1,    4.9218446020e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,3,2,16
0.000000000000,0.100000001490,0.200000002980,0.500000000000
1.000000000000,0.000000000000
1.000000000000,1.000000000000
*/
#[test]
fn hilbert_lowpass_n128_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.1f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.2f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 1f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        6.7471976280e0,     2.7314832690e1,     3.5928295140e1,    -2.6039844510e1,
        -1.8718376160e2,    -3.1466748050e2,    -9.3799743650e1,     6.6390417480e2,
        1.4881719970e3,     1.1442590330e3,    -1.3323229980e3,    -4.8160878910e3,
        -5.4747509770e3,     4.8691210940e2,     1.1488616210e4,     1.7832777340e4,
        7.5551621090e3,    -1.9989753910e4,    -4.4621593750e4,    -3.5061152340e4,
        2.0901800780e4,     8.9344390620e4,     1.0091352340e5,     9.3754375000e3,
        -1.4287407810e5,    -2.2395376560e5,    -1.1126700000e5,     1.7065292190e5,
        4.0685781250e5,     3.3399700000e5,    -1.0422837500e5,    -6.1316906250e5,
        -7.1169525000e5,    -1.5162631250e5,     7.4831475000e5,     1.2263941250e6,
        6.8724193750e5,    -6.6152325000e5,    -1.7727152500e6,    -1.5339325000e6,
        1.8304350000e5,     2.1482675000e6,     2.6063880000e6,     8.0354350000e5,
        -2.0914316250e6,    -3.6738665000e6,    -2.2766810000e6,     1.3698685000e6,
        4.3909110000e6,     4.0210332500e6,     1.0580450000e5,    -4.3963870000e6,
        -5.6416955000e6,    -2.1930315000e6,     3.4541682500e6,     6.6646090000e6,
        4.5053770000e6,    -1.5780250000e6,    -6.6975675000e6,    -6.5028285000e6,
        -9.2254900000e5,     5.5888660000e6,     7.6622290000e6,     3.5107570000e6,
    ];
    let expected_deviations = vec![
        7.0052212480e-1,    7.0052212480e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.8828125000e-4,    4.3945312500e-3,    1.0742187500e-2,    1.7578125000e-2,
        2.4414062500e-2,    3.1250000000e-2,    3.8085937500e-2,    4.4921875000e-2,
        5.1269531250e-2,    5.8105468750e-2,    6.4453125000e-2,    7.0312500000e-2,
        7.6171875000e-2,    8.1542968750e-2,    8.6425781250e-2,    9.1308593750e-2,
        9.4726562500e-2,    9.7656250000e-2,    9.9121093750e-2,    1.0000000150e-1,
        2.0000000300e-1,    2.0097656550e-1,    2.0341797170e-1,    2.0683594050e-1,
        2.1123047170e-1,    2.1611328420e-1,    2.2148437800e-1,    2.2734375300e-1,
        2.3369140920e-1,    2.4003906550e-1,    2.4638672170e-1,    2.5322264430e-1,
        2.6005858180e-1,    2.6738280060e-1,    2.7421873810e-1,    2.8154295680e-1,
        2.8886717560e-1,    2.9619139430e-1,    3.0351561310e-1,    3.1083983180e-1,
        3.1816405060e-1,    3.2548826930e-1,    3.3281248810e-1,    3.4062498810e-1,
        3.4794920680e-1,    3.5527342560e-1,    3.6308592560e-1,    3.7041014430e-1,
        3.7822264430e-1,    3.8554686310e-1,    3.9335936310e-1,    4.0068358180e-1,
        4.0849608180e-1,    4.1582030060e-1,    4.2363280060e-1,    4.3144530060e-1,
        4.3876951930e-1,    4.4658201930e-1,    4.5439451930e-1,    4.6171873810e-1,
        4.6953123810e-1,    4.7685545680e-1,    4.8466795680e-1,    4.9248045680e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n3_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        1.9019217790e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        3.3333334330e-1,    6.6666668650e-1,    3.3333334330e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n3_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        1.9019217790e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        3.3333334330e-1,    6.6666668650e-1,    3.3333334330e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n3_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        1.9019217790e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        3.3333334330e-1,    6.6666668650e-1,    3.3333334330e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n3_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        1.9019217790e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        3.3333334330e-1,    6.6666668650e-1,    3.3333334330e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n3_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        1.9019217790e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        3.3333334330e-1,    6.6666668650e-1,    3.3333334330e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
3,3,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n3_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(3, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        1.9019217790e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        3.3333334330e-1,    6.6666668650e-1,    3.3333334330e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n4_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        5.1640961320e-2,    2.1031600240e-1,
    ];
    let expected_deviations = vec![
        3.1735008960e-1,    6.3470017910e-1,    3.1735008960e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n4_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        5.9558212760e-2,    2.0175990460e-1,
    ];
    let expected_deviations = vec![
        3.2446596030e-1,    6.4893192050e-1,    3.2446596030e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    3.7999999520e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n4_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        5.9558212760e-2,    2.0175990460e-1,
    ];
    let expected_deviations = vec![
        3.2446596030e-1,    6.4893192050e-1,    3.2446596030e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    3.7999999520e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n4_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        5.9558212760e-2,    2.0175990460e-1,
    ];
    let expected_deviations = vec![
        3.2446596030e-1,    6.4893192050e-1,    3.2446596030e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    3.7999999520e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n4_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        5.9558212760e-2,    2.0175990460e-1,
    ];
    let expected_deviations = vec![
        3.2446596030e-1,    6.4893192050e-1,    3.2446596030e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    3.7999999520e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
4,3,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n4_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(4, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        5.9558212760e-2,    2.0175990460e-1,
    ];
    let expected_deviations = vec![
        3.2446596030e-1,    6.4893192050e-1,    3.2446596030e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.7000000180e-1,    3.3000001310e-1,    3.7999999520e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n11_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        8.5080064830e-2,    1.2181410190e-1,   -2.0892080660e-1,   -9.7210831940e-2,
        2.7695780990e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.2071456020e-1,    2.4142912030e-1,    1.2071456020e-1,
    ];
    let expected_extremal_frequencies = vec![
        5.0000000750e-2,    2.1999999880e-1,    2.7000001070e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.3000000720e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n11_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        7.9987883570e-2,    1.0919684920e-1,   -2.1337063610e-1,   -7.2943598030e-2,
        3.0976003410e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.1947672810e-1,    2.3895345630e-1,    1.1947672810e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.3333335070e-2,    1.3333334030e-1,    1.7000000180e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.4666665790e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n11_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        7.6160520320e-2,    1.0844547300e-1,   -2.0751068000e-1,   -7.4056133630e-2,
        3.0133080480e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.2387523800e-1,    2.4775047600e-1,    1.2387523800e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.6666670590e-2,    1.2666668000e-1,    2.7333328130e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.3999993800e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n11_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        7.6207578180e-2,    1.0848072170e-1,   -2.0745301250e-1,   -7.4271380900e-2,
        3.0105191470e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.2394772470e-1,    2.4789544940e-1,    1.2394772470e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.3750002980e-2,    1.2500001490e-1,    2.6999995110e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.3624994160e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n11_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        7.6126784090e-2,    1.0835477710e-1,   -2.0744411650e-1,   -7.4067316950e-2,
        3.0126941200e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.2397852540e-1,    2.4795705080e-1,    1.2397852540e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.4444445520e-2,    1.2777778510e-1,    2.7000012990e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.3833348160e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
11,3,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n11_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(11, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        7.6132349670e-2,    1.0840066520e-1,   -2.0738163590e-1,   -7.4088901280e-2,
        3.0119258170e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.2399457400e-1,    2.4798914790e-1,    1.2399457400e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.5945938680e-2,    1.2702707950e-1,    2.7135136720e-1,    3.3000001310e-1,
        3.7999999520e-1,    4.3675696850e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n12_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        1.0148525470e-3,    1.3810154800e-1,   -4.2693674560e-2,   -2.3439928890e-1,
        1.5488831700e-1,    2.3446805780e-1,
    ];
    let expected_deviations = vec![
        1.0771293940e-1,    2.1542587880e-1,    1.0771293940e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.2500000000e-1,    1.7000000180e-1,    2.1999999880e-1,    2.6166665550e-1,
        3.3000001310e-1,    3.7999999520e-1,    4.2166665200e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n12_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        -1.7492720870e-5,    1.2838177380e-1,   -4.4858485460e-2,   -2.3999787870e-1,
        1.5238027270e-1,    2.2613224390e-1,
    ];
    let expected_deviations = vec![
        1.1173610390e-1,    2.2347220780e-1,    1.1173610390e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.1111111190e-1,    1.7000000180e-1,    2.1999999880e-1,    2.7555555110e-1,
        3.3000001310e-1,    3.7999999520e-1,    4.3555557730e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n12_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        6.5283727600e-4,    1.3068677480e-1,   -4.2927540840e-2,   -2.3614099620e-1,
        1.5337604280e-1,    2.2826468940e-1,
    ];
    let expected_deviations = vec![
        1.1210712790e-1,    2.2421425580e-1,    1.1210712790e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.2222222240e-1,    1.7000000180e-1,    2.1999999880e-1,    2.7555561070e-1,
        3.3000001310e-1,    3.7999999520e-1,    4.2444455620e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n12_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        6.5795029510e-4,    1.3095131520e-1,   -4.2923688890e-2,   -2.3597502710e-1,
        1.5337859090e-1,    2.2853282090e-1,
    ];
    let expected_deviations = vec![
        1.1200077090e-1,    2.2400154170e-1,    1.1200077090e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.2500001490e-1,    1.7000000180e-1,    2.1999999880e-1,    2.7729168530e-1,
        3.3000001310e-1,    3.7999999520e-1,    4.2687508460e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n12_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        6.4578512680e-4,    1.3075524570e-1,   -4.2935281990e-2,   -2.3601388930e-1,
        1.5332953630e-1,    2.2837530080e-1,
    ];
    let expected_deviations = vec![
        1.1212012920e-1,    2.2424025830e-1,    1.1212012920e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.2499994780e-1,    1.7000000180e-1,    2.1999999880e-1,    2.7555552120e-1,
        3.3000001310e-1,    3.7999999520e-1,    4.2629611490e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
12,3,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n12_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(12, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        6.6309620160e-4,    1.3078762590e-1,   -4.2894743380e-2,   -2.3600286250e-1,
        1.5337738390e-1,    2.2837837040e-1,
    ];
    let expected_deviations = vec![
        1.1210513110e-1,    2.2421026230e-1,    1.1210513110e-1,
    ];
    let expected_extremal_frequencies = vec![
        1.2387382240e-1,    1.7000000180e-1,    2.1999999880e-1,    2.7405402060e-1,
        3.3000001310e-1,    3.7999999520e-1,    4.2729726430e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n16_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        -2.7958168650e-3,    8.1366524100e-3,    1.2464383620e-3,    1.4397120480e-1,
        -4.0892757480e-2,   -2.2451120620e-1,    1.5345579390e-1,    2.3679095510e-1,
    ];
    let expected_deviations = vec![
        1.0674790290e-1,    2.1349580590e-1,    1.0674790290e-1,
    ];
    let expected_extremal_frequencies = vec![
        3.1250000000e-2,    1.2500000000e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.8249999880e-1,    3.3000001310e-1,    3.7999999520e-1,    4.1124999520e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n16_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        -1.5485560290e-3,    5.9543275270e-3,    7.9025141890e-4,    1.4213883880e-1,
        -4.3051585560e-2,   -2.2511547800e-1,    1.5072363620e-1,    2.3806822300e-1,
    ];
    let expected_deviations = vec![
        1.0826438670e-1,    2.1652877330e-1,    1.0826438670e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.1666667910e-2,    1.2500000000e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.8249999880e-1,    3.3000001310e-1,    3.7999999520e-1,    4.2166668180e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n16_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        -1.7096669180e-3,    6.8232025950e-3,    8.4979133680e-4,    1.4120654760e-1,
        -4.2882762850e-2,   -2.2419254480e-1,    1.5008953210e-1,    2.3713529110e-1,
    ];
    let expected_deviations = vec![
        1.0925100740e-1,    2.1850201490e-1,    1.0925100740e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.5833334330e-2,    1.2916670740e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.7416661380e-1,    3.3000001310e-1,    3.7999999520e-1,    4.2166662220e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n16_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        -1.7274932470e-3,    6.8330420180e-3,    8.4711704400e-4,    1.4121329780e-1,
        -4.2855732140e-2,   -2.2416357700e-1,    1.5011608600e-1,    2.3712465170e-1,
    ];
    let expected_deviations = vec![
        1.0925470290e-1,    2.1850940590e-1,    1.0925470290e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.2968750000e-2,    1.2890625000e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.7468749880e-1,    3.3000001310e-1,    3.7999999520e-1,    4.2296874520e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n16_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        -1.7458240040e-3,    6.7760339010e-3,    8.1680994480e-4,    1.4120420810e-1,
        -4.2852304880e-2,   -2.2414086760e-1,    1.5013492110e-1,    2.3714976010e-1,
    ];
    let expected_deviations = vec![
        1.0927081850e-1,    2.1854163710e-1,    1.0927081850e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.5138902960e-2,    1.2673614920e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.7555534240e-1,    3.3000001310e-1,    3.7999999520e-1,    4.2166650300e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
16,3,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n16_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(16, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        -1.7433619360e-3,    6.7734932530e-3,    8.1907538700e-4,    1.4121051130e-1,
        -4.2851217090e-2,   -2.2414314750e-1,    1.5013803540e-1,    2.3715499040e-1,
    ];
    let expected_deviations = vec![
        1.0926641520e-1,    2.1853283050e-1,    1.0926641520e-1,
    ];
    let expected_extremal_frequencies = vec![
        4.3918918820e-2,    1.2668915090e-1,    1.7000000180e-1,    2.1999999880e-1,
        2.7405425910e-1,    3.3000001310e-1,    3.7999999520e-1,    4.2222988610e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n47_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        -4.2875455690e-3,   -1.3024343930e-4,    7.4382349850e-3,   -1.8478679700e-3,
        -1.6678662510e-3,   -2.3888400760e-3,   -9.3934237960e-3,    1.3832369820e-2,
        1.2827303260e-2,   -1.6799060630e-2,   -2.9182110450e-3,   -8.4048863500e-3,
        -1.1061532420e-3,    4.6484947200e-2,   -1.3558921400e-2,   -4.8747919500e-2,
        1.5251904730e-2,   -9.1761853550e-3,    5.1756050440e-2,    7.7620908620e-2,
        -1.8698312340e-1,   -7.8235968950e-2,    3.0279129740e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.3681572640e-3,    4.7363145280e-3,    2.3681572640e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0869565420e-2,    3.2608695330e-2,    5.4347828030e-2,    7.6086953280e-2,
        8.6956515910e-2,    1.0869564120e-1,    1.3043476640e-1,    1.4130432900e-1,
        1.5217389170e-1,    1.7000000180e-1,    2.1999999880e-1,    2.5260868670e-1,
        2.7434781190e-1,    2.9608693720e-1,    3.0695649980e-1,    3.1782606240e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.9086955790e-1,    4.0173912050e-1,
        4.1260868310e-1,    4.3434780840e-1,    4.5608693360e-1,    4.7782605890e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n47_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        -4.1174893270e-3,   -1.1414060140e-3,    6.4768819140e-3,   -1.3407034570e-3,
        3.6492268550e-4,   -2.0954690410e-3,   -1.1764832770e-2,    1.1793521230e-2,
        1.4253500850e-2,   -1.3911210000e-2,   -2.6257422290e-3,   -1.0691380130e-2,
        -3.1902745830e-3,    4.7036401930e-2,   -1.0639603250e-2,   -4.7287732360e-2,
        1.3072527950e-2,   -1.1812273410e-2,    5.1716580990e-2,    8.0446481700e-2,
        -1.8481002750e-1,   -8.0109529200e-2,    2.9953926800e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        3.6452559290e-3,    7.2905118580e-3,    3.6452559290e-3,
    ];
    let expected_extremal_frequencies = vec![
        7.2463769470e-3,    2.8985507790e-2,    5.0724633040e-2,    7.2463758290e-2,
        9.4202883540e-2,    1.1594200880e-1,    1.3043476640e-1,    1.5217389170e-1,
        1.5942026670e-1,    1.7000000180e-1,    2.1999999880e-1,    2.2724637390e-1,
        2.6347824930e-1,    2.8521737460e-1,    3.0695649980e-1,    3.2144925000e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8724637030e-1,    4.0173912050e-1,
        4.1623187070e-1,    4.3797099590e-1,    4.5971012120e-1,    4.8869562150e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n47_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        -3.7357213440e-3,   -1.8072654490e-3,    5.7804938410e-3,   -7.1770721120e-4,
        1.2595418370e-3,   -2.2839258890e-3,   -1.2491545640e-2,    1.1418004520e-2,
        1.4683376070e-2,   -1.3178605590e-2,   -2.5276360100e-3,   -1.1010547170e-2,
        -4.1078748180e-3,    4.6806558970e-2,   -9.3196742240e-3,   -4.6565577390e-2,
        1.2037731710e-2,   -1.2866223230e-2,    5.1912125200e-2,    8.1709966060e-2,
        -1.8424670400e-1,   -8.1032238900e-2,    2.9859384890e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        4.2759054340e-3,    8.5518108680e-3,    4.2759054340e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0144928470e-2,    3.1884063040e-2,    5.2173923700e-2,    7.2463758290e-2,
        9.4202838840e-2,    1.1304337530e-1,    1.3333319130e-1,    1.5072445570e-1,
        1.6376790400e-1,    1.7000000180e-1,    2.1999999880e-1,    2.2724635900e-1,
        2.4608689550e-1,    2.8376832600e-1,    3.0695691700e-1,    3.2289907340e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8579714300e-1,    4.0029001240e-1,
        4.1913074250e-1,    4.3942075970e-1,    4.6260935070e-1,    4.8724722860e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n47_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        -3.7216637280e-3,   -1.8269321880e-3,    5.7640424930e-3,   -7.0123188200e-4,
        1.2708089780e-3,   -2.2848346270e-3,   -1.2507061470e-2,    1.1401759460e-2,
        1.4699134980e-2,   -1.3159774240e-2,   -2.5270609190e-3,   -1.1011540890e-2,
        -4.1379001920e-3,    4.6801321210e-2,   -9.2762485150e-3,   -4.6558894220e-2,
        1.2011158280e-2,   -1.2884045020e-2,    5.1906019450e-2,    8.1743150950e-2,
        -1.8423134090e-1,   -8.1064775590e-2,    2.9858025910e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        4.2865448630e-3,    8.5730897260e-3,    4.2865448630e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0869564490e-2,    3.1249994410e-2,    5.2989121530e-2,    7.3369555180e-2,
        9.3749985100e-2,    1.1413041500e-1,    1.3315218690e-1,    1.5081532300e-1,
        1.6440235080e-1,    1.7000000180e-1,    2.1999999880e-1,    2.2679351270e-1,
        2.4581535160e-1,    2.8385865690e-1,    3.0695635080e-1,    3.2326060530e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8543474670e-1,    4.0038031340e-1,
        4.1804325580e-1,    4.3978226180e-1,    4.6287995580e-1,    4.8733633760e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n47_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        -3.7242469840e-3,   -1.8222854710e-3,    5.7667233050e-3,   -7.0212257560e-4,
        1.2687330600e-3,   -2.2894367580e-3,   -1.2502714060e-2,    1.1407723650e-2,
        1.4694452290e-2,   -1.3166314920e-2,   -2.5263442660e-3,   -1.1009702460e-2,
        -4.1322419420e-3,    4.6799793840e-2,   -9.2853968960e-3,   -4.6556893740e-2,
        1.2016091500e-2,   -1.2881213800e-2,    5.1910713320e-2,    8.1734210250e-2,
        -1.8423540890e-1,   -8.1056833270e-2,    2.9858410360e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        4.2884098370e-3,    8.5768196730e-3,    4.2884098370e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0265701450e-2,    3.1400948760e-2,    5.2536252890e-2,    7.3067627850e-2,
        9.3598939480e-2,    1.1352638900e-1,    1.3285006580e-1,    1.5036228300e-1,
        1.6425128280e-1,    1.7000000180e-1,    2.1999999880e-1,    2.2724643350e-1,
        2.4657025930e-1,    2.8400933740e-1,    3.0695581440e-1,    3.2325989010e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8603854180e-1,    3.9992719890e-1,
        4.1864669320e-1,    4.3978160620e-1,    4.6212422850e-1,    4.8688226940e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
47,3,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n47_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(47, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        -3.7257049700e-3,   -1.8205387750e-3,    5.7660043240e-3,   -7.0362456610e-4,
        1.2718590440e-3,   -2.2907417730e-3,   -1.2503085660e-2,    1.1410205620e-2,
        1.4691345390e-2,   -1.3164999890e-2,   -2.5235526260e-3,   -1.1014313440e-2,
        -4.1300929150e-3,    4.6801354740e-2,   -9.2889750380e-3,   -4.6553470190e-2,
        1.2016001160e-2,   -1.2883996590e-2,    5.1915347580e-2,    8.1732481720e-2,
        -1.8423691390e-1,   -8.1052571540e-2,    2.9858005050e-1,     0.0000000000e0,
    ];
    let expected_deviations = vec![
        4.2886626910e-3,    8.5773253810e-3,    4.2886626910e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0575793680e-2,    3.1727366150e-2,    5.2291452880e-2,    7.2855539620e-2,
        9.3419626360e-2,    1.1398371310e-1,    1.3278506700e-1,    1.5041120350e-1,
        1.6451211270e-1,    1.7000000180e-1,    2.1999999880e-1,    2.2705045340e-1,
        2.4643920360e-1,    2.8404247760e-1,    3.0695703630e-1,    3.2340851430e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8587552310e-1,    3.9997679000e-1,
        4.1877847910e-1,    4.3934282660e-1,    4.6225738530e-1,    4.8693460230e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n48_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        -1.1134617960e-3,   -4.0120380000e-3,    3.6471756180e-3,    3.0782162210e-3,
        -1.3674228680e-3,    1.3736193070e-3,   -1.0015378710e-2,   -2.1229432900e-3,
        2.0600700750e-2,   -3.1476672740e-3,   -1.1348270810e-2,   -9.1012008490e-4,
        -1.6957942400e-2,    2.7134532110e-2,    3.2301738860e-2,   -4.7877665610e-2,
        -1.4143575910e-2,    7.7988747510e-3,   -2.5266148150e-3,    1.0527963940e-1,
        -4.2931605130e-2,   -2.1863876280e-1,    1.5120080110e-1,    2.3798161750e-1,
    ];
    let expected_deviations = vec![
        2.8177842030e-3,    5.6355684060e-3,    2.8177842030e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0416666980e-2,    3.1250000000e-2,    5.2083335820e-2,    7.2916664180e-2,
        9.3749992550e-2,    1.1458332090e-1,    1.3541665670e-1,    1.4583332840e-1,
        1.5625000000e-1,    1.7000000180e-1,    2.1999999880e-1,    2.3041667040e-1,
        2.4083334210e-1,    2.7208331230e-1,    3.0333328250e-1,    3.1374993920e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.9041665200e-1,    4.0083330870e-1,
        4.2166662220e-1,    4.3208327890e-1,    4.6333324910e-1,    4.8416656260e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n48_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        -1.4412783780e-3,   -4.4335536660e-3,    3.7075176370e-3,    2.9828562400e-3,
        -1.3791574170e-3,    1.7955743240e-3,   -1.0169630870e-2,   -2.6943124830e-3,
        2.0805662500e-2,   -2.6936032810e-3,   -1.1185985990e-2,   -1.2364289720e-3,
        -1.7199371010e-2,    2.6953067630e-2,    3.2448355110e-2,   -4.7308605160e-2,
        -1.4101311560e-2,    7.2974208740e-3,   -2.7308538560e-3,    1.0530731830e-1,
        -4.2518485340e-2,   -2.1832875910e-1,    1.5072512630e-1,    2.3750850560e-1,
    ];
    let expected_deviations = vec![
        3.6222611090e-3,    7.2445222180e-3,    3.6222611090e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.3888888990e-2,    3.4722223880e-2,    5.5555555970e-2,    6.9444447760e-2,
        9.0277791020e-2,    1.1111113430e-1,    1.3194447760e-1,    1.4583337310e-1,
        1.5972226860e-1,    1.7000000180e-1,    2.1999999880e-1,    2.2694444660e-1,
        2.4777778980e-1,    2.8250002860e-1,    3.0333337190e-1,    3.1722226740e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8694444300e-1,    4.0083333850e-1,
        4.1472223400e-1,    4.3555557730e-1,    4.5638892050e-1,    4.7722226380e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n48_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        -1.0702564610e-3,   -4.8450222240e-3,    3.2488510480e-3,    3.3094743270e-3,
        -1.2476942500e-3,    2.1848124450e-3,   -9.9444380030e-3,   -3.9875479420e-3,
        2.0748328420e-2,   -1.5238290650e-3,   -1.1379895730e-2,   -1.1249487290e-3,
        -1.7488390210e-2,    2.5747869160e-2,    3.3539336170e-2,   -4.6264797450e-2,
        -1.4966279270e-2,    7.0322025570e-3,   -3.3020209520e-3,    1.0544602570e-1,
        -4.0855608880e-2,   -2.1888706090e-1,    1.4961552620e-1,    2.3768176140e-1,
    ];
    let expected_deviations = vec![
        4.2570433580e-3,    8.5140867160e-3,    4.2570433580e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.1111111380e-2,    3.1944446270e-2,    5.2777778360e-2,    7.3611140250e-2,
        9.3055635690e-2,    1.1388902370e-1,    1.3333351910e-1,    1.5000022950e-1,
        1.6388915480e-1,    1.7000000180e-1,    2.1999999880e-1,    2.2694446150e-1,
        2.4500006440e-1,    2.8249982000e-1,    3.0611073970e-1,    3.2277727130e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8555550580e-1,    3.9944428210e-1,
        4.1749969120e-1,    4.3694397810e-1,    4.5777714250e-1,    4.7861030700e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n48_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        -1.0563245740e-3,   -4.8448089510e-3,    3.2266259660e-3,    3.3155628480e-3,
        -1.2338503730e-3,    2.1903009620e-3,   -9.9436696620e-3,   -4.0112910790e-3,
        2.0745508370e-2,   -1.5008868650e-3,   -1.1378241700e-2,   -1.1198758150e-3,
        -1.7501339320e-2,    2.5720387700e-2,    3.3570837230e-2,   -4.6238303180e-2,
        -1.4989929270e-2,    7.0195980370e-3,   -3.3158697190e-3,    1.0545662050e-1,
        -4.0816571560e-2,   -2.1890275180e-1,    1.4959029850e-1,    2.3768481610e-1,
    ];
    let expected_deviations = vec![
        4.2562577870e-3,    8.5125155750e-3,    4.2562577870e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0416666050e-2,    3.1250003730e-2,    5.2083317190e-2,    7.2916656730e-2,
        9.3750029800e-2,    1.1328131710e-1,    1.3281255960e-1,    1.5104165670e-1,
        1.6406244040e-1,    1.7000000180e-1,    2.1999999880e-1,    2.2651039060e-1,
        2.4604156610e-1,    2.8250011800e-1,    3.0593779680e-1,    3.2286500930e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8520836830e-1,    3.9953139420e-1,
        4.1776070000e-1,    4.3729209900e-1,    4.5812559130e-1,    4.7895908360e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,3,36
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n48_grid36() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 36);

    let expected_impulse_response = vec![
        -1.0590815220e-3,   -4.8496099190e-3,    3.2301109750e-3,    3.3144031190e-3,
        -1.2314276540e-3,    2.1871111820e-3,   -9.9501637740e-3,   -3.9997519930e-3,
        2.0747931670e-2,   -1.5119900930e-3,   -1.1375842620e-2,   -1.1190585790e-3,
        -1.7498223110e-2,    2.5725632910e-2,    3.3558327700e-2,   -4.6237919480e-2,
        -1.4980333860e-2,    7.0123169570e-3,   -3.3109821380e-3,    1.0546087470e-1,
        -4.0829773990e-2,   -2.1889902650e-1,    1.4959698920e-1,    2.3768131430e-1,
    ];
    let expected_deviations = vec![
        4.2665554210e-3,    8.5331108420e-3,    4.2665554210e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0416666050e-2,    3.1249986960e-2,    5.2083373070e-2,    7.2916693990e-2,
        9.3749947850e-2,    1.1400450020e-1,    1.3310165700e-1,    1.5046270190e-1,
        1.6435153780e-1,    1.7000000180e-1,    2.1999999880e-1,    2.2694441680e-1,
        2.4546286460e-1,    2.8192105890e-1,    3.0622652170e-1,    3.2300886510e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8578701020e-1,    3.9967584610e-1,
        4.1761559250e-1,    4.3729144330e-1,    4.5812469720e-1,    4.7895795110e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
48,3,3,37
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n48_grid37() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(48, FilterType::HilbertTransform, &bands, 37);

    let expected_impulse_response = vec![
        -1.0592138860e-3,   -4.8491684720e-3,    3.2311920080e-3,    3.3136599230e-3,
        -1.2324526910e-3,    2.1879011770e-3,   -9.9496925250e-3,   -4.0009263900e-3,
        2.0747274160e-2,   -1.5112934630e-3,   -1.1375172060e-2,   -1.1203954930e-3,
        -1.7498994250e-2,    2.5726241990e-2,    3.3558264370e-2,   -4.6238549050e-2,
        -1.4979826290e-2,    7.0129372180e-3,   -3.3111944790e-3,    1.0546059160e-1,
        -4.0829978880e-2,   -2.1889841560e-1,    1.4959642290e-1,    2.3768083750e-1,
    ];
    let expected_deviations = vec![
        4.2664972130e-3,    8.5329944270e-3,    4.2664972130e-3,
    ];
    let expected_extremal_frequencies = vec![
        1.0698197410e-2,    3.1531520190e-2,    5.2364841100e-2,    7.3198162020e-2,
        9.3468420210e-2,    1.1373867840e-1,    1.3288290800e-1,    1.5090115370e-1,
        1.6441483800e-1,    1.7000000180e-1,    2.1999999880e-1,    2.2731991110e-1,
        2.4533815680e-1,    2.8193688390e-1,    3.0671131610e-1,    3.2303991910e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8563054800e-1,    3.9970692990e-1,
        4.1772469880e-1,    4.3743163350e-1,    4.5770162340e-1,    4.7909772400e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n81_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        1.0597098480e-5,   -8.4706101920e-5,   -9.1991511000e-5,   -2.0441012750e-4,
        4.8682687340e-4,    5.2076048450e-4,   -9.3889085110e-4,   -3.4791810320e-4,
        1.9148201680e-4,    1.3403469350e-4,    2.1879940760e-3,   -1.0575575290e-3,
        -3.7103467620e-3,    1.9609010780e-3,    1.4739816540e-3,    1.1071173940e-3,
        2.3913562760e-3,   -8.2234889270e-3,   -2.3551902270e-3,    1.0584164410e-2,
        -4.3307081800e-4,   -6.7709479480e-5,   -3.5926788110e-3,   -1.5332412910e-2,
        1.5634559090e-2,    1.6880653800e-2,   -1.6305562110e-2,   -2.6130923070e-3,
        -1.2530419980e-2,   -3.6959175490e-3,    5.1180236040e-2,   -1.1327080430e-2,
        -4.9229130150e-2,    1.3020446520e-2,   -1.2617461380e-2,    5.3646609190e-2,
        8.1278257070e-2,   -1.8650171160e-1,   -7.9877577720e-2,    2.9933458570e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        1.4398846540e-4,    2.8797693080e-4,    1.4398846540e-4,
    ];
    let expected_extremal_frequencies = vec![
        6.2500000930e-3,    1.8750000750e-2,    3.1250000000e-2,    4.3750002980e-2,
        5.6250005960e-2,    6.8750008940e-2,    8.1250011920e-2,    9.3750014900e-2,
        1.0625001790e-1,    1.1875002090e-1,    1.3125000890e-1,    1.3750000300e-1,
        1.4999999110e-1,    1.5624998510e-1,    1.6249997910e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2624999280e-1,    2.3249998690e-1,    2.3874998090e-1,
        2.5124996900e-1,    2.6374995710e-1,    2.7624994520e-1,    2.8249993920e-1,
        2.9499992730e-1,    3.0749991540e-1,    3.1374990940e-1,    3.1999990340e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8624998930e-1,    3.9249998330e-1,
        3.9874997740e-1,    4.1124996540e-1,    4.1749995950e-1,    4.2999994750e-1,
        4.4249993560e-1,    4.5499992370e-1,    4.6749991180e-1,    4.7999989990e-1,
        4.9249988790e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n81_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        -6.3372790460e-5,   -1.2299719670e-4,    4.7627880120e-5,   -2.6759126920e-4,
        4.3368773190e-4,    6.2205916040e-4,   -1.1336163150e-3,   -2.9707874640e-4,
        5.1091855860e-4,   -1.3213968490e-4,    2.1324488330e-3,   -9.2659366780e-4,
        -4.0292865600e-3,    2.3666429330e-3,    1.7495182110e-3,    3.9627461230e-4,
        2.5372598320e-3,   -8.0207074060e-3,   -2.6026463600e-3,    1.1284490120e-2,
        -7.3363102270e-4,   -9.9213281650e-4,   -2.8415946290e-3,   -1.5165681020e-2,
        1.5455504880e-2,    1.7452793200e-2,   -1.7337355760e-2,   -2.9539645180e-3,
        -1.1148042980e-2,   -4.0331413040e-3,    5.0950564440e-2,   -1.1142337690e-2,
        -5.0385463980e-2,    1.3842673970e-2,   -1.1441547420e-2,    5.2462209020e-2,
        8.1232465800e-2,   -1.8650011720e-1,   -8.0455049870e-2,    3.0085244770e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.1744934200e-4,    4.3489868400e-4,    2.1744934200e-4,
    ];
    let expected_extremal_frequencies = vec![
        4.1666668840e-3,    1.6666667540e-2,    2.9166666790e-2,    4.1666667910e-2,
        5.4166667160e-2,    6.6666670140e-2,    7.9166680570e-2,    9.1666691010e-2,
        1.0416670140e-1,    1.1666671190e-1,    1.2916670740e-1,    1.4166669550e-1,
        1.5000002090e-1,    1.5833334620e-1,    1.6250000890e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2416666150e-1,    2.2833332420e-1,    2.4083331230e-1,
        2.4916663770e-1,    2.6166662570e-1,    2.7416661380e-1,    2.8666660190e-1,
        2.9916659000e-1,    3.1166657810e-1,    3.1999990340e-1,    3.2416656610e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8416665790e-1,    3.8833332060e-1,
        3.9666664600e-1,    4.0916663410e-1,    4.2166662220e-1,    4.2999994750e-1,
        4.4249993560e-1,    4.5499992370e-1,    4.6749991180e-1,    4.7999989990e-1,
        4.9249988790e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n81_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        -1.7697093430e-5,   -1.4639567230e-4,   -2.0641116860e-5,   -3.3819524110e-4,
        4.5955387760e-4,    8.1565452270e-4,   -1.0916803730e-3,   -5.1133648960e-4,
        3.8278009740e-4,   -1.3047528050e-4,    2.3887776770e-3,   -6.0692301490e-4,
        -4.3557723980e-3,    1.9047727110e-3,    1.8411995840e-3,    7.2045507840e-4,
        3.0201149640e-3,   -8.1022409720e-3,   -3.5466058180e-3,    1.1138504370e-2,
        -1.1686875950e-5,   -4.9909204240e-4,   -2.7146365030e-3,   -1.6112685200e-2,
        1.4580842110e-2,    1.8433509390e-2,   -1.6345595940e-2,   -3.0622158670e-3,
        -1.1855726130e-2,   -5.2424985920e-3,    5.1374007020e-2,   -9.3552861360e-3,
        -5.0367839630e-2,    1.2772220190e-2,   -1.2350745500e-2,    5.2274748680e-2,
        8.2977242770e-2,   -1.8556745350e-1,   -8.1944361330e-2,    2.9980570080e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.9506752620e-4,    5.9013505230e-4,    2.9506752620e-4,
    ];
    let expected_extremal_frequencies = vec![
        6.6666672940e-3,    1.9166665150e-2,    3.1666655090e-2,    4.4166643170e-2,
        5.6666631250e-2,    6.9166623060e-2,    8.1666611140e-2,    9.4166599210e-2,
        1.0666658730e-1,    1.1833324280e-1,    1.2999990580e-1,    1.4166656140e-1,
        1.5166655180e-1,    1.6083320980e-1,    1.6749987010e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2249999640e-1,    2.2999998930e-1,    2.3999997970e-1,
        2.5083330270e-1,    2.6249995830e-1,    2.7499994640e-1,    2.8749993440e-1,
        2.9916659000e-1,    3.0999991300e-1,    3.1999990340e-1,    3.2749989630e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8249999280e-1,    3.8916665320e-1,
        3.9833331110e-1,    4.0833330150e-1,    4.1999995710e-1,    4.3166661260e-1,
        4.4416660070e-1,    4.5583325620e-1,    4.6833324430e-1,    4.8083323240e-1,
        4.9333322050e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
81,3,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n81_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(81, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        -1.9098037230e-5,   -1.4594876850e-4,   -1.8179578550e-5,   -3.3822102710e-4,
        4.5806498380e-4,    8.1436167240e-4,   -1.0920449860e-3,   -5.0771259700e-4,
        3.8495246550e-4,   -1.3394850250e-4,    2.3863378450e-3,   -6.0706166550e-4,
        -4.3528508390e-3,    1.9099437630e-3,    1.8384479920e-3,    7.1384559850e-4,
        3.0200968030e-3,   -8.0988714470e-3,   -3.5412586290e-3,    1.1139722540e-2,
        -2.0422623490e-5,   -5.0293002280e-4,   -2.7092257510e-3,   -1.6107099130e-2,
        1.4583524320e-2,    1.8427006900e-2,   -1.6354335470e-2,   -3.0584123450e-3,
        -1.1847244580e-2,   -5.2390978670e-3,    5.1370091740e-2,   -9.3651255590e-3,
        -5.0369307400e-2,    1.2781919910e-2,   -1.2344235550e-2,    5.2271470430e-2,
        8.2967817780e-2,   -1.8557184930e-1,   -8.1937178970e-2,    2.9981470110e-1,
        0.0000000000e0,
    ];
    let expected_deviations = vec![
        2.9426702530e-4,    5.8853405060e-4,    2.9426702530e-4,
    ];
    let expected_extremal_frequencies = vec![
        6.2500005590e-3,    1.8749998880e-2,    3.1249986960e-2,    4.3749976900e-2,
        5.7031214240e-2,    6.9531232120e-2,    8.1250026820e-2,    9.3750074510e-2,
        1.0625012220e-1,    1.1875016990e-1,    1.3046896460e-1,    1.4140650630e-1,
        1.5234404800e-1,    1.6093783080e-1,    1.6718785460e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2312501070e-1,    2.3015628760e-1,    2.3953132330e-1,
        2.5125008820e-1,    2.6296865940e-1,    2.7468723060e-1,    2.8718703990e-1,
        2.9890561100e-1,    3.1062418220e-1,    3.1999903920e-1,    3.2703018190e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8234370950e-1,    3.8859361410e-1,
        3.9796847110e-1,    4.0890580420e-1,    4.1984313730e-1,    4.3156170850e-1,
        4.4406151770e-1,    4.5656132700e-1,    4.6827989820e-1,    4.8156094550e-1,
        4.9406075480e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n82_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        2.2040254410e-5,   -2.1634553560e-4,    7.0096226410e-5,    7.3784220150e-5,
        -1.3715040400e-4,    6.6967587920e-4,   -3.5723642210e-4,   -1.1488627640e-3,
        8.9074141580e-4,    1.5091251410e-4,    4.2845011920e-4,    1.3156040800e-3,
        -3.5822899080e-3,   -8.0869067460e-4,    4.1644596490e-3,   -6.1562581690e-4,
        1.4637387360e-3,   -1.8569733950e-3,   -8.3599463110e-3,    7.5896349740e-3,
        7.1441344920e-3,   -5.6616538200e-3,    1.0236052330e-3,   -1.0647592130e-2,
        -2.9492061590e-3,    2.6130586860e-2,   -4.4320002200e-3,   -1.6146784650e-2,
        8.4428628910e-4,   -1.6419367860e-2,    2.8461134060e-2,    3.4375656400e-2,
        -5.2604753520e-2,   -1.5594569030e-2,    1.2426983570e-2,   -2.0374413580e-3,
        1.0369757560e-1,   -4.4037953020e-2,   -2.2030322250e-1,    1.5380026400e-1,
        2.4106818440e-1,
    ];
    let expected_deviations = vec![
        1.3074783780e-4,    2.6149567570e-4,    1.3074783780e-4,
    ];
    let expected_extremal_frequencies = vec![
        6.0975607480e-3,    1.8292682250e-2,    3.0487803740e-2,    4.2682923380e-2,
        5.4878048600e-2,    6.7073173820e-2,    7.9268299040e-2,    9.1463424270e-2,
        1.0365854950e-1,    1.1585367470e-1,    1.2804879250e-1,    1.3414634760e-1,
        1.4634145800e-1,    1.5243901310e-1,    1.5853656830e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2609755400e-1,    2.3219510910e-1,    2.3829266430e-1,
        2.5048777460e-1,    2.6268288490e-1,    2.7487799530e-1,    2.8707310560e-1,
        2.9926821590e-1,    3.1146332620e-1,    3.1756088140e-1,    3.2365843650e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8609755040e-1,    3.9219510560e-1,
        3.9829266070e-1,    4.1048777100e-1,    4.1658532620e-1,    4.2878043650e-1,
        4.4097554680e-1,    4.5317065720e-1,    4.6536576750e-1,    4.7756087780e-1,
        4.8975598810e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n82_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        7.2972761700e-5,   -2.2263474240e-4,   -6.2228355090e-7,   -6.2067716500e-5,
        -1.4475670470e-4,    9.2332740310e-4,   -2.3517751830e-4,   -1.3407445510e-3,
        6.1828538310e-4,    1.7593294610e-5,    8.0538936890e-4,    1.8134522250e-3,
        -3.8121026010e-3,   -1.4197950950e-3,    3.8567979350e-3,   -2.1971086970e-4,
        2.3919080850e-3,   -1.7954331120e-3,   -9.3971509490e-3,    6.9092335180e-3,
        7.5349006800e-3,   -4.4327978980e-3,    1.6207576260e-3,   -1.1852283960e-2,
        -4.2660832410e-3,    2.6413466780e-2,   -2.9723364860e-3,   -1.5049616810e-2,
        -1.2894254180e-4,   -1.8370807170e-2,    2.8343493120e-2,    3.6053914580e-2,
        -5.1135689020e-2,   -1.6136424620e-2,    1.0174572470e-2,   -2.8235521170e-3,
        1.0543873160e-1,   -4.2237747460e-2,   -2.2045812010e-1,    1.5166974070e-1,
        2.3962928350e-1,
    ];
    let expected_deviations = vec![
        1.9796773270e-4,    3.9593546530e-4,    1.9796773270e-4,
    ];
    let expected_extremal_frequencies = vec![
        8.1300809980e-3,    2.0325202490e-2,    3.2520323990e-2,    4.4715445490e-2,
        5.6910566990e-2,    6.9105684760e-2,    8.1300795080e-2,    9.3495905400e-2,
        1.0569101570e-1,    1.1788612600e-1,    1.3008123640e-1,    1.3821130990e-1,
        1.5040642020e-1,    1.5853649380e-1,    1.6260153060e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2406503560e-1,    2.3219510910e-1,    2.4032518270e-1,
        2.5252029300e-1,    2.6471540330e-1,    2.7691051360e-1,    2.8910562400e-1,
        3.0130073430e-1,    3.0943080780e-1,    3.2162591810e-1,    3.2569095490e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8406503200e-1,    3.8813006880e-1,
        3.9626014230e-1,    4.0845525260e-1,    4.1658532620e-1,    4.2878043650e-1,
        4.4097554680e-1,    4.5317065720e-1,    4.6536576750e-1,    4.7756087780e-1,
        4.8975598810e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n82_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        1.4027992440e-4,   -1.9637237710e-4,   -8.8150161900e-5,   -1.6862862690e-4,
        -1.6316559050e-4,    1.0723625310e-3,   -8.2066515460e-5,   -1.4625662010e-3,
        3.6963337330e-4,   -4.8812304160e-5,    1.0161768880e-3,    2.1640670020e-3,
        -3.8452679290e-3,   -1.9195094940e-3,    3.6098586860e-3,    9.1592577520e-5,
        2.9335105790e-3,   -1.6318258130e-3,   -1.0041443630e-2,    6.2713841910e-3,
        7.8761624170e-3,   -3.6094104870e-3,    1.9706496970e-3,   -1.2473242360e-2,
        -5.2938070150e-3,    2.6504917070e-2,   -1.8001887950e-3,   -1.4434436340e-2,
        -7.3503097520e-4,   -1.9582549110e-2,    2.8019536290e-2,    3.7330526860e-2,
        -5.0052482630e-2,   -1.6707971690e-2,    8.8348407300e-3,   -3.4234970810e-3,
        1.0647773000e-1,   -4.0747366850e-2,   -2.2071623800e-1,    1.5014578400e-1,
        2.3886702950e-1,
    ];
    let expected_deviations = vec![
        2.4619797480e-4,    4.9239594950e-4,    2.4619797480e-4,
    ];
    let expected_extremal_frequencies = vec![
        6.5040653570e-3,    1.9512193280e-2,    3.2520312820e-2,    4.5528430490e-2,
        5.8536548170e-2,    7.0731662210e-2,    8.3739779890e-2,    9.5934890210e-2,
        1.0813000050e-1,    1.2032511090e-1,    1.3170722130e-1,    1.4227631690e-1,
        1.5284541250e-1,    1.6178849340e-1,    1.6747954490e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2243902090e-1,    2.2975608710e-1,    2.3951217530e-1,
        2.5008127090e-1,    2.6227638130e-1,    2.7447149160e-1,    2.8747960930e-1,
        2.9886171220e-1,    3.1024381520e-1,    3.1999990340e-1,    3.2731696960e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8243901730e-1,    3.8813006880e-1,
        3.9707314970e-1,    4.0682923790e-1,    4.1821134090e-1,    4.2878043650e-1,
        4.4097554680e-1,    4.5235764980e-1,    4.6455276010e-1,    4.7593486310e-1,
        4.8812997340e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
82,3,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n82_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(82, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        1.3988604770e-4,   -1.9609529410e-4,   -8.6859014120e-5,   -1.6964881800e-4,
        -1.6500137280e-4,    1.0727352930e-3,   -8.0467376390e-5,   -1.4601034350e-3,
        3.6800454840e-4,   -5.4194882980e-5,    1.0174575730e-3,    2.1687073170e-3,
        -3.8429731500e-3,   -1.9203351110e-3,    3.6009519830e-3,    8.9380802820e-5,
        2.9451155570e-3,   -1.6277701360e-3,   -1.0046221320e-2,    6.2629068270e-3,
        7.8693469990e-3,   -3.5956287760e-3,    1.9836276770e-3,   -1.2484167700e-2,
        -5.3054429590e-3,    2.6500236240e-2,   -1.7910376190e-3,   -1.4412838970e-2,
        -7.4266968290e-4,   -1.9606253130e-2,    2.8018871320e-2,    3.7341032180e-2,
        -5.0034116950e-2,   -1.6704734410e-2,    8.8042765860e-3,   -3.4324713050e-3,
        1.0649868850e-1,   -4.0733639150e-2,   -2.2071206570e-1,    1.5012386440e-1,
        2.3884476720e-1,
    ];
    let expected_deviations = vec![
        2.4633583960e-4,    4.9267167920e-4,    2.4633583960e-4,
    ];
    let expected_extremal_frequencies = vec![
        6.8597560750e-3,    1.9817071040e-2,    3.2774377610e-2,    4.5731682330e-2,
        5.8688987050e-2,    7.0884101090e-2,    8.3841405810e-2,    9.6036516130e-2,
        1.0823162650e-1,    1.1966454240e-1,    1.3185966010e-1,    1.4253038170e-1,
        1.5243890880e-1,    1.6158524160e-1,    1.6768279670e-1,    1.7000000180e-1,
        2.1999999880e-1,    2.2228658200e-1,    2.2990852590e-1,    2.3905485870e-1,
        2.5048777460e-1,    2.6192069050e-1,    2.7487799530e-1,    2.8707310560e-1,
        2.9926821590e-1,    3.1070113180e-1,    3.2060965900e-1,    3.2746940850e-1,
        3.3000001310e-1,    3.7999999520e-1,    3.8228657840e-1,    3.8838413360e-1,
        3.9676827190e-1,    4.0667679910e-1,    4.1810971500e-1,    4.2878043650e-1,
        4.4097554680e-1,    4.5240846280e-1,    4.6384137870e-1,    4.7603648900e-1,
        4.8823159930e-1,    5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,3,3,2
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n128_grid2() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 2);

    let expected_impulse_response = vec![
        1.0621565710e-6,   -4.0869426810e-6,    4.6370414570e-6,    1.2735723430e-5,
        -2.1584308340e-5,   -7.5195962380e-6,    1.8297832870e-5,   -9.3133803600e-6,
        5.0866758100e-5,   -1.5953080950e-5,   -1.3665339790e-4,    8.9188906710e-5,
        8.7710177470e-5,   -3.4011103710e-5,    9.7571391960e-5,   -3.0231295390e-4,
        -1.4194900000e-4,    5.9676612730e-4,   -5.4191608800e-5,   -2.1451573410e-4,
        -3.7852587410e-5,   -7.5257010760e-4,    8.2779640800e-4,    1.1382985390e-3,
        -1.3078140330e-3,   -3.2493023900e-4,   -1.8876182730e-4,   -3.0294054890e-4,
        3.0140744060e-3,   -6.6331378180e-4,   -3.6193206910e-3,    1.3091347650e-3,
        2.0097999370e-4,    1.8240150530e-3,    3.4214719200e-3,   -6.9909263400e-3,
        -2.5461784100e-3,    6.3166203910e-3,   -1.2036296540e-4,    3.5721848250e-3,
        -2.7543718460e-3,   -1.3045557770e-2,    9.2808939520e-3,    9.8842596640e-3,
        -5.4998248820e-3,    1.8447549080e-3,   -1.4341374860e-2,   -4.8868563030e-3,
        3.0439257620e-2,   -3.3579096200e-3,   -1.6723809760e-2,   -2.6860507210e-4,
        -1.9729316230e-2,    3.0300367620e-2,    3.7809621540e-2,   -5.3283456710e-2,
        -1.6617689280e-2,    1.0394629090e-2,   -2.7511343360e-3,    1.0660593960e-1,
        -4.2705331000e-2,   -2.2114472090e-1,    1.5195469560e-1,    2.3963263630e-1,
    ];
    let expected_deviations = vec![
        1.2273662830e-6,    2.4547325670e-6,    1.2273662830e-6,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    1.1718750000e-2,    1.9531250000e-2,    2.7343750000e-2,
        3.5156250000e-2,    4.2968750000e-2,    5.0781250000e-2,    5.8593750000e-2,
        6.6406250000e-2,    7.4218750000e-2,    8.2031250000e-2,    8.9843750000e-2,
        9.7656250000e-2,    1.0546875000e-1,    1.1328125000e-1,    1.2109375000e-1,
        1.2890625000e-1,    1.3281250000e-1,    1.4062500000e-1,    1.4843750000e-1,
        1.5234375000e-1,    1.5625000000e-1,    1.6015625000e-1,    1.6406250000e-1,
        1.7000000180e-1,    2.1999999880e-1,    2.2390624880e-1,    2.2781249880e-1,
        2.3171874880e-1,    2.3562499880e-1,    2.3953124880e-1,    2.4734374880e-1,
        2.5515624880e-1,    2.5906249880e-1,    2.6687499880e-1,    2.7468749880e-1,
        2.8249999880e-1,    2.9031249880e-1,    2.9421874880e-1,    3.0203124880e-1,
        3.0984374880e-1,    3.1374999880e-1,    3.1765624880e-1,    3.2156249880e-1,
        3.2546874880e-1,    3.3000001310e-1,    3.7999999520e-1,    3.8390624520e-1,
        3.8781249520e-1,    3.9171874520e-1,    3.9562499520e-1,    3.9953124520e-1,
        4.0734374520e-1,    4.1515624520e-1,    4.1906249520e-1,    4.2687499520e-1,
        4.3468749520e-1,    4.4249999520e-1,    4.5031249520e-1,    4.5812499520e-1,
        4.6593749520e-1,    4.7374999520e-1,    4.8156249520e-1,    4.8937499520e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,3,3,3
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n128_grid3() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 3);

    let expected_impulse_response = vec![
        -1.5841082760e-7,   -5.6038620640e-6,    9.1733545560e-6,    1.9307104590e-5,
        -2.8690123140e-5,   -1.5307326980e-5,    1.7154779930e-5,   -3.9728256520e-6,
        7.5758296590e-5,   -2.0069583120e-5,   -1.7941379340e-4,    8.7809363320e-5,
        1.1068968160e-4,   -1.0447183740e-6,    1.2528189110e-4,   -3.8481166120e-4,
        -2.0398823840e-4,    6.7633565050e-4,    8.3504128270e-6,   -1.8287703280e-4,
        -1.0400078460e-4,   -9.3086704150e-4,    8.9956750160e-4,    1.3447938250e-3,
        -1.2962019540e-3,   -4.0262413680e-4,   -4.0532046110e-4,   -3.7975900340e-4,
        3.3761975360e-3,   -5.0697918050e-4,   -3.8321912290e-3,    1.0856395820e-3,
        2.0986888560e-5,    2.1374812350e-3,    3.8930103180e-3,   -7.2348145770e-3,
        -2.9788750690e-3,    6.1666523110e-3,    7.2712602560e-5,    4.2069125920e-3,
        -2.7149883100e-3,   -1.3766063380e-2,    8.9951027180e-3,    1.0140083730e-2,
        -4.9058170990e-3,    2.2078452170e-3,   -1.5048218890e-2,   -5.5734170600e-3,
        3.0737098310e-2,   -2.6932405310e-3,   -1.6251198950e-2,   -7.1772234510e-4,
        -2.0714053880e-2,    3.0336301770e-2,    3.8651447740e-2,   -5.2714716640e-2,
        -1.6871009020e-2,    9.3839317560e-3,   -3.1029116360e-3,    1.0744354130e-1,
        -4.1931528600e-2,   -2.2124224900e-1,    1.5102426710e-1,    2.3898908500e-1,
    ];
    let expected_deviations = vec![
        2.7459884680e-6,    5.4919769350e-6,    2.7459884680e-6,
    ];
    let expected_extremal_frequencies = vec![
        5.2083334890e-3,    1.3020833950e-2,    2.0833332090e-2,    2.8645830230e-2,
        3.6458332090e-2,    4.4270835820e-2,    5.2083339540e-2,    5.9895843270e-2,
        6.7708335820e-2,    7.5520828370e-2,    8.3333320920e-2,    9.1145813470e-2,
        9.8958306010e-2,    1.0677079860e-1,    1.1458329110e-1,    1.2239578370e-1,
        1.3020828370e-1,    1.3541662690e-1,    1.4322914180e-1,    1.4843748510e-1,
        1.5625000000e-1,    1.6145834330e-1,    1.6406251490e-1,    1.6666668650e-1,
        1.7000000180e-1,    2.1999999880e-1,    2.2260417040e-1,    2.2520834210e-1,
        2.2781251370e-1,    2.3562502860e-1,    2.4083337190e-1,    2.4604171510e-1,
        2.5385421510e-1,    2.6166668530e-1,    2.6687499880e-1,    2.7468746900e-1,
        2.8249993920e-1,    2.8770825270e-1,    2.9552072290e-1,    3.0333319310e-1,
        3.0854150650e-1,    3.1635397670e-1,    3.2156229020e-1,    3.2416644690e-1,
        3.2677060370e-1,    3.3000001310e-1,    3.7999999520e-1,    3.8260415200e-1,
        3.8520830870e-1,    3.8781246540e-1,    3.9302077890e-1,    4.0083324910e-1,
        4.0604156260e-1,    4.1385403280e-1,    4.1906234620e-1,    4.2687481640e-1,
        4.3468728660e-1,    4.4249975680e-1,    4.5031222700e-1,    4.5812469720e-1,
        4.6593716740e-1,    4.7374963760e-1,    4.8156210780e-1,    4.9197873470e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
128,3,3,15
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n128_grid15() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 15);

    let expected_impulse_response = vec![
        -1.1132793530e-7,   -6.8347617340e-6,    1.0110943550e-5,    2.4950997610e-5,
        -3.2240779550e-5,   -2.0817635230e-5,    1.9336468540e-5,   -5.4961001300e-6,
        8.6639047370e-5,   -1.2742559190e-5,   -2.0604337620e-4,    8.4813640570e-5,
        1.2836650420e-4,    9.6736039270e-7,    1.4828270650e-4,   -4.0906301000e-4,
        -2.5861160250e-4,    7.2288385130e-4,    4.3657171770e-5,   -1.9258180690e-4,
        -9.5561845230e-5,   -1.0187245210e-3,    8.8690419220e-4,    1.4854249310e-3,
        -1.3076614120e-3,   -4.5646593210e-4,   -4.3341412670e-4,   -4.8275582960e-4,
        3.5179196860e-3,   -3.4946727100e-4,   -4.0028025400e-3,    1.0146219280e-3,
        1.9776634870e-6,    2.1364749410e-3,    4.1802776980e-3,   -7.2764996440e-3,
        -3.3060407730e-3,    6.2211556360e-3,    1.3423955530e-4,    4.3409815990e-3,
        -2.4933828970e-3,   -1.4143617820e-2,    8.7545150890e-3,    1.0452914050e-2,
        -4.8137139530e-3,    2.3391996510e-3,   -1.5114245940e-2,   -6.1077419670e-3,
        3.0901581050e-2,   -2.2311815990e-3,   -1.6325037930e-2,   -7.3779956440e-4,
        -2.1027643230e-2,    3.0042864380e-2,    3.9255682380e-2,   -5.2470654250e-2,
        -1.7229124900e-2,    9.2775132510e-3,   -3.3807922150e-3,    1.0761545600e-1,
        -4.1262242940e-2,   -2.2151310740e-1,    1.5056943890e-1,    2.3902776840e-1,
    ];
    let expected_deviations = vec![
        4.3256054600e-6,    8.6512109190e-6,    4.3256054600e-6,
    ];
    let expected_extremal_frequencies = vec![
        4.1666668840e-3,    1.1979170140e-2,    2.0312501120e-2,    2.8645826500e-2,
        3.6458320920e-2,    4.4791646300e-2,    5.2604138850e-2,    6.0416631400e-2,
        6.8749956790e-2,    7.6562449340e-2,    8.4374941890e-2,    9.2187434430e-2,
        9.9999926980e-2,    1.0729158670e-1,    1.1510407920e-1,    1.2239573900e-1,
        1.2968745830e-1,    1.3697922230e-1,    1.4375014600e-1,    1.5052106980e-1,
        1.5625031290e-1,    1.6197955610e-1,    1.6614627840e-1,    1.6875047980e-1,
        1.7000000180e-1,    2.1999999880e-1,    2.2104167940e-1,    2.2416672110e-1,
        2.2833344340e-1,    2.3406268660e-1,    2.4031277000e-1,    2.4708369370e-1,
        2.5385451320e-1,    2.6062524320e-1,    2.6791679860e-1,    2.7520835400e-1,
        2.8197908400e-1,    2.8927063940e-1,    2.9604136940e-1,    3.0333292480e-1,
        3.0958282950e-1,    3.1583273410e-1,    3.2156181340e-1,    3.2572841640e-1,
        3.2885336880e-1,    3.3000001310e-1,    3.7999999520e-1,    3.8104164600e-1,
        3.8364577290e-1,    3.8833320140e-1,    3.9354145530e-1,    3.9979135990e-1,
        4.0604126450e-1,    4.1333281990e-1,    4.2010355000e-1,    4.2739510540e-1,
        4.3520748620e-1,    4.4301986690e-1,    4.5083224770e-1,    4.5864462850e-1,
        4.6645700930e-1,    4.7479021550e-1,    4.8312342170e-1,    4.9145662780e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}

/*
(CONVERGE FAILURE)
128,3,3,16
0.000000000000,0.170000001788,0.219999998808,0.330000013113
0.379999995232,0.500000000000
0.000000000000,1.000000000000,0.000000000000
2.000000000000,1.000000000000,2.000000000000
*/
#[test]
fn hilbert_bandpass_n128_grid16() {
    let mut bands = vec![];
    bands.push(Band {
        lower_edge: 0f32,
        upper_edge: 0.17f32,
        desired_value: 0f32,
        weight: 2f32,
    });
    bands.push(Band {
        lower_edge: 0.22f32,
        upper_edge: 0.33f32,
        desired_value: 1f32,
        weight: 1f32,
    });
    bands.push(Band {
        lower_edge: 0.38f32,
        upper_edge: 0.5f32,
        desired_value: 0f32,
        weight: 2f32,
    });

    let pm_output = design_filter(128, FilterType::HilbertTransform, &bands, 16);

    let expected_impulse_response = vec![
        -1.2289892480e-7,   -6.8939666560e-6,    1.0129610020e-5,    2.5075252780e-5,
        -3.2207783080e-5,   -2.0904073610e-5,    1.9123321180e-5,   -5.5519954000e-6,
        8.7083295510e-5,   -1.2482305460e-5,   -2.0642972960e-4,    8.4297498690e-5,
        1.2815881930e-4,    1.7803686210e-6,    1.4932148040e-4,   -4.0981310300e-4,
        -2.6009290010e-4,    7.2259386070e-4,    4.4915504990e-5,   -1.9042336500e-4,
        -9.6039228080e-5,   -1.0221141860e-3,    8.8589510410e-4,    1.4878872320e-3,
        -1.3043222720e-3,   -4.5605300690e-4,   -4.3867956260e-4,   -4.8629572850e-4,
        3.5221176220e-3,   -3.4381658770e-4,   -4.0018754080e-3,    1.0083147790e-3,
        -5.2941031750e-6,    2.1408488970e-3,    4.1904994290e-3,   -7.2748456150e-3,
        -3.3139348960e-3,    6.2109953720e-3,    1.3625115390e-4,    4.3560666960e-3,
        -2.4880804120e-3,   -1.4154649340e-2,    8.7417801840e-3,    1.0452475400e-2,
        -4.7966912390e-3,    2.3511948530e-3,   -1.5126999470e-2,   -6.1255227770e-3,
        3.0900232490e-2,   -2.2143479440e-3,   -1.6307525340e-2,   -7.4744061570e-4,
        -2.1052081140e-2,    3.0039099980e-2,    3.9273746310e-2,   -5.2451215680e-2,
        -1.7233062540e-2,    9.2500504110e-3,   -3.3911746000e-3,    1.0763537880e-1,
        -4.1241250930e-2,   -2.2151370350e-1,    1.5054461360e-1,    2.3901018500e-1,
    ];
    let expected_deviations = vec![
        4.3289132920e-6,    8.6578265840e-6,    4.3289132920e-6,
    ];
    let expected_extremal_frequencies = vec![
        3.9062500000e-3,    1.2207031250e-2,    2.0507812500e-2,    2.8320312500e-2,
        3.6621093750e-2,    4.4433593750e-2,    5.2734375000e-2,    6.0546875000e-2,
        6.8359375000e-2,    7.6660156250e-2,    8.4472656250e-2,    9.2285156250e-2,
        9.9609375000e-2,    1.0742187500e-1,    1.1523437500e-1,    1.2255859380e-1,
        1.2988281250e-1,    1.3720703120e-1,    1.4404296880e-1,    1.5039062500e-1,
        1.5625000000e-1,    1.6162109380e-1,    1.6601562500e-1,    1.6894531250e-1,
        1.7000000180e-1,    2.1999999880e-1,    2.2097656130e-1,    2.2390624880e-1,
        2.2878906130e-1,    2.3416015510e-1,    2.4050781130e-1,    2.4685546760e-1,
        2.5369140510e-1,    2.6052734260e-1,    2.6785156130e-1,    2.7517578010e-1,
        2.8201171760e-1,    2.8933593630e-1,    2.9617187380e-1,    3.0300781130e-1,
        3.0984374880e-1,    3.1570312380e-1,    3.2156249880e-1,    3.2595703010e-1,
        3.2888671760e-1,    3.3000001310e-1,    3.7999999520e-1,    3.8097655770e-1,
        3.8390624520e-1,    3.8830077650e-1,    3.9367187020e-1,    3.9953124520e-1,
        4.0636718270e-1,    4.1320312020e-1,    4.2003905770e-1,    4.2785155770e-1,
        4.3517577650e-1,    4.4298827650e-1,    4.5080077650e-1,    4.5861327650e-1,
        4.6642577650e-1,    4.7472655770e-1,    4.8302733900e-1,    4.9132812020e-1,
        5.0000000000e-1,
    ];

    test_inputs(&pm_output, &bands);
    test_impulse_response(&pm_output, &expected_impulse_response);
    test_deviations(&pm_output, &expected_deviations);
    test_extremal_frequencies(&pm_output, &expected_extremal_frequencies);
}