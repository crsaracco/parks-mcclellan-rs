use super::*;

// Since some of the calculations in the original Parks-McClellan algorithm
// use floats (f32s), we want to check against the number of sig-figs of
// floats, not doubles. Floats give sig-figs up to 0.0001%.
// See: https://stackoverflow.com/a/12815366/6622465
// However, there's probably some more math rounding errors and such inherent
// in the algorithm, so we'll remove one sig-fig (to 0.001%). If we hit this
// in all of the tests, it's pretty likely that *we're* more precise than the
// Fortran code because we're using f64s everywhere.
// ... Although it could also be that we're doing something slightly wrong
// somewhere, so it might be worth double-checking eventually.
const FLOATING_POINT_TOLERANCE: f64 = std::f32::EPSILON as f64; // roughly 0.000011920929 %

#[cfg(test)]
mod tests {
    use super::*;

    fn test_f64_error_within_tolerance(measured: f64, expected: f64) {
        if expected == 0.0 {
            assert!(measured.abs() <= FLOATING_POINT_TOLERANCE);
            return;
        }

        let error = (measured - expected) / expected;
        assert!(
            error.abs() <= FLOATING_POINT_TOLERANCE,
            format!("measured: {}, expected: {}, error: {}", measured, expected, error)
        );
    }

    // This is the MultipleBand example from the fortran code's docs.
    fn multiband_01() -> ParksMcclellanOutput {
        let mut bands = vec![];
        bands.push(Band{
            lower_edge: 0.0,
            upper_edge: 0.1,
            desired_value: 0.0,
            weight: 10.0,
        });
        bands.push(Band{
            lower_edge: 0.2,
            upper_edge: 0.35,
            desired_value: 1.0,
            weight: 1.0,
        });
        bands.push(Band{
            lower_edge: 0.425,
            upper_edge: 0.5,
            desired_value: 0.0,
            weight: 10.0,
        });
        design(32, JType::MultipleBand, &bands, 16)
    }

    // This is the Differentiator example from the fortran code's docs.
    fn differentiator_01() -> ParksMcclellanOutput {
        let mut bands = vec![];
        bands.push(Band{
            lower_edge: 0.0,
            upper_edge: 0.5,
            desired_value: 1.0,
            weight: 1.0,
        });
        design(32, JType::Differentiator, &bands, 20)
    }

    #[test]
    fn test_multiband_01_inputs() {
        let output = multiband_01();

        assert_eq!(output.filter_length, 32);

        // These floats should all be *exact*, since they're our inputs.
        assert_eq!(output.lower_band_edges[0], 0.0);
        assert_eq!(output.lower_band_edges[1], 0.2);
        assert_eq!(output.lower_band_edges[2], 0.425);

        assert_eq!(output.upper_band_edges[0], 0.1);
        assert_eq!(output.upper_band_edges[1], 0.35);
        assert_eq!(output.upper_band_edges[2], 0.5);

        assert_eq!(output.desired_values[0], 0.0);
        assert_eq!(output.desired_values[1], 1.0);
        assert_eq!(output.desired_values[2], 0.0);

        assert_eq!(output.weightings[0], 10.0);
        assert_eq!(output.weightings[1], 1.0);
        assert_eq!(output.weightings[2], 10.0);
    }

    #[test]
    fn test_multiband_01_impulse_response() {
        let output = multiband_01();
        test_f64_error_within_tolerance(output.impulse_response[0].into(), -0.5753406789E-02);
        test_f64_error_within_tolerance(output.impulse_response[1].into(), 0.9902641177E-03);
        test_f64_error_within_tolerance(output.impulse_response[2].into(), 0.7573357783E-02);
        test_f64_error_within_tolerance(output.impulse_response[3].into(), -0.6514116656E-02);
        test_f64_error_within_tolerance(output.impulse_response[4].into(), 0.1396050397E-01);
        test_f64_error_within_tolerance(output.impulse_response[5].into(), 0.2295166254E-02);
        test_f64_error_within_tolerance(output.impulse_response[6].into(), -0.1999405399E-01);
        test_f64_error_within_tolerance(output.impulse_response[7].into(), 0.7136960514E-02);
        test_f64_error_within_tolerance(output.impulse_response[8].into(), -0.3965735435E-01);
        test_f64_error_within_tolerance(output.impulse_response[9].into(), 0.1126007363E-01);
        test_f64_error_within_tolerance(output.impulse_response[10].into(), 0.6623364240E-01);
        test_f64_error_within_tolerance(output.impulse_response[11].into(), -0.1049720775E-01);
        test_f64_error_within_tolerance(output.impulse_response[12].into(), 0.8513614535E-01);
        test_f64_error_within_tolerance(output.impulse_response[13].into(), -0.1202498823E+00);
        test_f64_error_within_tolerance(output.impulse_response[14].into(), -0.2967857718E+00);
        test_f64_error_within_tolerance(output.impulse_response[15].into(), 0.3041091561E+00);
    }

    #[test]
    fn test_multiband_01_deviations() {
        let output = multiband_01();
        test_f64_error_within_tolerance(output.deviations[0].into(), 0.00151311723);
        test_f64_error_within_tolerance(output.deviations[1].into(), 0.01513117179);
        test_f64_error_within_tolerance(output.deviations[2].into(), 0.00151311723);
    }

    #[test]
    fn test_multiband_01_deviation_dbs() {
        let output = multiband_01();
        let deviation_dbs = output.deviation_dbs.unwrap();
        test_f64_error_within_tolerance(deviation_dbs[0].into(), -56.40254974365);
        test_f64_error_within_tolerance(deviation_dbs[1].into(), 0.13044279814);
        test_f64_error_within_tolerance(deviation_dbs[2].into(), -56.40254974365);
    }

    #[test]
    fn test_multiband_01_extremal_frequencies() {
        let output = multiband_01();
        test_f64_error_within_tolerance(output.extremal_frequencies[0].into(), 0.000000000);
        test_f64_error_within_tolerance(output.extremal_frequencies[1].into(), 0.027343750);
        test_f64_error_within_tolerance(output.extremal_frequencies[2].into(), 0.052734375);
        test_f64_error_within_tolerance(output.extremal_frequencies[3].into(), 0.076171875);
        test_f64_error_within_tolerance(output.extremal_frequencies[4].into(), 0.093750000);
        test_f64_error_within_tolerance(output.extremal_frequencies[5].into(), 0.100000001);
        test_f64_error_within_tolerance(output.extremal_frequencies[6].into(), 0.200000003);
        test_f64_error_within_tolerance(output.extremal_frequencies[7].into(), 0.219531253);
        test_f64_error_within_tolerance(output.extremal_frequencies[8].into(), 0.252734363);
        test_f64_error_within_tolerance(output.extremal_frequencies[9].into(), 0.283984363);
        test_f64_error_within_tolerance(output.extremal_frequencies[10].into(), 0.313281238);
        test_f64_error_within_tolerance(output.extremal_frequencies[11].into(), 0.338671863);
        test_f64_error_within_tolerance(output.extremal_frequencies[12].into(), 0.349999994);
        test_f64_error_within_tolerance(output.extremal_frequencies[13].into(), 0.425000012);
        test_f64_error_within_tolerance(output.extremal_frequencies[14].into(), 0.432812512);
        test_f64_error_within_tolerance(output.extremal_frequencies[15].into(), 0.450390637);
        test_f64_error_within_tolerance(output.extremal_frequencies[16].into(), 0.479687512);
    }

    #[test]
    fn test_differentiator_01_inputs() {
        let output = differentiator_01();

        assert_eq!(output.filter_length, 32);

        // These floats should all be *exact*, since they're our inputs.
        assert_eq!(output.lower_band_edges[0], 0.0);
        assert_eq!(output.upper_band_edges[0], 0.5);
        assert_eq!(output.desired_values[0], 1.0);
        assert_eq!(output.weightings[0], 1.0);
    }

    #[test]
    fn test_differentiator_01_impulse_response() {
        let output = differentiator_01();
        test_f64_error_within_tolerance(output.impulse_response[0].into(), -0.6257381756E-03);
        test_f64_error_within_tolerance(output.impulse_response[1].into(), 0.8548910846E-03);
        test_f64_error_within_tolerance(output.impulse_response[2].into(), -0.4238411202E-03);
        test_f64_error_within_tolerance(output.impulse_response[3].into(), 0.3989883116E-03);
        test_f64_error_within_tolerance(output.impulse_response[4].into(), -0.4351424286E-03);
        test_f64_error_within_tolerance(output.impulse_response[5].into(), 0.5008199369E-03);
        test_f64_error_within_tolerance(output.impulse_response[6].into(), -0.5971363280E-03);
        test_f64_error_within_tolerance(output.impulse_response[7].into(), 0.7333782269E-03);
        test_f64_error_within_tolerance(output.impulse_response[8].into(), -0.9301801911E-03);
        test_f64_error_within_tolerance(output.impulse_response[9].into(), 0.1226891065E-02);
        test_f64_error_within_tolerance(output.impulse_response[10].into(), -0.1701220521E-02);
        test_f64_error_within_tolerance(output.impulse_response[11].into(), 0.2527206670E-02);
        test_f64_error_within_tolerance(output.impulse_response[12].into(), -0.4160370212E-02);
        test_f64_error_within_tolerance(output.impulse_response[13].into(), 0.8130028844E-02);
        test_f64_error_within_tolerance(output.impulse_response[14].into(), -0.2253980562E-01);
        test_f64_error_within_tolerance(output.impulse_response[15].into(), 0.2026662678E+00);
    }

    #[test]
    fn test_differentiator_01_deviations() {
        let output = differentiator_01();
        test_f64_error_within_tolerance(output.deviations[0].into(), 0.00619235216);
    }

    #[test]
    fn test_differentiator_01_extremal_frequencies() {
        let output = differentiator_01();
        test_f64_error_within_tolerance(output.extremal_frequencies[0].into(), 0.001562500);
        test_f64_error_within_tolerance(output.extremal_frequencies[1].into(), 0.032812502);
        test_f64_error_within_tolerance(output.extremal_frequencies[2].into(), 0.065624975);
        test_f64_error_within_tolerance(output.extremal_frequencies[3].into(), 0.098437443);
        test_f64_error_within_tolerance(output.extremal_frequencies[4].into(), 0.131249934);
        test_f64_error_within_tolerance(output.extremal_frequencies[5].into(), 0.165625066);
        test_f64_error_within_tolerance(output.extremal_frequencies[6].into(), 0.198437691);
        test_f64_error_within_tolerance(output.extremal_frequencies[7].into(), 0.231250316);
        test_f64_error_within_tolerance(output.extremal_frequencies[8].into(), 0.264062941);
        test_f64_error_within_tolerance(output.extremal_frequencies[9].into(), 0.296875566);
        test_f64_error_within_tolerance(output.extremal_frequencies[10].into(), 0.329688191);
        test_f64_error_within_tolerance(output.extremal_frequencies[11].into(), 0.362500817);
        test_f64_error_within_tolerance(output.extremal_frequencies[12].into(), 0.395313442);
        test_f64_error_within_tolerance(output.extremal_frequencies[13].into(), 0.428126067);
        test_f64_error_within_tolerance(output.extremal_frequencies[14].into(), 0.459376186);
        test_f64_error_within_tolerance(output.extremal_frequencies[15].into(), 0.485938787);
        test_f64_error_within_tolerance(output.extremal_frequencies[16].into(), 0.500000000);
    }
}