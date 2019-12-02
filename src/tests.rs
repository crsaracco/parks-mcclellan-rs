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
const FLOATING_POINT_TOLERANCE: f64 = 0.00001; // 0.001%

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
        test_f64_error_within_tolerance(output.impulse_response[0], -0.5753406789E-02);
        test_f64_error_within_tolerance(output.impulse_response[1], 0.9902641177E-03);
        test_f64_error_within_tolerance(output.impulse_response[2], 0.7573357783E-02);
        test_f64_error_within_tolerance(output.impulse_response[3], -0.6514116656E-02);
        test_f64_error_within_tolerance(output.impulse_response[4], 0.1396050397E-01);
        test_f64_error_within_tolerance(output.impulse_response[5], 0.2295166254E-02);
        test_f64_error_within_tolerance(output.impulse_response[6], -0.1999405399E-01);
        test_f64_error_within_tolerance(output.impulse_response[7], 0.7136960514E-02);
        test_f64_error_within_tolerance(output.impulse_response[8], -0.3965735435E-01);
        test_f64_error_within_tolerance(output.impulse_response[9], 0.1126007363E-01);
        test_f64_error_within_tolerance(output.impulse_response[10], 0.6623364240E-01);
        test_f64_error_within_tolerance(output.impulse_response[11], -0.1049720775E-01);
        test_f64_error_within_tolerance(output.impulse_response[12], 0.8513614535E-01);
        test_f64_error_within_tolerance(output.impulse_response[13], -0.1202498823E+00);
        test_f64_error_within_tolerance(output.impulse_response[14], -0.2967857718E+00);
        test_f64_error_within_tolerance(output.impulse_response[15], 0.3041091561E+00);
    }

    #[test]
    fn test_multiband_01_deviations() {
        let output = multiband_01();
        test_f64_error_within_tolerance(output.deviations[0], 0.00151311723);
        test_f64_error_within_tolerance(output.deviations[1], 0.01513117179);
        test_f64_error_within_tolerance(output.deviations[2], 0.00151311723);
    }

    #[test]
    fn test_multiband_01_deviation_dbss() {
        let output = multiband_01();
        let deviation_dbs = output.deviation_dbs.unwrap();
        test_f64_error_within_tolerance(deviation_dbs[0], -56.40254974365);
        test_f64_error_within_tolerance(deviation_dbs[1], 0.13044279814);
        test_f64_error_within_tolerance(deviation_dbs[2], -56.40254974365);
    }

    #[test]
    fn test_multiband_01_extremal_frequencies() {
        let output = multiband_01();
        test_f64_error_within_tolerance(output.extremal_frequencies[0], 0.0000000);
        test_f64_error_within_tolerance(output.extremal_frequencies[1], 0.0273438);
        test_f64_error_within_tolerance(output.extremal_frequencies[2], 0.0527344);
        test_f64_error_within_tolerance(output.extremal_frequencies[3], 0.0761719);
        test_f64_error_within_tolerance(output.extremal_frequencies[4], 0.0937500);
        test_f64_error_within_tolerance(output.extremal_frequencies[5], 0.1000000);
        test_f64_error_within_tolerance(output.extremal_frequencies[6], 0.2000000);
        test_f64_error_within_tolerance(output.extremal_frequencies[7], 0.2195313);
        test_f64_error_within_tolerance(output.extremal_frequencies[8], 0.2527344);
        test_f64_error_within_tolerance(output.extremal_frequencies[9], 0.2839844);
        test_f64_error_within_tolerance(output.extremal_frequencies[10], 0.3132812);
        test_f64_error_within_tolerance(output.extremal_frequencies[11], 0.3386719);
        test_f64_error_within_tolerance(output.extremal_frequencies[12], 0.3500000);
        test_f64_error_within_tolerance(output.extremal_frequencies[13], 0.4250000);
        test_f64_error_within_tolerance(output.extremal_frequencies[14], 0.4328125);
        test_f64_error_within_tolerance(output.extremal_frequencies[15], 0.4503906);
        test_f64_error_within_tolerance(output.extremal_frequencies[16], 0.4796875);
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
        test_f64_error_within_tolerance(output.impulse_response[0], -0.6257381756E-03);
        test_f64_error_within_tolerance(output.impulse_response[1], 0.8548910846E-03);
        test_f64_error_within_tolerance(output.impulse_response[2], -0.4238411202E-03);
        test_f64_error_within_tolerance(output.impulse_response[3], 0.3989883116E-03);
        test_f64_error_within_tolerance(output.impulse_response[4], -0.4351424286E-03);
        test_f64_error_within_tolerance(output.impulse_response[5], 0.5008199369E-03);
        test_f64_error_within_tolerance(output.impulse_response[6], -0.5971363280E-03);
        test_f64_error_within_tolerance(output.impulse_response[7], 0.7333782269E-03);
        test_f64_error_within_tolerance(output.impulse_response[8], -0.9301801911E-03);
        test_f64_error_within_tolerance(output.impulse_response[9], 0.1226891065E-02);
        test_f64_error_within_tolerance(output.impulse_response[10], -0.1701220521E-02);
        test_f64_error_within_tolerance(output.impulse_response[11], 0.2527206670E-02);
        test_f64_error_within_tolerance(output.impulse_response[12], -0.4160370212E-02);
        test_f64_error_within_tolerance(output.impulse_response[13], 0.8130028844E-02);
        test_f64_error_within_tolerance(output.impulse_response[14], -0.2253980562E-01);
        test_f64_error_within_tolerance(output.impulse_response[15], 0.2026662678E+00);
    }

    #[test]
    fn test_differentiator_01_deviations() {
        let output = differentiator_01();
        test_f64_error_within_tolerance(output.deviations[0], 0.00619235216);
    }

    #[test]
    fn test_differentiator_01_extremal_frequencies() {
        let output = differentiator_01();
        test_f64_error_within_tolerance(output.extremal_frequencies[0], 0.0015625);
        test_f64_error_within_tolerance(output.extremal_frequencies[1], 0.0328125);
        test_f64_error_within_tolerance(output.extremal_frequencies[2], 0.0656250);
        test_f64_error_within_tolerance(output.extremal_frequencies[3], 0.0984374);
        test_f64_error_within_tolerance(output.extremal_frequencies[4], 0.1312499);
        test_f64_error_within_tolerance(output.extremal_frequencies[5], 0.1656251);
        test_f64_error_within_tolerance(output.extremal_frequencies[6], 0.1984377);
        test_f64_error_within_tolerance(output.extremal_frequencies[7], 0.2312503);
        test_f64_error_within_tolerance(output.extremal_frequencies[8], 0.2640629);
        test_f64_error_within_tolerance(output.extremal_frequencies[9], 0.2968756);
        test_f64_error_within_tolerance(output.extremal_frequencies[10], 0.3296882);
        test_f64_error_within_tolerance(output.extremal_frequencies[11], 0.3625008);
        test_f64_error_within_tolerance(output.extremal_frequencies[12], 0.3953134);
        test_f64_error_within_tolerance(output.extremal_frequencies[13], 0.4281261);
        test_f64_error_within_tolerance(output.extremal_frequencies[14], 0.4593762);
        test_f64_error_within_tolerance(output.extremal_frequencies[15], 0.4859388);
        test_f64_error_within_tolerance(output.extremal_frequencies[16], 0.5000000);
    }
}