mod pm;
use pm::{FilterType, design_filter};

/// A filter band.
pub struct Band {
    /// The lowest frequency ("edge") of the band.
    pub lower_edge: f32,

    /// The highest frequency ("edge") of the band.
    pub upper_edge: f32,

    /// The desired frequency response for this band.
    /// To make a pass band, set this to `1.0`.
    /// To make a stop band, set this to `0.0`.
    pub desired_value: f32,

    /// The "weight" of this band. The higher the weight (relative to other bands), the harder
    /// the algorithm tries to match this band in the filter's frequency response.
    pub weight: f32,
}

/// Design a multiple-band filter (low-pass, high-pass, band-pass, notch, etc).
///
///  - `filter_length` is the length of the impulse response of the resulting filter.
///    Must be in range [3, 128].
///  - `bands` holds information about the desired characteristics of the frequency response bands
///    for the filter. Must be in range [1, 10].
///  - `grid_density` is the "density" of the "dense grid". A higher number here will produce a more
///    accurate filter, at the expense of compute time. Must be greater than 0.
///
/// Example: pass-band filter with a pass-band of [0.2, 0.35] * sampling_frequency
/// ```
/// use parks_mcclellan::{Band, design_multiple_band_filter};
///
/// let mut bands = vec![];
/// // Stop frequencies between [0.0, 0.1]
/// bands.push(Band {
///     lower_edge: 0.0,
///     upper_edge: 0.1,
///     desired_value: 0.0,
///     weight: 10.0,
/// });
/// // Pass through frequencies between [0.2, 0.35]
/// // Note that [0.1, 0.2] is a transition band.
/// bands.push(Band {
///     lower_edge: 0.2,
///     upper_edge: 0.35,
///     desired_value: 1.0,
///     weight: 1.0,
/// });
/// // Stop frequencies between [0.425, 0.5]
/// // Note that [0.35, 0.425] is a transition band.
/// bands.push(Band {
///    lower_edge: 0.425,
///    upper_edge: 0.5,
///    desired_value: 0.0,
///    weight: 10.0,
/// });
///
/// let impulse_response = design_multiple_band_filter(32, &bands, 16);
/// ```
pub fn design_multiple_band_filter(
    filter_length: usize,
    bands: &[Band],
    grid_density: usize
) -> Vec<f32> {
    let pm_output = design_filter(filter_length, FilterType::MultipleBand, bands, grid_density);
    pm_output.impulse_response
}

/// Design a digital [Differentiator](https://en.wikipedia.org/wiki/Differentiator) filter.
///
///  - `filter_length` is the length of the impulse response of the resulting filter.
///    Must be in range [3, 128].
///  - `bands` holds information about the desired characteristics of the frequency response bands
///    for the filter. Must be in range [1, 10].
///  - `grid_density` is the "density" of the "dense grid". A higher number here will produce a more
///    accurate filter, at the expense of compute time. Must be greater than 0.
///
/// Example: full-band differentiator
/// ```
/// use parks_mcclellan::{Band, design_differentiator_filter};
///
/// let mut bands = vec![];
/// bands.push(Band {
///     lower_edge: 0.0,
///     upper_edge: 0.5,
///     desired_value: 1.0,
///     weight: 1.0,
/// });
///
/// let impulse_response = design_differentiator_filter(32, &bands, 20);
/// ```
pub fn design_differentiator_filter(
    filter_length: usize,
    bands: &[Band],
    grid_density: usize
) -> Vec<f32> {
    let pm_output = design_filter(filter_length, FilterType::Differentiator, bands, grid_density);
    pm_output.impulse_response
}

/// Design a [Hilbert Transform filter](https://ccrma.stanford.edu/~jos/mdft/Analytic_Signals_Hilbert_Transform.html).
///
///  - `filter_length` is the length of the impulse response of the resulting filter.
///    Must be in range [3, 128].
///  - `bands` holds information about the desired characteristics of the frequency response bands
///    for the filter. Must be in range [1, 10].
///  - `grid_density` is the "density" of the "dense grid". A higher number here will produce a more
///    accurate filter, at the expense of compute time. Must be greater than 0.
///
/// Example: full-band hilbert transformer
/// ```
/// use parks_mcclellan::{Band, design_hilbert_transform_filter};
///
/// let mut bands = vec![];
/// bands.push(Band {
///     lower_edge: 0.0,
///     upper_edge: 0.5,
///     desired_value: 1.0,
///     weight: 1.0,
/// });
///
/// let impulse_response = design_hilbert_transform_filter(32, &bands, 20);
/// ```
pub fn design_hilbert_transform_filter(
    filter_length: usize,
    bands: &[Band],
    grid_density: usize,
) -> Vec<f32> {
    let pm_output = design_filter(filter_length, FilterType::HilbertTransform, bands, grid_density);
    pm_output.impulse_response
}