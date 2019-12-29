# parks-mcclellan-rs

A Rust implementation of the [Parks-McClellan filter design algorithm](https://en.wikipedia.org/wiki/Parks%E2%80%93McClellan_filter_design_algorithm) for finding optimal Chebyshev FIR filters.

The construction of this filter is currently **not** suitable for use in a real-time DSP thread, so construct the filter before using it.

Includes more than 450 integration-level tests comparing the output to code from the canonical implementation. Not guaranteed to be 100% byte-for-byte exact to the canonical implementation for all inputs, but it should be pretty dang close.

## TODO

 - More refactoring
 - Remove the 128-sample limit on impulse response length
 - Remove the 10-band limit
 - Write a bunch of unit-level tests
 - Make the maximum number of iterations a parameter, or at least make it a higher limit
 - Write some "functional" tests (testing stop-band attenuation, for example)
 - Fix possibilities of dividing by zero
 - Write a 100% `f64` version
