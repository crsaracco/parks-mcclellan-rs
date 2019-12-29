# parks-mcclellan-rs

A Rust implementation of the [Parks-McClellan filter design algorithm](https://en.wikipedia.org/wiki/Parks%E2%80%93McClellan_filter_design_algorithm) for finding optimal Chebyshev FIR filters.

## TODO

 - [ ] Finish refactoring everything
 - [ ] Make separate `design_filter` functions for the different types of filter, and hide `FilterType`
 - [ ] Remove everything in ParksMcClellanOutput that's unnecessary (maybe just leave the impulse response?)
