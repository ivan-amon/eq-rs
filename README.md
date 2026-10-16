# eq-rs

An audio equalizer written in Rust, built from biquad filters.

By Iván Amón.

## Status

Work in progress. The Rust side does offline filtering: it reads a WAV file, applies the filters
and writes a new WAV file. Real-time processing is not supported for the moment.

## Repository

- [`src/`](src/): the Rust code.
- [`docs/biquads.md`](docs/biquads.md): the math of every filter, and the notation used across the repo.
- [`analysis/`](analysis/): Python scripts that check the filter formulas and plot them in the
  s-plane, the z-plane and the frequency domain. See its [README](analysis/README.md).
