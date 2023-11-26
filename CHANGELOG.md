# Changelog

This file tracks the notable changes to the 2decomp&fft project (starting from v2.0).

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

When a release is made the "Unreleased" section should be renamed to that version and date
_e.g._ vX.Y - YYYY-MM-DD and a new "Unreleased" section started above.

## Unreleased

### Added

- Add profiling to transpose subroutines when `p_row` or `p_col` is equal to one. See [PR #251](https://github.com/2decomp-fft/2decomp-fft/pull/251)
- Add multigrid capability for the FFT (experimental support limited to fftw3_f03 currently). Based on the work of [Lionel Gelebart](https://github.com/LionelGelebart) and coworkers, see [AMITEX_FFTP](https://amitexfftp.github.io/AMITEX/index.html) and [2decomp15_mg](https://github.com/LionelGelebart/2decomp15_mg).

### Fixed

- Fuse transpose CPU and GPU memory buffers to reduce memory usage. See [PR 271](https://github.com/2decomp-fft/2decomp-fft/pull/271)

### Changed
### Deprecated
### Removed

## v2.0.1

- See JOSS paper (under submission)

### Added
### Fixed
### Changed
### Deprecated
### Removed

## v2.0 - 2023-08-30
