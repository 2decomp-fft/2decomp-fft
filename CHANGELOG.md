# Changelog

This file tracks the notable changes to the 2decomp&fft project (starting from v2.0).

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

When a release is made the "Unreleased" section should be renamed to that version and date
_e.g._ vX.Y - YYYY-MM-DD and a new "Unreleased" section started above.

## Unreleased

### Added

- Add the possibility to skip c2c transforms along one or several dimensions [PR #339](https://github.com/2decomp-fft/2decomp-fft/pull/339), [PR #340](https://github.com/2decomp-fft/2decomp-fft/pull/340)
- Add profiling to transpose subroutines when `p_row` or `p_col` is equal to one. See [PR #251](https://github.com/2decomp-fft/2decomp-fft/pull/251)
- The FFT can be performed on various grids for all backends. In addition, the fftw_f03 FFT backend also supports in-place r2c and in-place c2r transforms (experimental support currently). Based on the work of [Lionel Gelebart](https://github.com/LionelGelebart) and coworkers, see [AMITEX_FFTP](https://amitexfftp.github.io/AMITEX/index.html) and [2decomp15_mg](https://github.com/LionelGelebart/2decomp15_mg). See [PR 302](https://github.com/2decomp-fft/2decomp-fft/pull/302).
- The external code can now provide [hints for MPI IO](https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node315.htm#Node316). Non-blocking MPI IO is also added. See [PR #269](https://github.com/2decomp-fft/2decomp-fft/pull/269)
- Add the possibility to allocate 3D integers and logical arrays. See [PR 345](https://github.com/2decomp-fft/2decomp-fft/pull/345)

### Fixed

- Fix 3D complex-to-complex FFT for the MKL backend. See [PR 341](https://github.com/2decomp-fft/2decomp-fft/pull/341)
- Fuse transpose CPU and GPU memory buffers to reduce memory usage. See [PR 271](https://github.com/2decomp-fft/2decomp-fft/pull/271)

### Changed

- Major update in the IO routines. The new IO interface is not backward compatible. The examples are updated accordingly. See [PR 344](https://github.com/2decomp-fft/2decomp-fft/pull/344)
- The transpose subroutines have been moved into submodules rather than `#include` files. See [PR #206](https://github.com/2decomp-fft/2decomp-fft/pull/206)

### Deprecated

- Subroutines `init_coarser_mesh_stat*` and `fine_to_coarse*` will be removed. See [PR #316](https://github.com/2decomp-fft/2decomp-fft/pull/316).
- FFT backend fftw3 will be removed. See [PR #346](https://github.com/2decomp-fft/2decomp-fft/pull/346)

### Removed

- Removed unused `decomp_info` objects `phg`, `ph1`, `ph2`, `ph3` and `ph4`. See [PR #315](https://github.com/2decomp-fft/2decomp-fft/pull/315).
- Remove unused variable `real2_type`. See [PR #314](https://github.com/2decomp-fft/2decomp-fft/pull/314).

## v2.0.1

- See JOSS paper ([10.21105/joss.05813](https://doi.org/10.21105/joss.05813))

### Added
### Fixed
### Changed
### Deprecated
### Removed

## v2.0 - 2023-08-30
