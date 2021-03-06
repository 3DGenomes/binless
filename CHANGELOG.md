# Changelog
All notable changes to *binless* will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
for versions 0.x of binless, minor releases might break backwards compatibility.

## [0.15.1]
### Fixed
- bug arising when a particular locus does not have reads in at least one
  contact category

## [0.15.0]
### Added
- less errors regarding documentation when compiling
- New `chromosome_binless.R` script to normalize whole chromosomes
- New virtual 4C plot, example in `fast_binless.R`

## [0.14.0]
### Changed
- Tutorials now use the smaller SEMA3C dataset, which can be run quickly on a
  laptop
### Removed
- Old article/ folder

## [0.13.0]
### Added
- use 64 bit integers in gfl_graph_fl to expand the index limit in lasso.
- use long as StorageIndex to expand the limit of elements in Eigen SparseMatrix in fast binless.
- possibility to dump fast binless output to csv file to limit the memory used in the creation of DataFrame.
  To reduce the file size the csv will be compressed in gzip format using BOOST libraries if USE_BOOST is     activated in Makevars file.
- possibility to prevent computation of patch numbers in fast binless. This can
  be slow in very large datasets.
- maxdiag argument in fast_binless allows to fit only the first
  counter-diagonals, for speed purposes
- `get_ligation_ratio` helper

### Changed
- changed plot appearance for binless matrices, by default the median value is
  white

### Fixed
- fast binless can be run on matrices of up to 65535 bins. That means human
  chromosome 1 at 3.8kb resolution
- bug causing failed assertion during init on rare occasions in fast binless
- considerable speedup of init phase in fast binless
- installation link
- bug causing bin_data to fail if passed argument was not named cs
- segfaults caused by invalid Rcpp::wrap applied to unsigned (not supported)
- several warnings due to unsigned to int cast
- updated readme to indicate dependencies

## [0.12.0]
### Added
- First part of decay in fast binless is not forced to decrease. Use `free_diag`
 argument to decide (default 10kb)
- patchno column to fast binless output
- allow to pass a vector of lambda2 values to fast binless
- allow to pass lambda1 to get thresholdings similar to optimized binless
- allow to use final output of optimized binless as input to fast binless

### Changed
- Exposures are now fit separately
- Use `nobs` column in fast binless
- `ncounts` becomes `nobs` in optimized binless
- Use same column names in fast and optimized binless outputs
- fast binless model now mimicks that of slow binless
  - using binwise averages (breaks backwards compatibility)
  - using negative binomial with fixed dispersion
- require `pos1` and `pos2` columns in fast binless input
- fits in optimized binless are performed by group and are more efficient
  memorywise
- drop stan dependency and fit dispersion manually on a subset of rows, taking
  the median value
- Switch to BIC to estimate lambda2
- reorder arguments in fast_binless

### Fixed
- sigma parameter was ignored in `GeneralizedAdditiveModel`
- output matrix properly returns factor labels for name, bin1 and bin2 in fast binless
- bug causing weights to be twice too small
- centering bug on first iteration of fast binless
- bugs in the BIC calculation of optimized binless, causing lambda2 estimates to
  be too low and unstable
- Bug in fast binless when entire rows had nobs == 0

## [0.11.0]
### Changed
- spline base construction now migrated to C++ side
- fast binless now has smooth decay like optimized binless

### Fixed
- better storage and reporting of residuals
- correctly report when fused lasso does not converge
- provide means to increase nperf during binning
- updated documentation

### Removed
- requirement that rows and counter diagonals be nonzero in fast binless

## [0.10.2] - 2017-12-07
### Fixed
- bug causing decay to be flat in the first few basis functions

## [0.10.1] - 2017-11-29
### Fixed
- bug causing failure when a count was observed at the farthest distance

## [0.10.0] - 2017-11-28
### Added
- New function `zoom_csnorm` to take parts of a CSnorm object
- Remove large count outliers (possible PCR duplicates)
- Tutorial on arrow plots
### Changed
- Biases are not constrained while doing normalization
- Do not remove bad counter diagonals except very close to main diagonal
  (smoothing is good enough)
- Reorganize source code
### Fixed
- Better display and error messages for groupings
- Compress cts matrix during normalization for lower memory footprint
- Improved convergence: default bf_per_kb = 50 (avoids some oscillations)
- enlarge decay bins at border to avoid NAs
- pixels in plots are properly aligned with bin borders
- report missing data in binned matrices (fixes empty lines in binless matrices)
### Removed
- Deprecated code (exact model, stan spline helper functions and R graph calls)

## [0.9.0] - 2017-11-02
### Added
- Optimized binless
- Tutorials

## [0.2.0] - 2017-10-19
### Changed
- Propose only a fast and approximate binless implementation

## 0.1.0 - 2017-08-15
### Added
- Initial commit


[0.15.1]: ../../compare/v0.15.0...HEAD
[0.15.0]: ../../compare/v0.14.0...v0.15.0
[0.14.0]: ../../compare/v0.13.0...v0.14.0
[0.13.0]: ../../compare/v0.12.0...v0.13.0
[0.12.0]: ../../compare/v0.11.0...v0.12.0
[0.11.0]: ../../compare/v0.10.2...v0.11.0
[0.10.2]: ../../compare/v0.10.1...v0.10.2
[0.10.1]: ../../compare/v0.10.0...v0.10.1
[0.10.0]: ../../compare/v0.9.0...v0.10.0
[0.9.0]: ../../compare/v0.2.0...v0.9.0
[0.2.0]: ../../compare/v0.1.0...v0.2.0

