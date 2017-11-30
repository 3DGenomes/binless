# Changelog
All notable changes to *binless* will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixed
- better storage and reporting of residuals

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


[Unreleased]: ../../compare/v0.10.1...HEAD
[0.10.1]: ../../compare/v0.10.0...v0.10.1
[0.10.0]: ../../compare/v0.9.0...v0.10.0
[0.9.0]: ../../compare/v0.2.0...v0.9.0
[0.2.0]: ../../compare/v0.1.0...v0.2.0

