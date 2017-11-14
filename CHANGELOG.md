# Changelog
All notable changes to *binless* will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- New function `zoom_csnorm` to take parts of a CSnorm object
- Remove large count outliers (possible PCR duplicates)
### Fixed
- Better display and error messages for groupings
- Compress cts matrix during normalization for lower memory footprint
- Improved convergence: default bf_per_kb = 50 (avoids some oscillations)
- enlarge decay bins at border to avoid NAs

## [0.9.0] - 2017-11-02
### Added
- Optimized binless
- Tutorials

## 0.2.0 - 2017-10-19
### Changed
- Propose only a fast and approximate binless implementation

## 0.1.0 - 2017-08-15
### Added
- Initial commit


[Unreleased]: https://github.com/3DGenomes/csnorm/compare/v0.9.0...HEAD
[0.9.0]: https://github.com/3DGenomes/csnorm/compare/v0.2.0...v0.9.0

