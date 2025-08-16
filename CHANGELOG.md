# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] -- 2025-08-16

### Added

- **Trait-based architecture**: New `Spline` trait enables specialized implementations for different types.
- **Error handling**: Functions now return `Result<T, SplineError>` instead of panicking.
  - `SplineError::InvalidKnotLength` for insufficient knots.
  - `SplineError::InvalidKnotPattern` for incorrect knot patterns.
- **Clone support**: Changed from `Copy` to `Clone` requirement, enabling support for types like `nalgebra::Vector3` etc.
- **Binary search optimization** in `spline_inverse`:
  - Automatically switches to binary search for splines with >12 segments
  - Provides O(log n) complexity instead of O(n).
  - Up to 7x speedup for 1,000+ segments.
- **Improved segment range detection**: Evaluates actual spline values at boundaries instead of just checking control points.
- **Prelude module**: Convenient `use uniform_cubic_splines::prelude::*` for common imports.
- Integration examples with nalgebra vector types.
- Comprehensive benchmarking examples.

### Changed

- **BREAKING**: All spline functions now return `Result` types for better error handling.
- **BREAKING**: Trait bounds changed from `Copy` to `Clone` for broader type support.
- **BREAKING**: Minimum supported Rust version updated to 1.83.
- `spline_inverse` algorithm significantly diverged from Open Shading Language implementation with major optimizations.
- Monotonicity check now performed once at the beginning of `spline_inverse`
- Adaptive root-finding uses Illinois modification of Regula Falsi with automatic bisection fallback.

### Performance

Current benchmarks on AMD Ryzen 7 6800H:

- `spline()`: ~2.4 ns for all segment counts.
- `spline_inverse()`:
  - 10 segments: 19 ns
  - 100 segments: 62 ns (was ~277 ns in v0.3)
  - 1000 segments: ~313 ns (7x faster than linear search)

## [0.3.3] - Previous version

### Changed

- Major performance optimizations from original OSL port.
- `spline()` function about twice as fast as original.
- `spline_inverse()` up to 3.5 times faster.

### Performance

- 10 points: 15.5 ns (spline), 95 ns (spline_inverse)
- 50 points: 16.4 ns (spline), 181.2 ns (spline_inverse)
- 100 points: 16.8 ns (spline), 277.7 ns (spline_inverse)
- 500 points: 20.6 ns (spline), 1.275 µs (spline_inverse)
