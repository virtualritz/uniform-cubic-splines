#![no_std]
#![cfg_attr(has_f16, feature(f16))] // FIXME: remove it when it's stablized
#![cfg_attr(has_f128, feature(f128))] // FIXME: remove it when it's stablized
#![feature(test)]
//! Uniform cubic spline interpolation & inversion.
//!
//! This crate supports the following types of splines:
//! * [B-spline](https://en.wikipedia.org/wiki/B-spline)
//! * [Bezier](https://en.wikipedia.org/wiki/Composite_B%C3%A9zier_curve)
//! * [Catmull-Rom](https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Catmull%E2%80%93Rom_spline)
//! * [Hermite](https://en.wikipedia.org/wiki/Cubic_Hermite_spline)
//! * Linear
//! * Power
//!
//! The crate uses generics to allow interpolation of any type for which
//! certain traits are defined.
//!
//! I.e. you can use this crate to interpolate splines in 1D, 2D, 3D, etc.
//!
//! ## Version 0.4 Changes
//!
//! * **Trait-based architecture**: New [`Spline`] trait enables specialized
//!   implementations.
//! * **Error handling**: Functions return [`Result<T,
//!   SplineError>`](SplineResult) instead of panicking.
//! * **Clone-based design**: Uses `Clone` instead of `Copy` to support more
//!   types like `nalgebra::Vector3`.
//! * **Extensibility**: External crates can implement [`Spline`] for their own
//!   types.
//!
//! ## Example
//!
//! Using a combination of [`spline_inverse()`] and [`spline()`] it is possible
//! to compute a full *spline-with-non-uniform-abscissæ*:
//!
//! ```
//! use uniform_cubic_splines::prelude::*;
//!
//! // We want to evaluate the spline at knot value 0.3.
//! let x = 0.3;
//!
//! // The first and last points are never interpolated.
//! let knot_spacing = [0.0, 0.0, 0.1, 0.3, 1.0, 1.0];
//! let knots = [0.0, 0.0, 1.3, 4.2, 3.2, 3.2];
//!
//! let v = spline_inverse::<CatmullRom, _>(x, &knot_spacing).unwrap();
//! let y = spline::<CatmullRom, _, _>(v, &knots).unwrap();
//!
//! assert!(y - 4.2 < 1e-6);
//! ```
//!
//! ## Background
//!
//! The code was originally a Rust port of the resp. implementations found in
//! the [Open Shading Language](https://github.com/imageworks/OpenShadingLanguage)
//! C++ source. However, it has since diverged significantly with extensive
//! optimizations:
//!
//! * **Binary search optimization**: `spline_inverse` uses binary search for
//!   splines with >12 segments, achieving O(log n) complexity instead of O(n)
//! * **Improved range detection**: Evaluates actual spline values at segment
//!   boundaries rather than just checking control points
//! * **Adaptive root-finding**: Illinois-modified Regula Falsi with automatic
//!   bisection fallback
//!
//! If you come from a background of computer graphics/shading languages used in
//! offline rendering this crate should feel familiar while providing better
//! performance.
//!
//! ## `no-std`
//!
//! The crate does not depend on the standard library (i.e. is marked `no_std`).
//!
//! ## Cargo Features
#![doc = document_features::document_features!()]

use core::cmp::PartialOrd;
use core::num::NonZeroU16;
use core::ops::{Add, Div, Mul, Sub};
use lerp::Lerp;
use num_traits::{FromPrimitive, One, Zero};

pub mod prelude {
    //! Convenience re-exports.
    pub use crate::basis::*;
    pub use crate::error::*;
    pub use crate::*;
}

#[macro_use]
mod basis_macros;
pub mod basis;
pub use basis::*;

mod error;
pub use error::{SplineError, SplineResult};

mod spline_trait;
pub use spline_trait::Spline;

/// Options for [`spline_inverse_with()`] function.
#[derive(Clone, Copy, Debug)]
pub struct SplineInverseOptions<T> {
    /// Controls the max. number of iterations of the
    /// [regular falsi](https://en.wikipedia.org/wiki/Regula_falsi) algorithm.
    ///
    /// If `None`, defaults to `32`.
    pub max_iterations: Option<NonZeroU16>,
    /// Controls the cutoff precision that is used to determine when the result
    /// is a good enough approximation, even if the specified number of
    /// `max_iterations` was *not* reached yet.
    ///
    /// If `None`, defaults depend on the type:
    /// - `f16`: `1.0e-3`
    /// - `f32`: `1.0e-6`
    /// - `f64`: `1.0e-10`
    /// - `f128`: `1.0e-20`
    /// - Other types: `1.0e-6`
    pub precision: Option<T>,
}

impl<T> Default for SplineInverseOptions<T> {
    fn default() -> Self {
        Self {
            max_iterations: None,
            precision: None,
        }
    }
}

/// As `x` varies from `0` to `1`, this function returns the value of a cubic
/// interpolation of uniformly spaced `knots`.
///
/// The input value `x` will be clamped to the range `[0, 1]`.
///
/// Depending on the chosen [`Basis`] the length of the `knots` parameter has
/// certain constraints. If these constraints are not honored, the function
/// returns an error.
///
/// # Errors
///
/// Returns [`SplineError::InvalidKnotLength`] if the knot vector has too few
/// knots. Returns [`SplineError::InvalidKnotPattern`] if the knot vector length
/// doesn't match the required pattern for the basis type.
///
/// Use the [`is_len_ok()`](crate::basis::Basis::is_len_ok()) helper to check if
/// a knot slice you want to feed to this function has the correct length.
///
/// # Examples
///
/// ```
/// use uniform_cubic_splines::{basis::CatmullRom, spline};
///
/// //                 0.0  0.25 0.5  0.75 1.0
/// let knots = [-0.4, 0.0, 0.4, 0.5, 0.9, 1.0, 1.9];
///
/// assert_eq!(0.4, spline::<CatmullRom, _, _>(0.25f64, &knots).unwrap());
/// ```
pub fn spline<B, T, U>(x: T, knots: &[U]) -> SplineResult<U>
where
    B: Basis<T>,
    U: Spline<Scalar = T>,
{
    U::spline::<B>(x, knots)
}

/// Evaluates a spline for a single segment defined by 4 control points.
/// This is the performance-critical inner loop of the `spline` function.
#[inline]
pub fn spline_segment<B, T, U>(x: T, cv: &[U]) -> U
where
    B: Basis<T>,
    U: Spline<Scalar = T>,
{
    U::spline_segment::<B>(x, cv)
}

/// Computes the inverse of a single spline segment.
///
/// This returns the value `x` in [0, 1] for which `spline_segment(x, cv)` would
/// return `y`. Returns `None` if no solution exists within the segment.
///
/// The `cv` parameter must contain exactly 4 control values.
///
/// # Examples
///
/// ```
/// use uniform_cubic_splines::{
///     basis::Linear, spline_inverse_segment_with, SplineInverseOptions,
/// };
///
/// let cv = [0.0, 0.25, 0.75, 1.0];
/// let y = 0.5;
///
/// let x = spline_inverse_segment_with::<Linear, _>(
///     y,
///     &cv,
///     &SplineInverseOptions::default(),
/// );
/// assert!(x.is_some());
/// ```
pub fn spline_inverse_segment<B, T>(y: T, cv: &[T]) -> Option<T>
where
    B: Basis<T>,
    T: Spline<Scalar = T>,
{
    T::spline_inverse_segment::<B>(y, cv, &SplineInverseOptions::default())
}

/// Computes the inverse of a single spline segment with custom options.
///
/// This returns the value `x` in [0, 1] for which `spline_segment(x, cv)` would
/// return `y`. Returns `None` if no solution exists within the segment.
///
/// The `cv` parameter must contain exactly 4 control values.
pub fn spline_inverse_segment_with<B, T>(
    y: T,
    cv: &[T],
    options: &SplineInverseOptions<T>,
) -> Option<T>
where
    B: Basis<T>,
    T: Spline<Scalar = T>,
{
    T::spline_inverse_segment::<B>(y, cv, options)
}

/// Computes the inverse of the [`spline()`] function.
///
/// This returns the value `x` for which `spline(x)` would return `y`.
///
/// Results are undefined if the `knots` do not specifiy a monotonic (only
/// increasing or only decreasing) set of values.
///
/// Depending on the choosen [`Basis`] the length of the `knots` parameter has
/// certain constraints. If these constraints are not honored the code will
/// produce undefined results in a `release` build. See below.
///
/// If no solution can be found the function returns `None`.
///
/// The underlying algorithm uses the
/// [regular falsi](https://en.wikipedia.org/wiki/Regula_falsi) method to find
/// the solution. If you need more control over the algorithm use the
/// [`spline_inverse_with()`] function.
///
/// # Panics
///
/// If the `monotonic_check` feature is enabled this will panic if the `knots`
/// slice is not monotonic. Note: disabling `monotonic_check` can improve
/// performance by ~5-10% but will produce undefined results for non-monotonic
/// knots. If the `knots` slice has the wrong length this will panic with a
/// resp. error message when the code is built with debug assertions *enabled*.
///
/// Use the [`is_len_ok()`](crate::basis::Basis::is_len_ok()) helper to check if
/// a knot slice you want to feed to this function has the correct length.
///
/// # Examples
///
/// ```
/// # use uniform_cubic_splines::prelude::*;
/// let knots = [0.0, 0.0, 0.5, 0.5];
///
/// assert_eq!(Some(0.5), spline_inverse::<Linear, _>(0.25f64, &knots));
/// ```
pub fn spline_inverse<B, T>(y: T, knots: &[T]) -> Option<T>
where
    B: Basis<T>,
    T: Spline<Scalar = T>,
{
    T::spline_inverse::<B>(y, knots, &SplineInverseOptions::default())
}

/// Computes the inverse of the [`spline()`] function with control over
/// iterations & precision via resp. [`SplineInverseOptions`].
///
/// See the [`spline_inverse()`] function for details; including panics.
///
/// # Examples
///
/// ```
/// # use core::num::NonZeroU16;
/// # use uniform_cubic_splines::prelude::*;
/// let knots = [0.0, 0.0, 0.5, 0.5];
///
/// // Use custom max iterations
/// assert_eq!(
///     Some(0.5),
///     spline_inverse_with::<Linear, _>(
///         0.25f64,
///         &knots,
///         &SplineInverseOptions {
///             max_iterations: Some(NonZeroU16::new(16).unwrap()),
///             precision: None, // Use default precision for f64 (1e-10)
///         }
///     )
/// );
///
/// // Or use custom precision with f32
/// let knots_f32 = [0.0f32, 0.0, 0.5, 0.5];
/// assert_eq!(
///     Some(0.5),
///     spline_inverse_with::<Linear, _>(
///         0.25f32,
///         &knots_f32,
///         &SplineInverseOptions {
///             max_iterations: None,  // Use default (32)
///             precision: Some(1e-8), // Custom precision
///         }
///     )
/// );
/// ```
pub fn spline_inverse_with<B, T>(
    y: T,
    knots: &[T],
    options: &SplineInverseOptions<T>,
) -> Option<T>
where
    B: Basis<T>,
    T: Spline<Scalar = T>,
{
    // Delegate to the trait implementation which handles validation
    T::spline_inverse::<B>(y, knots, options)
}

/// Returns `true` if a `knots` slice you want to feed into [`spline()`] has the
/// correct length for the choosen [`Basis`].
#[deprecated(
    since = "0.3.0",
    note = "Use the resp. basis's `is_len_ok()` function."
)]
pub fn is_len_ok<B, T>(len: usize) -> bool
where
    B: Basis<T>,
{
    B::is_len_ok(len)
}

#[inline]
pub(crate) fn invert<T>(
    function: &impl Fn(T) -> T,
    y: T,
    x_min: T,
    x_max: T,
    max_iterations: u16,
    epsilon: T,
) -> Option<T>
where
    T: Clone
        + PartialOrd
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + FromPrimitive
        + Lerp<T>,
{
    // Optimized root-finding with adaptive strategy.
    let mut v0 = function(x_min.clone());
    let mut v1 = function(x_max.clone());

    let increasing = v0.clone() < v1.clone();
    let (vmin, vmax) = if increasing {
        (v0.clone(), v1.clone())
    } else {
        (v1.clone(), v0.clone())
    };

    // Early exit if y is outside the range.
    if !(vmin <= y.clone() && y.clone() <= vmax) {
        return None;
    }

    // Check if we're already at the boundaries.
    let two = T::from_f64(2.0).unwrap();
    let diff0 = y.clone() - v0.clone();
    let abs_diff0 = if diff0 < T::zero() {
        T::zero() - diff0
    } else {
        diff0
    };
    if abs_diff0 < epsilon.clone() * two.clone() {
        return Some(x_min);
    }
    let diff1 = y.clone() - v1.clone();
    let abs_diff1 = if diff1 < T::zero() {
        T::zero() - diff1
    } else {
        diff1
    };
    if abs_diff1 < epsilon.clone() * two {
        return Some(x_max);
    }

    // Adaptive strategy: Illinois modification of Regula Falsi to prevent
    // stalling, with automatic fallback to bisection if convergence is
    // slow.
    let mut x_min = x_min;
    let mut x_max = x_max;
    let mut last_side = 0i8; // Track which side was updated last.
    let mut stuck_count = 0u8; // Count iterations stuck on same side.

    // Dynamic switching threshold based on problem difficulty.
    let rf_iterations = (3 * max_iterations) / 4;

    for iters in 0..max_iterations {
        // Adaptive epsilon based on current interval.
        let interval_diff = x_max.clone() - x_min.clone();
        let interval_size = if interval_diff < T::zero() {
            T::zero() - interval_diff
        } else {
            interval_diff
        };
        let threshold = interval_size.clone() * T::from_f64(1e-12).unwrap();
        let adaptive_epsilon = if epsilon.clone() > threshold {
            epsilon.clone()
        } else {
            threshold
        };

        // Check for convergence with adaptive tolerance.
        if interval_size < adaptive_epsilon {
            return Some(x_min.lerp(x_max, T::from_f64(0.5).unwrap()));
        }

        // Compute interpolation factor.
        let t = if iters < rf_iterations && stuck_count < 3 {
            // Modified Regula Falsi with Illinois modification.
            let mut t_rf = (y.clone() - v0.clone()) / (v1.clone() - v0.clone());

            // Apply Illinois modification if we're stuck on one side.
            if stuck_count > 0 {
                if last_side > 0 {
                    // Stuck on right side, reduce left weight.
                    v0 = v0 * T::from_f64(0.5).unwrap();
                } else {
                    // Stuck on left side, reduce right weight.
                    v1 = v1 * T::from_f64(0.5).unwrap();
                }
                t_rf = (y.clone() - v0.clone()) / (v1.clone() - v0.clone());
            }

            // Ensure t is in valid range, otherwise bisect.
            if t_rf > T::from_f64(0.01).unwrap()
                && t_rf < T::from_f64(0.99).unwrap()
            {
                t_rf
            } else {
                T::from_f64(0.5).unwrap()
            }
        } else {
            // Bisection for guaranteed convergence.
            T::from_f64(0.5).unwrap()
        };

        let x = x_min.clone().lerp(x_max.clone(), t);
        let v = function(x.clone());

        // Update interval and track convergence.
        let side = if (v.clone() < y.clone()) == increasing {
            x_min = x.clone();
            v0 = v.clone();
            1i8
        } else {
            x_max = x.clone();
            v1 = v.clone();
            -1i8
        };

        // Track if we're stuck on same side.
        if side == last_side {
            stuck_count += 1;
        } else {
            stuck_count = 0;
        }
        last_side = side;

        // Check value convergence.
        let v_diff = v.clone() - y.clone();
        let abs_v_diff = if v_diff < T::zero() {
            T::zero() - v_diff
        } else {
            v_diff
        };
        if abs_v_diff < epsilon {
            return Some(x);
        }
    }

    // Return best guess after max iterations.
    Some(x_min.lerp(x_max, T::from_f64(0.5).unwrap()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn len_validation_bspline() {
        // B-spline has STEP=1, so any length >= 4 works.
        assert!(<Bspline as Basis<f64>>::is_len_ok(4));
        assert!(<Bspline as Basis<f64>>::is_len_ok(7));
        assert!(!<Bspline as Basis<f64>>::is_len_ok(3));
    }

    #[test]
    fn len_validation_catmull_rom() {
        // Catmull-Rom has STEP=1, so any length >= 4 works.
        assert!(<CatmullRom as Basis<f64>>::is_len_ok(4));
        assert!(<CatmullRom as Basis<f64>>::is_len_ok(8));
        assert!(!<CatmullRom as Basis<f64>>::is_len_ok(2));
    }

    #[test]
    fn len_validation_linear() {
        // Linear has STEP=1, so any length >= 4 works.
        assert!(<Linear as Basis<f64>>::is_len_ok(4));
        assert!(<Linear as Basis<f64>>::is_len_ok(10));
        assert!(!<Linear as Basis<f64>>::is_len_ok(1));
    }

    #[test]
    fn len_validation_bezier() {
        // Bezier has STEP=3, so (len-4) must be divisible by 3.
        assert!(<Bezier as Basis<f64>>::is_len_ok(4)); // (4-4)/3 = 0, valid
        assert!(<Bezier as Basis<f64>>::is_len_ok(7)); // (7-4)/3 = 1, valid
        assert!(<Bezier as Basis<f64>>::is_len_ok(10)); // (10-4)/3 = 2, valid
        assert!(!<Bezier as Basis<f64>>::is_len_ok(5)); // (5-4)/3 has remainder
        assert!(!<Bezier as Basis<f64>>::is_len_ok(8)); // (8-4)/3 has remainder
    }

    #[test]
    fn len_validation_hermite() {
        // Hermite has STEP=2, so (len-4) must be divisible by 2.
        assert!(<Hermite as Basis<f64>>::is_len_ok(4)); // (4-4)/2 = 0, valid
        assert!(<Hermite as Basis<f64>>::is_len_ok(6)); // (6-4)/2 = 1, valid
        assert!(<Hermite as Basis<f64>>::is_len_ok(8)); // (8-4)/2 = 2, valid
        assert!(!<Hermite as Basis<f64>>::is_len_ok(5)); // (5-4)/2 has remainder
        assert!(!<Hermite as Basis<f64>>::is_len_ok(7)); // (7-4)/2 has
                                                         // remainder
    }

    #[test]
    fn len_validation_power() {
        // Power has STEP=4, so (len-4) must be divisible by 4.
        assert!(<Power as Basis<f64>>::is_len_ok(4)); // (4-4)/4 = 0, valid
        assert!(<Power as Basis<f64>>::is_len_ok(8)); // (8-4)/4 = 1, valid
        assert!(<Power as Basis<f64>>::is_len_ok(12)); // (12-4)/4 = 2, valid
        assert!(!<Power as Basis<f64>>::is_len_ok(5)); // (5-4)/4 has remainder
        assert!(!<Power as Basis<f64>>::is_len_ok(6)); // (6-4)/4 has remainder
        assert!(!<Power as Basis<f64>>::is_len_ok(7)); // (7-4)/4 has remainder
    }

    #[test]
    fn spline_evaluation_smoke_test() {
        // Basic smoke test to ensure spline evaluation works.
        let knots = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0];
        let result = spline::<CatmullRom, _, _>(0.5, &knots);
        // Just check it doesn't panic.
        assert!(result.is_ok());
    }

    #[test]
    fn bezier_invalid_length_returns_error() {
        let knots = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0]; // 8 knots, invalid for Bezier.
        let result = spline::<Bezier, _, _>(0.5, &knots);
        assert!(result.is_err());
    }

    #[test]
    fn hermite_invalid_length_returns_error() {
        let knots = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0]; // 7 knots, invalid for Hermite.
        let result = spline::<Hermite, _, _>(0.5, &knots);
        assert!(result.is_err());
    }

    #[test]
    fn power_invalid_length_returns_error() {
        let knots = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0, 0.5]; // 9 knots, invalid for Power (not 4n+4).
        let result = spline::<Power, _, _>(0.5, &knots);
        assert!(result.is_err());
    }

    #[test]
    fn bspline_too_few_knots_returns_error() {
        let knots = [0.0, 1.0, 2.0]; // Only 3 knots.
        let result = spline::<Bspline, _, _>(0.5, &knots);
        assert!(result.is_err());
    }

    #[test]
    fn bezier_valid_lengths_work() {
        // Test that valid Bezier lengths work. STEP=3, so (len-4) must be
        // divisible by 3.
        let knots4 = [0.0, 0.0, 1.0, 4.0]; // (4-4)/3 = 0, valid.
        let result = spline::<Bezier, _, _>(0.5, &knots4);
        assert!(result.is_ok());

        let knots7 = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0]; // (7-4)/3 = 1, valid.
        let result = spline::<Bezier, _, _>(0.5, &knots7);
        assert!(result.is_ok());

        let knots10 = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0, 2.0, 3.0]; // (10-4)/3 = 2, valid.
        let result = spline::<Bezier, _, _>(0.5, &knots10);
        assert!(result.is_ok());
    }

    #[test]
    fn hermite_valid_lengths_work() {
        // Test that valid Hermite lengths work. STEP=2, so (len-4) must be
        // divisible by 2.
        let knots4 = [0.0, 0.0, 1.0, 4.0]; // (4-4)/2 = 0, valid.
        let result = spline::<Hermite, _, _>(0.5, &knots4);
        assert!(result.is_ok());

        let knots6 = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0]; // (6-4)/2 = 1, valid.
        let result = spline::<Hermite, _, _>(0.5, &knots6);
        assert!(result.is_ok());

        let knots8 = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0]; // (8-4)/2 = 2, valid.
        let result = spline::<Hermite, _, _>(0.5, &knots8);
        assert!(result.is_ok());
    }

    #[test]
    fn power_valid_lengths_work() {
        // Test that valid Power lengths work. STEP=4, so (len-4) must be
        // divisible by 4.
        let knots4 = [0.0, 0.0, 1.0, 4.0]; // (4-4)/4 = 0, valid.
        let result = spline::<Power, _, _>(0.5, &knots4);
        assert!(result.is_ok());

        let knots8 = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0]; // (8-4)/4 = 1, valid.
        let result = spline::<Power, _, _>(0.5, &knots8);
        assert!(result.is_ok());

        let knots12 =
            [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0]; // (12-4)/4 = 2, valid.
        let result = spline::<Power, _, _>(0.5, &knots12);
        assert!(result.is_ok());
    }
}
