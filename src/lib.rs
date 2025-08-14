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
//! let y = spline::<CatmullRom, _, _>(v, &knots);
//!
//! assert!(y - 4.2 < 1e-6);
//! ```
//!
//! ## Background
//!
//! The code is a Rust port of the resp. implementations found in the
//! [Open Shading Language](https://github.com/imageworks/OpenShadingLanguage)
//! C++ source but was optimized quite a bit afterwards.
//!
//! If you come from a background of computer graphics/shading languages used in
//! offline rendering this crate should feel like home.
//!
//! ## `no-std`
//!
//! The crate does not depend on the standard library (i.e. is marked `no_std`).
//!
//! ## Cargo Features
#![doc = document_features::document_features!()]

use core::{
    num::NonZeroU16,
    ops::{Add, Mul},
};
use lerp::Lerp;
use num_traits::{cast::FromPrimitive, float::Float, identities::Zero};

pub mod prelude {
    //! Convenience re-exports.
    pub use crate::basis::*;
    pub use crate::*;
}

#[macro_use]
mod basis_macros;
pub mod basis;
pub use basis::*;

/// Options for [`spline_inverse_with()`] function.
#[derive(Clone, Copy, Debug)]
pub struct SplineInverseOptions<T> {
    /// Controls the max. number of iterations of the
    /// [regular falsi](https://en.wikipedia.org/wiki/Regula_falsi) algorithm.
    ///
    /// The default is `32`.
    pub max_iterations: NonZeroU16,
    /// Controls the cutoff precision that is used to determine when the result
    /// is a good enough approximation, even if the specified number of
    /// `max_iterations` was *not* reached yet.
    ///
    /// The default is `1.0e-6`.
    pub precision: T,
}

impl<T> Default for SplineInverseOptions<T>
where
    T: FromPrimitive,
{
    fn default() -> Self {
        Self {
            max_iterations: NonZeroU16::new(32).unwrap(),
            precision: T::from_f64(1.0e-6).unwrap(),
        }
    }
}

#[cfg(debug_assertions)]
macro_rules! len_check {
    ($knots:ident) => {
        // UX
        if $knots.len() < 4 + B::EXTRA_KNOTS {
            panic!(
                "{} curve must have at least {} knots. Found: {}.",
                B::NAME,
                4 + B::EXTRA_KNOTS,
                $knots.len()
            );
        } else if (B::EXTRA_KNOTS != 0)
            && (($knots.len() - B::EXTRA_KNOTS) % 4 != 0)
        {
            panic!(
                "{} curve must have 4×𝘯+{} knots. Found: {}.",
                B::NAME,
                B::EXTRA_KNOTS,
                $knots.len()
            );
        }
    };
}

/// As `x` varies from `0` to `1`, this function returns the value of a cubic
/// interpolation of uniformly spaced `knots`.
///
/// The input value `x` will be clamped to the range `[0, 1]`.
///
/// Depending on the choosen [`Basis`] the length of the `knots` parameter has
/// certain constraints. If these constraints are not honored the code will
/// produce undefined results in a `release` build. See below.
///
/// # Panics
///
/// If the `knots` slice has the wrong length this will panic with a resp. error
/// message when the code is built with debug assertions *enabled*.
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
/// assert_eq!(0.4, spline::<CatmullRom, _, _>(0.25f64, &knots));
/// ```
pub fn spline<B, T, U>(x: T, knots: &[U]) -> U
where
    B: Basis<T>,
    T: Float + FromPrimitive,
    U: Add<Output = U> + Clone + Mul<T, Output = U> + Zero,
{
    #[cfg(debug_assertions)]
    len_check!(knots);

    let number_of_segments = ((knots.len() - 4) / B::STEP) + 1;
    let num_seg_t = T::from_usize(number_of_segments).unwrap();

    // Scale x to find the segment and the interpolation parameter within it.
    let x_scaled = x.clamp(T::zero(), T::one()) * num_seg_t;
    let mut segment = x_scaled.to_usize().unwrap_or(0);

    // Clamp segment to the last valid index.
    segment = segment.min(number_of_segments - 1);

    // The interpolation parameter, t, for the segment [0, 1].
    let t = x_scaled - T::from_usize(segment).unwrap();

    let start = segment * B::STEP;
    let cv = &knots[start..start + 4];

    spline_segment::<B, T, U>(t, cv)
}

/// Evaluates a spline for a single segment defined by 4 control points.
/// This is the performance-critical inner loop of the `spline` function.
#[inline]
fn spline_segment<B, T, U>(x: T, cv: &[U]) -> U
where
    B: Basis<T>,
    T: Float,
    U: Add<Output = U> + Clone + Mul<T, Output = U> + Zero,
{
    let m = B::MATRIX;
    let c0 = &cv[0];
    let c1 = &cv[1];
    let c2 = &cv[2];
    let c3 = &cv[3];

    // Calculate the polynomial coefficients by transforming the control points
    // with the basis matrix. This is equivalent to `w = M * C`.
    // The matrix rows are grouped by powers of x for Horner's method.
    let w0 = c0.clone() * m[0][0]
        + c1.clone() * m[0][1]
        + c2.clone() * m[0][2]
        + c3.clone() * m[0][3];
    let w1 = c0.clone() * m[1][0]
        + c1.clone() * m[1][1]
        + c2.clone() * m[1][2]
        + c3.clone() * m[1][3];
    let w2 = c0.clone() * m[2][0]
        + c1.clone() * m[2][1]
        + c2.clone() * m[2][2]
        + c3.clone() * m[2][3];
    let w3 = c0.clone() * m[3][0]
        + c1.clone() * m[3][1]
        + c2.clone() * m[3][2]
        + c3.clone() * m[3][3];

    // Evaluate the polynomial `((w0*x + w1)*x + w2)*x + w3` using Horner's
    // method.
    // This is highly efficient and should allow the compiler to generate fused
    // multiply-add (FMA) instructions.
    ((w0 * x + w1) * x + w2) * x + w3
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
/// slice is not monotonic.
/// If the `knots` slice has the wrong length this will panic with a resp. error
/// message when the code is built with debug assertions *enabled*.
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
    T: Float + FromPrimitive,
{
    spline_inverse_with::<B, T>(y, knots, &SplineInverseOptions::default())
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
/// assert_eq!(
///     Some(0.5),
///     spline_inverse_with::<Linear, _>(
///         0.25f64,
///         &knots,
///         &SplineInverseOptions {
///             max_iterations: NonZeroU16::new(16).unwrap(),
///             ..Default::default()
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
    T: Float + FromPrimitive,
{
    #[cfg(debug_assertions)]
    len_check!(knots);

    #[cfg(feature = "monotonic_check")]
    if !knots.is_sorted() {
        panic!("The knots array fed to spline_inverse() is not monotonic.");
    }

    // Account for out-of-range inputs by clamping to the spline's boundary
    // values.
    let low_index = if B::STEP == 1 { 1 } else { 0 };
    let high_index = if B::STEP == 1 {
        knots.len() - 2
    } else {
        knots.len() - 1
    };

    // If increasing ...
    if knots[1] < knots[knots.len() - 2] {
        if y <= knots[low_index] {
            return Some(T::zero());
        }
        if y >= knots[high_index] {
            return Some(T::one());
        }
    } else {
        if y >= knots[low_index] {
            return Some(T::zero());
        }
        if y <= knots[high_index] {
            return Some(T::one());
        }
    }

    let number_of_segments = (knots.len() - 4) / B::STEP + 1;
    let inv_num_segments =
        T::one() / T::from_usize(number_of_segments).unwrap();

    // Search each segment for the value.
    for s in 0..number_of_segments {
        let start = s * B::STEP;
        let cv = &knots[start..start + 4];

        // This closure is cheap as it only operates on the 4 control points
        // of the current segment. `invert` will search for `x_local` in [0, 1].
        let spline_on_segment =
            |x_local: T| spline_segment::<B, T, T>(x_local, cv);

        if let Some(x_local) = invert(
            &spline_on_segment,
            y,
            T::zero(),
            T::one(),
            options.max_iterations.into(),
            options.precision,
        ) {
            // Convert the local solution `x_local` in [0, 1] back to the global
            // coordinate space, also in [0, 1].
            let s_t = T::from_usize(s).unwrap();
            return Some((s_t + x_local) * inv_num_segments);
        }
    }

    None
}

/// Returns `true` if a `knots` slice you want to feed into [`spline()`] has the
/// correct length for the choosen [`Basis`].
#[deprecated(
    since = "0.3.0",
    note = "Use the resp. basis's `is_len_ok()` function."
)]
pub fn is_len_ok<B>(len: usize) -> bool
where
    B: Basis<f32>,
{
    if 0 == B::EXTRA_KNOTS {
        4 <= len
    } else {
        4 + B::EXTRA_KNOTS <= len && 0 == (len - B::EXTRA_KNOTS) % 4
    }
}

#[inline]
fn invert<T>(
    function: &impl Fn(T) -> T,
    y: T,
    x_min: T,
    x_max: T,
    max_iterations: u16,
    epsilon: T,
) -> Option<T>
where
    T: Float + FromPrimitive + Lerp<T>,
{
    // Use the Regula Falsi method, falling back to bisection if it
    // struggles to converge. This is a robust approach for root-finding.
    let mut v0 = function(x_min);
    let mut v1 = function(x_max);

    let mut x = x_min;
    let increasing = v0 < v1;

    let (vmin, vmax) = if increasing { (v0, v1) } else { (v1, v0) };

    // If y is outside the range of this segment, there's no solution here.
    if !(vmin <= y && y <= vmax) {
        return None;
    }

    // Already close enough at the boundaries.
    if (v0 - v1).abs() < epsilon {
        return Some(x);
    }

    // Switch to bisection if Regula Falsi hasn't converged after 3/4 of the
    // maximum number of iterations.
    // See, e.g., "Numerical Recipes" for the basic ideas behind both methods.
    let rf_iterations = (3 * max_iterations) / 4;

    let mut x_min = x_min;
    let mut x_max = x_max;

    for iters in 0..max_iterations {
        // Interpolation factor.
        let t = if iters < rf_iterations {
            // Regula falsi.
            let t_rf = (y - v0) / (v1 - v0);
            if t_rf > T::zero() && t_rf < T::one() {
                t_rf
            } else {
                // RF convergence failure (e.g., due to curvature), bisect
                // instead.
                T::from_f64(0.5).unwrap()
            }
        } else {
            // Bisection.
            T::from_f64(0.5).unwrap()
        };
        x = x_min.lerp(x_max, t);

        let v = function(x);
        if (v < y) == increasing {
            x_min = x;
            v0 = v;
        } else {
            x_max = x;
            v1 = v;
        }
        // Check for convergence.
        if (x_max - x_min).abs() < epsilon || (v - y).abs() < epsilon {
            return Some(x);
        }
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn len_validation_bspline() {
        // B-spline requires >= 4 knots.
        assert!(<Bspline as Basis<f64>>::is_len_ok(4));
        assert!(<Bspline as Basis<f64>>::is_len_ok(5));
        assert!(<Bspline as Basis<f64>>::is_len_ok(10));
        assert!(<Bspline as Basis<f64>>::is_len_ok(100));
        assert!(!<Bspline as Basis<f64>>::is_len_ok(3));
        assert!(!<Bspline as Basis<f64>>::is_len_ok(2));
        assert!(!<Bspline as Basis<f64>>::is_len_ok(1));
        assert!(!<Bspline as Basis<f64>>::is_len_ok(0));
    }

    #[test]
    fn len_validation_catmull_rom() {
        // Catmull-Rom requires >= 4 knots.
        assert!(<CatmullRom as Basis<f64>>::is_len_ok(4));
        assert!(<CatmullRom as Basis<f64>>::is_len_ok(5));
        assert!(<CatmullRom as Basis<f64>>::is_len_ok(10));
        assert!(<CatmullRom as Basis<f64>>::is_len_ok(100));
        assert!(!<CatmullRom as Basis<f64>>::is_len_ok(3));
        assert!(!<CatmullRom as Basis<f64>>::is_len_ok(2));
        assert!(!<CatmullRom as Basis<f64>>::is_len_ok(1));
        assert!(!<CatmullRom as Basis<f64>>::is_len_ok(0));
    }

    #[test]
    fn len_validation_linear() {
        // Linear requires >= 4 knots.
        assert!(<Linear as Basis<f64>>::is_len_ok(4));
        assert!(<Linear as Basis<f64>>::is_len_ok(5));
        assert!(<Linear as Basis<f64>>::is_len_ok(10));
        assert!(<Linear as Basis<f64>>::is_len_ok(100));
        assert!(!<Linear as Basis<f64>>::is_len_ok(3));
        assert!(!<Linear as Basis<f64>>::is_len_ok(2));
        assert!(!<Linear as Basis<f64>>::is_len_ok(1));
        assert!(!<Linear as Basis<f64>>::is_len_ok(0));
    }

    #[test]
    fn len_validation_bezier() {
        // Bezier requires 4×n+3 knots (n >= 1).
        assert!(<Bezier as Basis<f64>>::is_len_ok(7)); // 4×1+3
        assert!(<Bezier as Basis<f64>>::is_len_ok(11)); // 4×2+3
        assert!(<Bezier as Basis<f64>>::is_len_ok(15)); // 4×3+3
        assert!(<Bezier as Basis<f64>>::is_len_ok(19)); // 4×4+3
        assert!(<Bezier as Basis<f64>>::is_len_ok(23)); // 4×5+3

        // Invalid lengths.
        assert!(!<Bezier as Basis<f64>>::is_len_ok(3));
        assert!(!<Bezier as Basis<f64>>::is_len_ok(4));
        assert!(!<Bezier as Basis<f64>>::is_len_ok(5));
        assert!(!<Bezier as Basis<f64>>::is_len_ok(6));
        assert!(!<Bezier as Basis<f64>>::is_len_ok(8));
        assert!(!<Bezier as Basis<f64>>::is_len_ok(9));
        assert!(!<Bezier as Basis<f64>>::is_len_ok(10));
        assert!(!<Bezier as Basis<f64>>::is_len_ok(12));
        assert!(!<Bezier as Basis<f64>>::is_len_ok(13));
        assert!(!<Bezier as Basis<f64>>::is_len_ok(14));
    }

    #[test]
    fn len_validation_hermite() {
        // Hermite requires 4×n+2 knots (n >= 1).
        assert!(<Hermite as Basis<f64>>::is_len_ok(6)); // 4×1+2
        assert!(<Hermite as Basis<f64>>::is_len_ok(10)); // 4×2+2
        assert!(<Hermite as Basis<f64>>::is_len_ok(14)); // 4×3+2
        assert!(<Hermite as Basis<f64>>::is_len_ok(18)); // 4×4+2
        assert!(<Hermite as Basis<f64>>::is_len_ok(22)); // 4×5+2

        // Invalid lengths.
        assert!(!<Hermite as Basis<f64>>::is_len_ok(2));
        assert!(!<Hermite as Basis<f64>>::is_len_ok(3));
        assert!(!<Hermite as Basis<f64>>::is_len_ok(4));
        assert!(!<Hermite as Basis<f64>>::is_len_ok(5));
        assert!(!<Hermite as Basis<f64>>::is_len_ok(7));
        assert!(!<Hermite as Basis<f64>>::is_len_ok(8));
        assert!(!<Hermite as Basis<f64>>::is_len_ok(9));
        assert!(!<Hermite as Basis<f64>>::is_len_ok(11));
        assert!(!<Hermite as Basis<f64>>::is_len_ok(12));
        assert!(!<Hermite as Basis<f64>>::is_len_ok(13));
    }

    #[test]
    fn len_validation_power() {
        // Power requires 4×n+4 knots (n >= 1).
        assert!(<Power as Basis<f64>>::is_len_ok(8)); // 4×1+4
        assert!(<Power as Basis<f64>>::is_len_ok(12)); // 4×2+4
        assert!(<Power as Basis<f64>>::is_len_ok(16)); // 4×3+4
        assert!(<Power as Basis<f64>>::is_len_ok(20)); // 4×4+4
        assert!(<Power as Basis<f64>>::is_len_ok(24)); // 4×5+4

        // Invalid lengths.
        assert!(!<Power as Basis<f64>>::is_len_ok(4));
        assert!(!<Power as Basis<f64>>::is_len_ok(5));
        assert!(!<Power as Basis<f64>>::is_len_ok(6));
        assert!(!<Power as Basis<f64>>::is_len_ok(7));
        assert!(!<Power as Basis<f64>>::is_len_ok(9));
        assert!(!<Power as Basis<f64>>::is_len_ok(10));
        assert!(!<Power as Basis<f64>>::is_len_ok(11));
        assert!(!<Power as Basis<f64>>::is_len_ok(13));
        assert!(!<Power as Basis<f64>>::is_len_ok(14));
        assert!(!<Power as Basis<f64>>::is_len_ok(15));
    }

    #[test]
    fn spline_evaluation_smoke_test() {
        // Basic smoke test to ensure spline evaluation works.
        let knots = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0];
        let result = spline::<CatmullRom, _, _>(0.5, &knots);
        // Just check it doesn't panic and returns a reasonable value.
        assert!(result.is_finite());
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "Bezier curve must have 4×𝘯+3 knots. Found: 8")]
    fn bezier_invalid_length_panics() {
        let knots = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0]; // 8 knots, invalid for Bezier.
        let _ = spline::<Bezier, _, _>(0.5, &knots);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "Hermite curve must have 4×𝘯+2 knots. Found: 7")]
    fn hermite_invalid_length_panics() {
        let knots = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0]; // 7 knots, invalid for Hermite.
        let _ = spline::<Hermite, _, _>(0.5, &knots);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "Power curve must have 4×𝘯+4 knots. Found: 9")]
    fn power_invalid_length_panics() {
        let knots = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0, 0.5]; // 9 knots, invalid for Power (not 4n+4).
        let _ = spline::<Power, _, _>(0.5, &knots);
    }

    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(
        expected = "B-spline curve must have at least 4 knots. Found: 3"
    )]
    fn bspline_too_few_knots_panics() {
        let knots = [0.0, 1.0, 2.0]; // Only 3 knots.
        let _ = spline::<Bspline, _, _>(0.5, &knots);
    }

    #[test]
    fn bezier_valid_lengths_work() {
        // Test that valid Bezier lengths don't panic.
        let knots7 = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0]; // 7 = 4×1+3.
        let result = spline::<Bezier, _, _>(0.5, &knots7);
        assert!(result.is_finite());

        let knots11 = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0, 2.0, 3.0, 4.0]; // 11 = 4×2+3.
        let result = spline::<Bezier, _, _>(0.5, &knots11);
        assert!(result.is_finite());
    }

    #[test]
    fn hermite_valid_lengths_work() {
        // Test that valid Hermite lengths don't panic.
        let knots6 = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0]; // 6 = 4×1+2.
        let result = spline::<Hermite, _, _>(0.5, &knots6);
        assert!(result.is_finite());

        let knots10 = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0, 2.0, 3.0]; // 10 = 4×2+2.
        let result = spline::<Hermite, _, _>(0.5, &knots10);
        assert!(result.is_finite());
    }

    #[test]
    fn power_valid_lengths_work() {
        // Test that valid Power lengths don't panic.
        let knots8 = [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0]; // 8 = 4×1+4.
        let result = spline::<Power, _, _>(0.5, &knots8);
        assert!(result.is_finite());

        let knots12 =
            [0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0, 1.0, 2.0, 3.0, 4.0, 5.0]; // 12 = 4×2+4.
        let result = spline::<Power, _, _>(0.5, &knots12);
        assert!(result.is_finite());
    }
}
