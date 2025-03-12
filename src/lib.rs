#![no_std]
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
//! ## Cargo Features
//!
//! * `monotonic_check` -- The [`spline_inverse()`] code will check if the knot
//!   vector is monotonic (on by default).
//!
//! The crate does not depend on the standard library (i.e. is marked `no_std`).
//!
//! ## Example
//!
//! Using a combination of [`spline()`] and [`spline_inverse()`] it is possible
//! to compute a full spline-with-non-uniform-abscissæ:
//!
//! ```
//! use uniform_cubic_splines::{basis::CatmullRom, spline, spline_inverse};
//!
//! // We want to evaluate the spline at knot value 0.3.
//! let x = 0.3;
//!
//! // The first and last points are never interpolated.
//! let knot_spacing = [0.0, 0.0, 0.1, 0.3, 1.0, 1.0];
//! let knots = [0.0, 0.0, 1.3, 4.2, 3.2, 3.2];
//!
//! let v =
//!     spline_inverse::<CatmullRom, _>(x, &knot_spacing, None, None).unwrap();
//! let y = spline::<CatmullRom, _, _>(v, &knots);
//!
//! assert!(y - 4.2 < 1e-6);
//! ```
//!
//! ## Background
//!
//! The code is a Rust port of the resp. implementations found in the
//! [Open Shading Language](https://github.com/imageworks/OpenShadingLanguage)
//! C++ source.
//!
//! If you come from a background of computer graphics/shading languages used in
//! offline rendering this crate should feel like home.
use core::ops::{Add, Mul};
use lerp::Lerp;
use num_traits::{
    cast::{AsPrimitive, FromPrimitive},
    float::Float,
    identities::{One, Zero},
};

#[macro_use]
mod basis_macros;
pub mod basis;
use basis::*;

/// As `x` varies from `0` to `1`, this function returns the value of a cubic
/// interpolation of uniformly spaced `knots`.
///
/// The input value `x` will be clamped to the range `[0, 1]`.
///
/// Depending on the choosen [`Basis`] the length of the `knots` parameter has
/// certain constraints. If these constraints are not honored the code will
/// produce undefined results in a `release` build.
///
/// # Panics
///
/// If the `knots` slice has the wrong length this will panic with a resp. error
/// message when the code is built with debug assertion enabled.
///
/// Use the [`is_len_ok()`] helper to check if a knot slice you want to feed to
/// this function has the correct length.
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
    T: AsPrimitive<usize> + Float + FromPrimitive + PartialOrd + One + Zero,
    U: Add<Output = U> + Clone + Mul<T, Output = U> + Zero,
{
    // UX
    #[cfg(debug_assertions)]
    if knots.len() < 4 + B::EXTRA_KNOTS {
        panic!(
            "{} curve must have at least {} knots. Found: {}.",
            B::NAME,
            4 + B::EXTRA_KNOTS,
            knots.len()
        );
    } else if (B::EXTRA_KNOTS != 0) && ((knots.len() - B::EXTRA_KNOTS) % 4 == 0)
    {
        panic!(
            "{} curve must have 4×𝘯+{} knots. Found: {}.",
            B::NAME,
            B::EXTRA_KNOTS,
            knots.len()
        );
    }

    let number_of_segments: usize = ((knots.len() - 4) / B::STEP) + 1;

    let mut x = clamp(x, Zero::zero(), One::one())
        * T::from_usize(number_of_segments).unwrap();

    let mut segment: usize = x.as_();

    let segment_bound = number_of_segments - 1;
    if segment > segment_bound {
        segment = segment_bound;
    }

    // x is the position along the segment.
    x = x - T::from_usize(segment).unwrap();

    let start = segment * B::STEP;

    // Get a slice for the segment.
    let cv = &knots[start..start + 4];

    B::MATRIX
        .iter()
        .map(|row| {
            cv.iter()
                .zip(row.iter())
                .fold(U::zero(), |total, (cv, basis)| {
                    total + cv.clone() * *basis
                })
        })
        .fold(Zero::zero(), |acc, elem| acc * x + elem)
}

/// Computes the inverse of the [`spline()`] function.
///
/// This returns the value `x` for which `spline(x)` would return `y`.
///
/// Results are undefined if the `knots` do not specifiy a monotonic (only
/// increasing or only decreasing) set of values.
///
/// If no solution can be found the function returns `None`.
///
/// The underlying algorithm uses the
/// [regular falsi](https://en.wikipedia.org/wiki/Regula_falsi) method to find
/// the solution.
///
/// The `iterations` parameter controls the max. number of iterations of this
/// this algorithm. If omitted, the default is `32`.
///
/// The `precision` parameter controls the cutoff precision that is used to
/// determine when the result is a good enough approximation, even if the
/// specified number of `iterations` was not reached yet. If omitted, the
/// default is `1e-6`.
///
/// # Panics
///
/// If the `monotonic_check` feature is enabled this will panic if the `knots`
/// slice is not monotonic.
///
/// # Examples
///
/// ```
/// use uniform_cubic_splines::{basis::Linear, spline_inverse};
///
/// let knots = [0.0, 0.0, 0.5, 0.5];
///
/// assert_eq!(
///     Some(0.5),
///     spline_inverse::<Linear, _>(0.25f64, &knots, None, None)
/// );
/// ```
pub fn spline_inverse<B, T>(
    y: T,
    knots: &[T],
    iterations: Option<usize>,
    precision: Option<T>,
) -> Option<T>
where
    B: Basis<T>,
    T: AsPrimitive<usize> + Float + FromPrimitive + PartialOrd + One + Zero,
{
    #[cfg(feature = "monotonic_check")]
    if !knots.is_sorted() {
        panic!("The knots array fed to spline_inverse() is not monotonic.");
    }

    // Account for out-of-range inputs;
    // just clamp to the values we have.
    let low_index: usize = if B::STEP == 1 { 1 } else { 0 };

    let high_index = if B::STEP == 1 {
        knots.len() - 2
    } else {
        knots.len() - 1
    };

    // If increasing ...
    if knots[1] < knots[knots.len() - 2] {
        if y <= knots[low_index] {
            return Some(Zero::zero());
        }
        if y >= knots[high_index] {
            return Some(One::one());
        }
    } else {
        if y >= knots[low_index] {
            return Some(Zero::zero());
        }
        if y <= knots[high_index] {
            return Some(One::one());
        }
    }

    let spline_function = |x| spline::<B, T, T>(x, knots);

    let number_of_segments = (knots.len() - 4) / B::STEP + 1;
    let number_of_segments_inverted = 1.0 / number_of_segments as f64;

    // Search each interval.
    let mut r0 = num_traits::Zero::zero();

    for s in 0..number_of_segments {
        let r1 =
            T::from_f64(number_of_segments_inverted * (s + 1) as f64).unwrap();

        if let Some(x) = invert(
            &spline_function,
            y,
            r0,
            r1,
            iterations.unwrap_or(32),
            precision.unwrap_or(T::from_f64(1.0e-6).unwrap()),
        ) {
            return Some(x);
        }

        // Start of next interval is end of this one.
        r0 = r1;
    }

    None
}

/// Returns `true` if a `knots` slice you want to feed into
/// [`spline()`] has the correct length for the choosen [`Basis`].
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
    function: &dyn Fn(T) -> T,
    y: T,
    x_min: T,
    x_max: T,
    max_iterations: usize,
    epsilon: T,
) -> Option<T>
where
    T: AsPrimitive<usize> + Float + FromPrimitive + PartialOrd + One + Zero,
{
    // Use the Regula Falsi method, falling back to bisection if it
    // hasn't converged after 3/4 of the maximum number of iterations.
    // See, e.g., "Numerical Recipes" for the basic ideas behind both
    // methods.
    let mut v0 = function(x_min);
    let mut v1 = function(x_max);

    let mut x = x_min;
    let increasing = v0 < v1;

    let vmin = if increasing { v0 } else { v1 };
    let vmax = if increasing { v1 } else { v0 };

    if !(vmin <= y && y <= vmax) {
        return None;
    }

    // Already close enough.
    if Float::abs(v0 - v1) < epsilon {
        return Some(x);
    }

    // How many times to try regula falsi.
    let rf_iterations = (3 * max_iterations) / 4;

    let mut x_min = x_min;
    let mut x_max = x_max;

    for iters in 0..max_iterations {
        // Interpolation factor.
        let mut t: T;
        if iters < rf_iterations {
            // Regula falsi.
            t = (y - v0) / (v1 - v0);
            if t <= num_traits::Zero::zero() || t >= num_traits::One::one() {
                // RF convergence failure -- bisect instead.
                t = T::from_f64(0.5).unwrap();
            }
        } else {
            // Bisection.
            t = T::from_f64(0.5).unwrap();
        }
        x = x_min.lerp(x_max, t);

        let v = function(x);
        if (v < y) == increasing {
            x_min = x;
            v0 = v;
        } else {
            x_max = x;
            v1 = v;
        }
        if Float::abs(x_max - x_min) < epsilon || Float::abs(v - y) < epsilon {
            return Some(x); // converged
        }
    }
    Some(x)
}

#[inline]
fn clamp<T>(value: T, min: T, max: T) -> T
where
    T: PartialOrd,
{
    if value < min {
        min
    } else if value > max {
        max
    } else {
        value
    }
}
