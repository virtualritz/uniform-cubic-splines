//! Trait-based spline evaluation system.
//!
//! This module provides the `Spline` trait which enables optimized
//! implementations for specific types while maintaining a generic fallback.

use crate::{Basis, SplineError, SplineInverseOptions, SplineResult};
use core::cmp::PartialOrd;
use core::num::NonZeroU16;
use core::ops::{Add, Div, Mul, Sub};
use lerp::Lerp;
use num_traits::{FromPrimitive, One, ToPrimitive, Zero};

/// Trait for types that can be interpolated using splines.
///
/// This trait enables specialized implementations for types that can benefit
/// from optimizations, while providing a generic implementation
/// for all suitable scalar types.
pub trait Spline: Sized {
    /// The scalar type used for interpolation parameters.
    type Scalar;

    /// Evaluates a spline at parameter `x` with the given `knots`.
    ///
    /// As `x` varies from 0 to 1, this returns the interpolated value.
    fn spline<B: Basis<Self::Scalar>>(
        x: Self::Scalar,
        knots: &[Self],
    ) -> SplineResult<Self>;

    /// Evaluates a single spline segment.
    ///
    /// This is the core operation that evaluates a cubic segment defined by
    /// 4 control points.
    fn spline_segment<B: Basis<Self::Scalar>>(
        x: Self::Scalar,
        cv: &[Self],
    ) -> Self;

    /// Computes the inverse of a spline.
    ///
    /// Returns the parameter `x` for which `spline(x)` would return `y`.
    fn spline_inverse<B: Basis<Self::Scalar>>(
        y: Self::Scalar,
        knots: &[Self],
        options: &SplineInverseOptions<Self::Scalar>,
    ) -> Option<Self::Scalar>;

    /// Computes the inverse of a single spline segment.
    ///
    /// Returns the parameter `x` in [0, 1] for which `spline_segment(x, cv)`
    /// would return `y`. Returns `None` if no solution exists within the
    /// segment.
    fn spline_inverse_segment<B: Basis<Self::Scalar>>(
        y: Self::Scalar,
        cv: &[Self],
        options: &SplineInverseOptions<Self::Scalar>,
    ) -> Option<Self::Scalar>;
}

// Helper to get default precision based on type size
fn default_precision<T: FromPrimitive>() -> T {
    // Use size_of to detect type, defaulting to 1e-6 for unknown types
    let size = core::mem::size_of::<T>();
    match size {
        2 => {
            T::from_f64(1.0e-3).unwrap_or_else(|| T::from_f64(1.0e-6).unwrap())
        } // f16
        4 => T::from_f64(1.0e-6).unwrap(), // f32
        8 => T::from_f64(1.0e-10).unwrap(), // f64
        16 => T::from_f64(1.0e-20)
            .unwrap_or_else(|| T::from_f64(1.0e-10).unwrap()), // f128
        _ => T::from_f64(1.0e-6).unwrap(), // default
    }
}

/// Generic implementation of Spline for scalar types that support the necessary
/// arithmetic operations. This works for f16, f32, f64, f128, and potentially
/// other numeric types that implement the required traits.
impl<T> Spline for T
where
    T: Clone
        + Zero
        + One
        + PartialOrd
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>
        + ToPrimitive
        + FromPrimitive
        + Lerp<T>,
{
    type Scalar = T;

    #[inline]
    fn spline<B: Basis<T>>(x: T, knots: &[T]) -> SplineResult<T> {
        // Validate knot length
        if knots.len() < 4 {
            return Err(SplineError::InvalidKnotLength {
                basis_name: B::NAME,
                min_knots: 4,
                actual: knots.len(),
            });
        } else if (knots.len() - 4) % B::STEP != 0 {
            return Err(SplineError::InvalidKnotPattern {
                basis_name: B::NAME,
                extra: B::STEP,
                actual: knots.len(),
            });
        }

        let number_of_segments = ((knots.len() - 4) / B::STEP) + 1;
        let num_seg_t = T::from_usize(number_of_segments).unwrap();

        // Scale x to find the segment and the interpolation parameter within
        // it. Manual clamp since we don't require Float trait.
        let x_clamped = if x < T::zero() {
            T::zero()
        } else if x > T::one() {
            T::one()
        } else {
            x
        };
        let x_scaled = x_clamped * num_seg_t;
        let mut segment = x_scaled.to_usize().unwrap_or(0);

        // Clamp segment to the last valid index.
        segment = segment.min(number_of_segments - 1);

        // The interpolation parameter, t, for the segment [0, 1].
        let t = x_scaled - T::from_usize(segment).unwrap();

        let start = segment * B::STEP;
        let cv = &knots[start..start + 4];

        Ok(Self::spline_segment::<B>(t, cv))
    }

    #[inline]
    #[allow(clippy::clone_on_copy)] // We use clone() for flexibility with non-Copy types
    fn spline_segment<B: Basis<T>>(x: T, cv: &[T]) -> T {
        let m = B::MATRIX;
        let c0 = &cv[0];
        let c1 = &cv[1];
        let c2 = &cv[2];
        let c3 = &cv[3];

        // Calculate the polynomial coefficients.
        // We use clone() here because T might not be Copy (e.g., custom vector
        // types)
        let w0 = c0.clone() * m[0][0].clone()
            + c1.clone() * m[0][1].clone()
            + c2.clone() * m[0][2].clone()
            + c3.clone() * m[0][3].clone();
        let w1 = c0.clone() * m[1][0].clone()
            + c1.clone() * m[1][1].clone()
            + c2.clone() * m[1][2].clone()
            + c3.clone() * m[1][3].clone();
        let w2 = c0.clone() * m[2][0].clone()
            + c1.clone() * m[2][1].clone()
            + c2.clone() * m[2][2].clone()
            + c3.clone() * m[2][3].clone();
        let w3 = c0.clone() * m[3][0].clone()
            + c1.clone() * m[3][1].clone()
            + c2.clone() * m[3][2].clone()
            + c3.clone() * m[3][3].clone();

        // Evaluate using Horner's method.
        ((w0 * x.clone() + w1) * x.clone() + w2) * x + w3
    }

    fn spline_inverse<B: Basis<T>>(
        y: T,
        knots: &[T],
        options: &SplineInverseOptions<T>,
    ) -> Option<T> {
        // Validate knot length
        if knots.len() < 4 || (knots.len() - 4) % B::STEP != 0 {
            // For inverse, we just return None on invalid input rather than
            // error
            return None;
        }

        // Extract options with defaults
        let max_iterations = options
            .max_iterations
            .unwrap_or(NonZeroU16::new(32).unwrap())
            .into();
        let precision = options
            .precision
            .clone()
            .unwrap_or_else(|| default_precision::<T>());

        // Perform monotonicity check once at the beginning
        #[cfg(feature = "monotonic_check")]
        {
            let increasing = knots[1] < knots[knots.len() - 2];
            for i in 1..knots.len() - 1 {
                let is_wrong_order = if increasing {
                    knots[i] > knots[i + 1]
                } else {
                    knots[i] < knots[i + 1]
                };
                if is_wrong_order {
                    // For inverse, return None instead of panicking
                    return None;
                }
            }
        }

        let number_of_segments = (knots.len() - 4) / B::STEP + 1;

        // Determine if knots are increasing or decreasing
        let increasing = knots[1] < knots[knots.len() - 2];

        // Handle boundary cases.
        let low_index = if B::STEP == 1 { 1 } else { 0 };
        let high_index = if B::STEP == 1 {
            knots.len() - 2
        } else {
            knots.len() - 1
        };

        if increasing {
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

        let inv_num_segments =
            T::one() / T::from_usize(number_of_segments).unwrap();

        // Optimization: Use binary search for > 12 segments, linear search
        // otherwise Benchmarking shows binary search becomes beneficial
        // around 10-20 segments Using 12 as threshold to be
        // conservative and avoid overhead for smaller counts
        if number_of_segments > 12 {
            // Binary search to find the segment containing y
            let mut left = 0;
            let mut right = number_of_segments - 1;

            while left <= right {
                let mid = (left + right) / 2;
                let start = mid * B::STEP;
                let cv = &knots[start..start + 4];

                // Get the actual range of this segment by evaluating at
                // boundaries
                let y0 = Self::spline_segment::<B>(T::zero(), cv);
                let y1 = Self::spline_segment::<B>(T::one(), cv);

                let (min_y, max_y) = if y0.clone() < y1.clone() {
                    (y0.clone(), y1.clone())
                } else {
                    (y1.clone(), y0.clone())
                };

                // Check if y is in this segment's range
                if y.clone() >= min_y && y.clone() <= max_y {
                    // Found potential segment, try to invert
                    let spline_on_segment =
                        |x_local: T| Self::spline_segment::<B>(x_local, cv);

                    if let Some(x_local) = crate::invert(
                        &spline_on_segment,
                        y.clone(),
                        T::zero(),
                        T::one(),
                        max_iterations,
                        precision.clone(),
                    ) {
                        let s_t = T::from_usize(mid).unwrap();
                        return Some((s_t + x_local) * inv_num_segments);
                    }
                }

                // Adjust search range based on segment start value
                if increasing {
                    if y.clone() < y0 {
                        right = mid.saturating_sub(1);
                    } else {
                        left = mid + 1;
                    }
                } else if y.clone() > y0 {
                    right = mid.saturating_sub(1);
                } else {
                    left = mid + 1;
                }

                if left > right {
                    break;
                }
            }

            // If binary search didn't find it, fall back to checking adjacent
            // segments This handles edge cases where segment
            // boundaries overlap
            let check_start = left.saturating_sub(1);
            let check_end = (right + 2).min(number_of_segments);

            for s in check_start..check_end {
                let start = s * B::STEP;
                let cv = &knots[start..start + 4];

                // Check actual segment range
                let y0 = Self::spline_segment::<B>(T::zero(), cv);
                let y1 = Self::spline_segment::<B>(T::one(), cv);
                let (min_y, max_y) = if y0.clone() < y1.clone() {
                    (y0, y1)
                } else {
                    (y1, y0)
                };

                if y.clone() >= min_y && y.clone() <= max_y {
                    let spline_on_segment =
                        |x_local: T| Self::spline_segment::<B>(x_local, cv);

                    if let Some(x_local) = crate::invert(
                        &spline_on_segment,
                        y.clone(),
                        T::zero(),
                        T::one(),
                        max_iterations,
                        precision.clone(),
                    ) {
                        let s_t = T::from_usize(s).unwrap();
                        return Some((s_t + x_local) * inv_num_segments);
                    }
                }
            }
        } else {
            // Linear search for small number of segments (<=12)
            // This is actually faster for small segment counts due to better
            // cache locality
            for s in 0..number_of_segments {
                let start = s * B::STEP;
                let cv = &knots[start..start + 4];

                // Determine the actual range of the segment
                let y0 = Self::spline_segment::<B>(T::zero(), cv);
                let y1 = Self::spline_segment::<B>(T::one(), cv);

                let (min_y, max_y) = if y0.clone() < y1.clone() {
                    (y0, y1)
                } else {
                    (y1, y0)
                };

                if y.clone() >= min_y && y.clone() <= max_y {
                    let spline_on_segment =
                        |x_local: T| Self::spline_segment::<B>(x_local, cv);

                    if let Some(x_local) = crate::invert(
                        &spline_on_segment,
                        y.clone(),
                        T::zero(),
                        T::one(),
                        max_iterations,
                        precision.clone(),
                    ) {
                        let s_t = T::from_usize(s).unwrap();
                        return Some((s_t + x_local) * inv_num_segments);
                    }
                }
            }
        }

        None
    }

    fn spline_inverse_segment<B: Basis<T>>(
        y: T,
        cv: &[T],
        options: &SplineInverseOptions<T>,
    ) -> Option<T> {
        // Extract options with defaults
        let max_iterations = options
            .max_iterations
            .unwrap_or(NonZeroU16::new(32).unwrap())
            .into();
        let precision = options
            .precision
            .clone()
            .unwrap_or_else(|| default_precision::<T>());

        // Create a closure for the segment spline function
        let spline_on_segment = |x: T| Self::spline_segment::<B>(x, cv);

        // Use the invert function to find the solution
        crate::invert(
            &spline_on_segment,
            y,
            T::zero(),
            T::one(),
            max_iterations,
            precision,
        )
    }
}
