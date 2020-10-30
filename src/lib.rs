#![no_std]
/// Cubic spline interpolation.
///
/// If you come from a background of Open or RenderMan shading language
/// this crate should feel like home.
///
/// The code is a Rust port of the resp. implementation found in the
/// Open Shading Language C++ source.
use core::ops::{Add, Mul};
use num_traits::{
    cast::{AsPrimitive, FromPrimitive},
    float::Float,
    identities::{One, Zero},
};

#[macro_use]
mod basis_macros;
pub mod basis;
use basis::*;

/// As `x` varies from `0` to `1`, this function returns the value
/// of a cubic interpolation of uniformly spaced `knots`.
/// The input value `x` will be clamped to the range `[0, 1]`.
///
/// Depending on the choosen [`Basis`] the length of the `knots`
/// parameter has certain constraints.
/// # Panics
/// If the `knots` slice has the wrong length this will panic when
/// the code is build with debug assertion enabled and produces
/// undefined behvior in a release build.
///
/// Use the [`is_len_ok()`] helper to check if a
/// knot slice you want to feed to this function has the correct
/// length.
/// # Examples
/// ```
/// use simple_spline::{spline, basis::CatmullRom};
///
/// //                 0.0  0.25 0.5  0.75 1.0
/// let knots = [-0.4, 0.0, 0.4, 0.5, 0.9, 1.0, 1.9];
///
/// assert!(spline::<CatmullRom, _, _>(0.25f64, &knots) == 0.4);
/// ```
pub fn spline<B, T, U>(x: T, knots: &[U]) -> U
where
    B: Basis<T>,
    T: AsPrimitive<usize> + Float + FromPrimitive + PartialOrd + One + Zero,
    U: Add<Output = U> + Copy + Mul<T, Output = U> + Zero,
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

    let mut x = x;
    x = *clamp(&mut x, num_traits::Zero::zero(), num_traits::One::one())
        * T::from_usize(number_of_segments).unwrap();

    let segment = x.floor();

    // x is the position along the segment.
    x = x - segment;

    let mut segment: usize = segment.as_();
    let segment_bound = number_of_segments - 1;
    if segment > segment_bound {
        segment = segment_bound;
    }

    let start = segment * B::STEP;

    // Get a slice for the segment.
    let cv = &knots[start..start + 4];

    let mut tk = (0..4)
        .map(|knot| {
            cv.iter()
                .zip(B::MATRIX[knot].iter())
                .fold(U::zero(), |total, (cv, basis)| total + *cv * *basis)
        });

    ((tk.next().unwrap() * x + tk.next().unwrap()) * x + tk.next().unwrap()) * x + tk.next().unwrap()
}

#[inline]
fn clamp<'a, T>(value: &'a mut T, min: T, max: T) -> &'a T
where
    T: PartialOrd,
{
    if *value < min {
        *value = min;
    } else if *value > max {
        *value = max;
    }
    value
}
