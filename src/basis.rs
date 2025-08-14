//! Spline basis types.
use crate::*;

/// A cubic spline basis.
///
/// Some basis types require a particular number of knot values,
/// where *n* is a positive integer (n ≥ 1):
///
/// `Bezier` splines require *4×n+3* values.
///
/// `Hermite` splines require *4×n+2* values.
///
/// `Power` splines require *4×n+4* values.
///
/// `B-spline`, `CatmullRom` and `Linear` splines may use any
/// number of values with a minimum of 4 knots.
pub trait Basis<T: Float> {
    const NAME: &'static str;
    const STEP: usize;
    const MATRIX: [[T; 4]; 4];
    const EXTRA_KNOTS: usize;
    /// Returns `true` if a `knots` slice you want to feed into
    /// [`spline()`] has the correct length for this `Basis`.
    fn is_len_ok(len: usize) -> bool {
        if Self::EXTRA_KNOTS == 0 {
            len >= 4
        } else {
            len >= 4 + Self::EXTRA_KNOTS && (len - Self::EXTRA_KNOTS) % 4 == 0
        }
    }
}

/// A *B-spline* basis.
///
/// It requires *≥4* values in the knot vector.
pub struct Bspline;

#[deprecated(since = "0.3.0", note = "Use `Bspline` instead")]
pub type BSpline = Bspline;

b_spline_basis!(f64);
b_spline_basis!(f32);

/// A *Bezier* spline basis.
///
/// It requires *4×n+3* values in the knot vector.
pub struct Bezier;

bezier_basis!(f64);
bezier_basis!(f32);

/// A *Catmull-Rom* spline basis.
///
/// It requires *≥4* values in the knot vector.
pub struct CatmullRom;

catmull_rom_basis!(f64);
catmull_rom_basis!(f32);

/// A *Hermite* spline basis.
///
/// It requires *4×n+2* values in the knot vector.
pub struct Hermite;

hermite_basis!(f64);
hermite_basis!(f32);

/// A *Linear* basis.
///
/// It requires *≥4* values in the knot vector.
///
/// To maintain consistency with the other spline types, `Linear`
/// splines will ignore the first and last data value, interpolating
/// piecewise linearly between *y₁* and *yₙ₋₂* (with *n* being the
/// knot array index, starting from `0`).
pub struct Linear;

linear_basis!(f64);
linear_basis!(f32);

/// A *Power* basis.
///
/// It requires *4×n+4* values in the knot vector.
pub struct Power;

power_basis!(f64);
power_basis!(f32);
