//! Spline basis types.
use crate::*;

/// A cubic spline basis.
///
/// Some basis types require a particular number of knot values.
///
/// `Bezier` splines require *4×n+3* values.
///
/// `Hermite` splines require *4×n+2* values.
///
/// `Power` splines require *4×n+4* values.
///
/// `B-spline`, `CatmullRom` and `Linear` splines may use any
/// number of values with *n≥4*.
pub trait Basis<T: Float> {
    const NAME: &'static str;
    const STEP: usize;
    const MATRIX: [[T; 4]; 4];
    const EXTRA_KNOTS: usize;
    /// Returns `true` if a `knots` slice you want to feed into
    /// [`spline()`] has the correct length for this `Basis`.
    fn is_len_ok(len: usize) -> bool {
        if 0 == Self::EXTRA_KNOTS {
            4 <= len
        } else {
            4 + Self::EXTRA_KNOTS <= len && 0 == (len - Self::EXTRA_KNOTS) % 4
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
#[cfg(has_f16)]
b_spline_basis!(f16);
#[cfg(has_f128)]
b_spline_basis!(f128);

/// A *Bezier* spline basis.
///
/// It requires *4×n+3* values in the knot vector.
pub struct Bezier;

bezier_basis!(f64);
bezier_basis!(f32);
#[cfg(has_f16)]
bezier_basis!(f16);
#[cfg(has_f128)]
bezier_basis!(f128);

/// A *Catmull-Rom* spline basis.
///
/// It requires *≥4* values in the knot vector.
pub struct CatmullRom;

catmull_rom_basis!(f64);
catmull_rom_basis!(f32);
#[cfg(has_f16)]
catmull_rom_basis!(f16);
#[cfg(has_f128)]
catmull_rom_basis!(f128);

/// A *Hermite* spline basis.
///
/// It requires *4×n+2* values in the knot vector.
pub struct Hermite;

hermite_basis!(f64);
hermite_basis!(f32);
#[cfg(has_f16)]
hermite_basis!(f16);
#[cfg(has_f128)]
hermite_basis!(f128);

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
#[cfg(has_f16)]
linear_basis!(f16);
#[cfg(has_f128)]
linear_basis!(f128);

/// A *Power* basis.
///
/// It requires *4×n+4* values in the knot vector.
pub struct Power;

power_basis!(f64);
power_basis!(f32);
#[cfg(has_f16)]
power_basis!(f16);
#[cfg(has_f128)]
power_basis!(f128);
