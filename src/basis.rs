use crate::*;

/// A cubic spline basis.
/// Some basis types require a particular number of knot values.
///
/// `Bezier` splines require `4n+3` values.
///
/// `Hermite` splines require `4n+2` values.
///
/// `B-sppline`, `CatmullRom` and `Linear` splines may use any
/// number of values with *n≥4*.
pub trait Basis<T: Float> {
    const NAME: &'static str;
    const STEP: usize;
    const MATRIX: [[T; 4]; 4];
    const EXTRA_KNOTS: usize;
}

/// A *B-spline* basis. It require ≥4 values in the knot vector.
pub struct BSpline;

b_spline_basis!(f64);
b_spline_basis!(f32);

/// A *Bezie* spline basis. It require *4×𝘯+3* values in the knot vector.
pub struct Bezier;

bezier_basis!(f64);
bezier_basis!(f32);

/// A *Catmull-Rom* spline basis. It require ≥4 values in the knot vector.
pub struct CatmullRom;

catmull_rom_basis!(f64);
catmull_rom_basis!(f32);

/// A *Hermite* spline basis. It require *4×𝘯+2* values in the knot vector.
pub struct Hermite;

hermite_basis!(f64);
hermite_basis!(f32);

/// A *Linear* basis. It require ≥4 values in the knot vector.
///
/// To maintain consistency with the other spline types, `Linear`
/// splines will ignore the first and last data value, interpolating
/// piecewise linearly between *y₁* and *yₙ₋₂*.
pub struct Linear;

linear_basis!(f64);
linear_basis!(f32);

/// A *Power* basis. It require *4×𝘯+4* values in the knot vector.
pub struct Power;

power_basis!(f64);
power_basis!(f32);
