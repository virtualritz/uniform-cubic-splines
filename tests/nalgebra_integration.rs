//! Integration tests demonstrating spline interpolation with nalgebra types.
//!
//! This shows how to create wrapper types to use nalgebra vectors with the
//! Spline trait.

use nalgebra::{Point3, Vector3};
use uniform_cubic_splines::{spline, CatmullRom, Spline, SplineResult};

/// Wrapper type for nalgebra Vector3 to implement Spline trait
#[derive(Clone, Copy, Debug)]
struct SplineVec3(Vector3<f32>);

impl From<Vector3<f32>> for SplineVec3 {
    fn from(v: Vector3<f32>) -> Self {
        SplineVec3(v)
    }
}

impl From<SplineVec3> for Vector3<f32> {
    fn from(v: SplineVec3) -> Self {
        v.0
    }
}

impl SplineVec3 {
    fn new(x: f32, y: f32, z: f32) -> Self {
        SplineVec3(Vector3::new(x, y, z))
    }
}

// Implement Spline trait for the wrapper type
impl Spline for SplineVec3 {
    type Scalar = f32;

    fn spline<B: uniform_cubic_splines::Basis<f32>>(
        x: f32,
        knots: &[Self],
    ) -> SplineResult<Self> {
        // Use component-wise spline evaluation
        let x_coords: Vec<f32> = knots.iter().map(|v| v.0.x).collect();
        let y_coords: Vec<f32> = knots.iter().map(|v| v.0.y).collect();
        let z_coords: Vec<f32> = knots.iter().map(|v| v.0.z).collect();

        let x_result = spline::<B, _, _>(x, &x_coords)?;
        let y_result = spline::<B, _, _>(x, &y_coords)?;
        let z_result = spline::<B, _, _>(x, &z_coords)?;

        Ok(SplineVec3(Vector3::new(x_result, y_result, z_result)))
    }

    fn spline_segment<B: uniform_cubic_splines::Basis<f32>>(
        x: f32,
        cv: &[Self],
    ) -> Self {
        use uniform_cubic_splines::spline_segment;

        let x_cv = [cv[0].0.x, cv[1].0.x, cv[2].0.x, cv[3].0.x];
        let y_cv = [cv[0].0.y, cv[1].0.y, cv[2].0.y, cv[3].0.y];
        let z_cv = [cv[0].0.z, cv[1].0.z, cv[2].0.z, cv[3].0.z];

        let x_result = spline_segment::<B, _, _>(x, &x_cv);
        let y_result = spline_segment::<B, _, _>(x, &y_cv);
        let z_result = spline_segment::<B, _, _>(x, &z_cv);

        SplineVec3(Vector3::new(x_result, y_result, z_result))
    }

    fn spline_inverse<B: uniform_cubic_splines::Basis<f32>>(
        _y: f32,
        _knots: &[Self],
        _options: &uniform_cubic_splines::SplineInverseOptions<f32>,
    ) -> Option<f32> {
        // Inverse doesn't make sense for vector-valued splines
        unimplemented!(
            "spline_inverse is not defined for vector-valued splines"
        )
    }

    fn spline_inverse_segment<B: uniform_cubic_splines::Basis<f32>>(
        _y: f32,
        _cv: &[Self],
        _options: &uniform_cubic_splines::SplineInverseOptions<f32>,
    ) -> Option<f32> {
        // Inverse doesn't make sense for vector-valued splines
        unimplemented!(
            "spline_inverse_segment is not defined for vector-valued splines"
        )
    }
}

/// Wrapper type for nalgebra Point3
#[derive(Clone, Copy, Debug)]
struct SplinePoint3(Point3<f32>);

impl From<Point3<f32>> for SplinePoint3 {
    fn from(p: Point3<f32>) -> Self {
        SplinePoint3(p)
    }
}

impl From<SplinePoint3> for Point3<f32> {
    fn from(p: SplinePoint3) -> Self {
        p.0
    }
}

impl SplinePoint3 {
    fn new(x: f32, y: f32, z: f32) -> Self {
        SplinePoint3(Point3::new(x, y, z))
    }
}

impl Spline for SplinePoint3 {
    type Scalar = f32;

    fn spline<B: uniform_cubic_splines::Basis<f32>>(
        x: f32,
        knots: &[Self],
    ) -> SplineResult<Self> {
        let x_coords: Vec<f32> = knots.iter().map(|p| p.0.x).collect();
        let y_coords: Vec<f32> = knots.iter().map(|p| p.0.y).collect();
        let z_coords: Vec<f32> = knots.iter().map(|p| p.0.z).collect();

        let x_result = spline::<B, _, _>(x, &x_coords)?;
        let y_result = spline::<B, _, _>(x, &y_coords)?;
        let z_result = spline::<B, _, _>(x, &z_coords)?;

        Ok(SplinePoint3(Point3::new(x_result, y_result, z_result)))
    }

    fn spline_segment<B: uniform_cubic_splines::Basis<f32>>(
        x: f32,
        cv: &[Self],
    ) -> Self {
        use uniform_cubic_splines::spline_segment;

        let x_cv = [cv[0].0.x, cv[1].0.x, cv[2].0.x, cv[3].0.x];
        let y_cv = [cv[0].0.y, cv[1].0.y, cv[2].0.y, cv[3].0.y];
        let z_cv = [cv[0].0.z, cv[1].0.z, cv[2].0.z, cv[3].0.z];

        let x_result = spline_segment::<B, _, _>(x, &x_cv);
        let y_result = spline_segment::<B, _, _>(x, &y_cv);
        let z_result = spline_segment::<B, _, _>(x, &z_cv);

        SplinePoint3(Point3::new(x_result, y_result, z_result))
    }

    fn spline_inverse<B: uniform_cubic_splines::Basis<f32>>(
        _y: f32,
        _knots: &[Self],
        _options: &uniform_cubic_splines::SplineInverseOptions<f32>,
    ) -> Option<f32> {
        unimplemented!("spline_inverse is not defined for point-valued splines")
    }

    fn spline_inverse_segment<B: uniform_cubic_splines::Basis<f32>>(
        _y: f32,
        _cv: &[Self],
        _options: &uniform_cubic_splines::SplineInverseOptions<f32>,
    ) -> Option<f32> {
        unimplemented!(
            "spline_inverse_segment is not defined for point-valued splines"
        )
    }
}

#[test]
fn test_vector3_spline() {
    let knots = vec![
        SplineVec3::new(0.0, 0.0, 0.0),
        SplineVec3::new(1.0, 0.5, 0.25),
        SplineVec3::new(2.0, 1.0, 1.0),
        SplineVec3::new(3.0, 0.5, 2.25),
        SplineVec3::new(4.0, 0.0, 4.0),
        SplineVec3::new(5.0, -0.5, 6.25),
        SplineVec3::new(6.0, -1.0, 9.0),
    ];

    // Evaluate spline at t=0.5
    let result = spline::<CatmullRom, _, _>(0.5, &knots).unwrap();
    let vec: Vector3<f32> = result.into();

    // The result should be a valid Vector3
    assert!(vec.x.is_finite());
    assert!(vec.y.is_finite());
    assert!(vec.z.is_finite());

    // Test at boundaries
    let result_start = spline::<CatmullRom, _, _>(0.0, &knots).unwrap();
    let result_end = spline::<CatmullRom, _, _>(1.0, &knots).unwrap();

    let vec_start: Vector3<f32> = result_start.into();
    let vec_end: Vector3<f32> = result_end.into();

    assert!(vec_start.x.is_finite());
    assert!(vec_end.x.is_finite());
}

#[test]
fn test_point3_spline() {
    let knots = vec![
        SplinePoint3::new(0.0, 0.0, 0.0),
        SplinePoint3::new(1.0, 2.0, 3.0),
        SplinePoint3::new(2.0, 1.0, 2.0),
        SplinePoint3::new(3.0, 3.0, 1.0),
        SplinePoint3::new(4.0, 2.0, 0.0),
        SplinePoint3::new(5.0, 1.0, 1.0),
    ];

    // Evaluate spline at various points
    for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let result = spline::<CatmullRom, _, _>(t, &knots).unwrap();
        let point: Point3<f32> = result.into();
        assert!(point.x.is_finite());
        assert!(point.y.is_finite());
        assert!(point.z.is_finite());
    }
}

#[test]
fn test_segment_evaluation() {
    let cv = [
        SplineVec3::new(0.0, 0.0, 0.0),
        SplineVec3::new(1.0, 1.0, 1.0),
        SplineVec3::new(2.0, 0.5, 2.0),
        SplineVec3::new(3.0, 0.0, 3.0),
    ];

    // Test segment evaluation at various t values
    for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let result = SplineVec3::spline_segment::<CatmullRom>(t, &cv);
        let vec: Vector3<f32> = result.into();
        assert!(vec.x.is_finite());
        assert!(vec.y.is_finite());
        assert!(vec.z.is_finite());
    }
}

#[test]
fn test_different_basis_types() {
    use uniform_cubic_splines::{Bezier, Bspline, Hermite, Linear};

    let knots = vec![
        SplineVec3::new(0.0, 0.0, 0.0),
        SplineVec3::new(1.0, 1.0, 1.0),
        SplineVec3::new(2.0, 0.0, 2.0),
        SplineVec3::new(3.0, -1.0, 3.0),
        SplineVec3::new(4.0, 0.0, 4.0),
        SplineVec3::new(5.0, 1.0, 5.0),
        SplineVec3::new(6.0, 0.0, 6.0),
    ];

    // Test with different basis types
    let _ = spline::<Bspline, _, _>(0.5, &knots).unwrap();
    let _ = spline::<CatmullRom, _, _>(0.5, &knots).unwrap();
    let _ = spline::<Linear, _, _>(0.5, &knots).unwrap();

    // Bezier has STEP=3, so it works with 4, 7, 10, 13... knots
    let bezier_knots = vec![
        SplineVec3::new(0.0, 0.0, 0.0),
        SplineVec3::new(1.0, 1.0, 1.0),
        SplineVec3::new(2.0, 0.0, 2.0),
        SplineVec3::new(3.0, -1.0, 3.0),
        SplineVec3::new(4.0, 0.0, 4.0),
        SplineVec3::new(5.0, 1.0, 5.0),
        SplineVec3::new(6.0, 0.0, 6.0),
    ];
    let _ = spline::<Bezier, _, _>(0.5, &bezier_knots).unwrap();

    // Hermite has STEP=2, so it works with 4, 6, 8, 10... knots
    let hermite_knots = vec![
        SplineVec3::new(0.0, 0.0, 0.0),
        SplineVec3::new(1.0, 1.0, 1.0),
        SplineVec3::new(2.0, 0.0, 2.0),
        SplineVec3::new(3.0, -1.0, 3.0),
        SplineVec3::new(4.0, 0.0, 4.0),
        SplineVec3::new(5.0, 1.0, 5.0),
    ];
    let _ = spline::<Hermite, _, _>(0.5, &hermite_knots).unwrap();
}

#[test]
fn test_convenience_conversions() {
    // Test that conversions between wrapper and nalgebra types work smoothly
    let nalgebra_vec = Vector3::new(1.0, 2.0, 3.0);
    let spline_vec: SplineVec3 = nalgebra_vec.into();
    let back_to_nalgebra: Vector3<f32> = spline_vec.into();

    assert_eq!(nalgebra_vec, back_to_nalgebra);

    let nalgebra_point = Point3::new(4.0, 5.0, 6.0);
    let spline_point: SplinePoint3 = nalgebra_point.into();
    let back_to_point: Point3<f32> = spline_point.into();

    assert_eq!(nalgebra_point, back_to_point);
}
