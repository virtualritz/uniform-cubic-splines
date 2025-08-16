//! Example exploring minimal trait bounds for spline interpolation
//!
//! This demonstrates what traits are actually needed vs what we currently
//! require.

use std::ops::{Add, Mul, Sub};

// What we actually need for spline_segment (the core operation)
trait MinimalSpline: Sized + Clone {
    type Scalar: Clone;

    // For scalar * vector (used in basis matrix multiplication)
    fn mul_scalar(&self, scalar: Self::Scalar) -> Self;

    // For vector + vector
    fn add(&self, other: &Self) -> Self;
}

// Example implementation for a simple 2D point
#[derive(Clone, Debug)]
struct Point2D {
    x: f64,
    y: f64,
}

impl MinimalSpline for Point2D {
    type Scalar = f64;

    fn mul_scalar(&self, scalar: f64) -> Self {
        Point2D {
            x: self.x * scalar,
            y: self.y * scalar,
        }
    }

    fn add(&self, other: &Self) -> Self {
        Point2D {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

// This is essentially what spline_segment does
fn spline_segment_minimal<T: MinimalSpline<Scalar = f64>>(
    x: f64,
    cv: &[T; 4],
    matrix: [[f64; 4]; 4],
) -> T {
    // Calculate polynomial coefficients using only mul_scalar and add
    let w0 = cv[0]
        .mul_scalar(matrix[0][0])
        .add(&cv[1].mul_scalar(matrix[0][1]))
        .add(&cv[2].mul_scalar(matrix[0][2]))
        .add(&cv[3].mul_scalar(matrix[0][3]));

    let w1 = cv[0]
        .mul_scalar(matrix[1][0])
        .add(&cv[1].mul_scalar(matrix[1][1]))
        .add(&cv[2].mul_scalar(matrix[1][2]))
        .add(&cv[3].mul_scalar(matrix[1][3]));

    let w2 = cv[0]
        .mul_scalar(matrix[2][0])
        .add(&cv[1].mul_scalar(matrix[2][1]))
        .add(&cv[2].mul_scalar(matrix[2][2]))
        .add(&cv[3].mul_scalar(matrix[2][3]));

    let w3 = cv[0]
        .mul_scalar(matrix[3][0])
        .add(&cv[1].mul_scalar(matrix[3][1]))
        .add(&cv[2].mul_scalar(matrix[3][2]))
        .add(&cv[3].mul_scalar(matrix[3][3]));

    // Horner's method: ((w0 * x + w1) * x + w2) * x + w3
    w0.mul_scalar(x)
        .add(&w1)
        .mul_scalar(x)
        .add(&w2)
        .mul_scalar(x)
        .add(&w3)
}

fn main() {
    println!("=== Minimal Trait Bounds Demo ===\n");

    // Example with our minimal trait
    let cv = [
        Point2D { x: 0.0, y: 0.0 },
        Point2D { x: 1.0, y: 2.0 },
        Point2D { x: 2.0, y: 1.0 },
        Point2D { x: 3.0, y: 0.0 },
    ];

    // Catmull-Rom matrix
    let matrix = [
        [-0.5, 1.5, -1.5, 0.5],
        [1.0, -2.5, 2.0, -0.5],
        [-0.5, 0.0, 0.5, 0.0],
        [0.0, 1.0, 0.0, 0.0],
    ];

    let result = spline_segment_minimal(0.5, &cv, matrix);
    println!(
        "Interpolated point at t=0.5: ({:.2}, {:.2})",
        result.x, result.y
    );

    println!("\nCurrent trait requirements:");
    println!("- Full implementation: Float + FromPrimitive + Zero + One + Clone + Lerp");
    println!("- Actually needed for spline_segment: Clone + Mul<Scalar> + Add");
    println!("\nThe Float requirement excludes vector types unnecessarily!");
}
