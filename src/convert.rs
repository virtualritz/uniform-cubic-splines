//! Basis conversion for control values.
//!
//! Uniform cubic splines expressed in different bases represent the same
//! underlying cubic polynomial with different control values (CVs). This
//! module provides functions to convert a set of CVs from one [`Basis`] to
//! another without changing the curve that is evaluated.
//!
//! The math is straightforward: evaluating a cubic spline segment computes
//! `y(t) = [t³ t² t 1] · (M · P)`, where `M` is the basis matrix and `P` is
//! the 4-element column vector of control values. For two bases to describe
//! the same curve, their product `M · P` must be equal, so the conversion
//! from an input basis `IN` to an output basis `OUT` is given by
//! `P_out = M_out⁻¹ · M_in · P_in`.
//!
//! # Non-invertible output bases
//!
//! The [`Linear`](crate::basis::Linear) basis matrix has rank 2 and is not
//! invertible. Using `Linear` as the *output* basis of these functions will
//! panic. `Linear` may still be used as the *input* basis: the resulting
//! cubic polynomial is linear, which is expressible in any of the other
//! bases exactly.

use crate::Basis;
use core::ops::{Add, Div, Mul, Sub};
use num_traits::{One, Zero};

#[cfg(feature = "alloc")]
extern crate alloc;
#[cfg(feature = "alloc")]
use alloc::vec::Vec;

/// Converts the four control values of a single spline segment from basis
/// `BasisIn` to basis `BasisOut`.
///
/// The output control values describe the exact same cubic polynomial as
/// the input, i.e. evaluating the returned CVs with
/// [`spline_segment()`](crate::spline_segment()) in `BasisOut` yields identical
/// values to evaluating `cv` in `BasisIn` for every parameter value.
///
/// # Panics
///
/// Panics if `BasisOut` is not invertible. Among the bases shipped with this
/// crate, only [`Linear`](crate::basis::Linear) has a non-invertible matrix
/// and so cannot be used as the output basis.
///
/// # Examples
///
/// ```
/// use uniform_cubic_splines::{
///     basis::{Bezier, CatmullRom},
///     convert::convert_segment,
///     spline_segment,
/// };
///
/// let cv_cr = [0.0_f64, 1.0, 2.0, 1.5];
/// let cv_bezier = convert_segment::<CatmullRom, Bezier, _>(&cv_cr);
///
/// // Same curve, different representation.
/// let t = 0.37;
/// let y_cr = spline_segment::<CatmullRom, _, _>(t, &cv_cr);
/// let y_bz = spline_segment::<Bezier, _, _>(t, &cv_bezier);
/// assert!((y_cr - y_bz).abs() < 1e-12);
/// ```
pub fn convert_segment<BasisIn, BasisOut, T>(cv: &[T; 4]) -> [T; 4]
where
    BasisIn: Basis<T>,
    BasisOut: Basis<T>,
    T: Clone
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>,
{
    let conversion = conversion_matrix::<BasisIn, BasisOut, T>();
    apply_matrix(&conversion, cv)
}

/// Converts a full knot vector of control values from basis `BasisIn` to basis
/// `BasisOut`, returning a freshly allocated [`Vec`].
///
/// The input must be a valid knot vector for `BasisIn`
/// (see [`Basis::is_len_ok()`]); the output will be a valid knot vector
/// for `BasisOut` that describes the same piecewise cubic curve. For `k` input
/// segments, the returned [`Vec`] has length `4 + (k - 1) · BasisOut::STEP`.
///
/// This function is only available with the `alloc` cargo feature enabled.
///
/// # Panics
///
/// Panics if `cv.len()` is not a valid knot length for `BasisIn`, or if
/// `BasisOut` is not invertible (see [`convert_segment()`]).
///
/// # Examples
///
/// ```
/// # #[cfg(feature = "alloc")] {
/// use uniform_cubic_splines::{
///     basis::{Basis, Bezier, CatmullRom},
///     convert::convert_spline,
/// };
///
/// let cvs = [0.0_f64, 1.0, 2.0, 1.5, 0.5, 0.0];
/// let out = convert_spline::<CatmullRom, Bezier, _>(&cvs);
/// assert!(<Bezier as Basis<f64>>::is_len_ok(out.len()));
/// # }
/// ```
#[cfg(feature = "alloc")]
pub fn convert_spline<BasisIn, BasisOut, T>(cv: &[T]) -> Vec<T>
where
    BasisIn: Basis<T>,
    BasisOut: Basis<T>,
    T: Clone
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>,
{
    assert!(
        BasisIn::is_len_ok(cv.len()),
        "input slice is not a valid knot vector for the input basis"
    );

    let conversion = conversion_matrix::<BasisIn, BasisOut, T>();
    let step_in = BasisIn::STEP;
    let step_out = BasisOut::STEP;
    let number_of_segments = (cv.len() - 4) / step_in + 1;

    let output_len = 4 + (number_of_segments - 1) * step_out;
    let mut out: Vec<T> = Vec::with_capacity(output_len);

    (0..number_of_segments).for_each(|s| {
        let start = s * step_in;
        let window: [T; 4] = [
            cv[start].clone(),
            cv[start + 1].clone(),
            cv[start + 2].clone(),
            cv[start + 3].clone(),
        ];
        let converted = apply_matrix(&conversion, &window);

        // First segment contributes all four CVs; subsequent segments
        // contribute the last `step_out` of their converted window so the
        // resulting knot vector matches `BasisOut`'s stride and the adjacent
        // windows share their boundary CVs by construction.
        if s == 0 {
            out.extend(converted);
        } else {
            out.extend(converted.into_iter().skip(4 - step_out));
        }
    });

    out
}

/// Computes `C = M_out⁻¹ · M_in`, panicking if `M_out` is not invertible.
fn conversion_matrix<BasisIn, BasisOut, T>() -> [[T; 4]; 4]
where
    BasisIn: Basis<T>,
    BasisOut: Basis<T>,
    T: Clone
        + Zero
        + One
        + Add<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Div<Output = T>,
{
    let m_out_inv = invert_4x4(BasisOut::MATRIX)
        .expect("output basis matrix is not invertible");
    mul_4x4(&m_out_inv, &BasisIn::MATRIX)
}

/// Applies a 4x4 matrix of scalars to a column vector of four CVs.
#[inline]
fn apply_matrix<T>(m: &[[T; 4]; 4], v: &[T; 4]) -> [T; 4]
where
    T: Clone + Add<Output = T> + Mul<Output = T>,
{
    [
        v[0].clone() * m[0][0].clone()
            + v[1].clone() * m[0][1].clone()
            + v[2].clone() * m[0][2].clone()
            + v[3].clone() * m[0][3].clone(),
        v[0].clone() * m[1][0].clone()
            + v[1].clone() * m[1][1].clone()
            + v[2].clone() * m[1][2].clone()
            + v[3].clone() * m[1][3].clone(),
        v[0].clone() * m[2][0].clone()
            + v[1].clone() * m[2][1].clone()
            + v[2].clone() * m[2][2].clone()
            + v[3].clone() * m[2][3].clone(),
        v[0].clone() * m[3][0].clone()
            + v[1].clone() * m[3][1].clone()
            + v[2].clone() * m[3][2].clone()
            + v[3].clone() * m[3][3].clone(),
    ]
}

/// Computes `a · b` for two 4x4 matrices.
fn mul_4x4<T>(a: &[[T; 4]; 4], b: &[[T; 4]; 4]) -> [[T; 4]; 4]
where
    T: Clone + Zero + Add<Output = T> + Mul<Output = T>,
{
    let row = |i: usize| -> [T; 4] {
        let cell = |j: usize| {
            a[i][0].clone() * b[0][j].clone()
                + a[i][1].clone() * b[1][j].clone()
                + a[i][2].clone() * b[2][j].clone()
                + a[i][3].clone() * b[3][j].clone()
        };
        [cell(0), cell(1), cell(2), cell(3)]
    };
    [row(0), row(1), row(2), row(3)]
}

/// Inverts a 4x4 matrix via Gauss-Jordan elimination with a simple row-swap
/// fallback when a zero pivot is encountered. Returns `None` for singular
/// matrices.
fn invert_4x4<T>(m: [[T; 4]; 4]) -> Option<[[T; 4]; 4]>
where
    T: Clone + Zero + One + Sub<Output = T> + Mul<Output = T> + Div<Output = T>,
{
    // Augmented matrix `[M | I]`, stored as 4 rows of 8 columns.
    let mut a: [[T; 8]; 4] = [
        [
            m[0][0].clone(),
            m[0][1].clone(),
            m[0][2].clone(),
            m[0][3].clone(),
            T::one(),
            T::zero(),
            T::zero(),
            T::zero(),
        ],
        [
            m[1][0].clone(),
            m[1][1].clone(),
            m[1][2].clone(),
            m[1][3].clone(),
            T::zero(),
            T::one(),
            T::zero(),
            T::zero(),
        ],
        [
            m[2][0].clone(),
            m[2][1].clone(),
            m[2][2].clone(),
            m[2][3].clone(),
            T::zero(),
            T::zero(),
            T::one(),
            T::zero(),
        ],
        [
            m[3][0].clone(),
            m[3][1].clone(),
            m[3][2].clone(),
            m[3][3].clone(),
            T::zero(),
            T::zero(),
            T::zero(),
            T::one(),
        ],
    ];

    for i in 0..4 {
        if a[i][i].is_zero() {
            // Find a row below with a non-zero entry in column `i`.
            let swap_row = ((i + 1)..4).find(|&k| !a[k][i].is_zero());
            if let Some(k) = swap_row {
                a.swap(i, k);
            } else {
                return None;
            }
        }

        let pivot = a[i][i].clone();
        (0..8).for_each(|j| {
            a[i][j] = a[i][j].clone() / pivot.clone();
        });

        (0..4).for_each(|k| {
            if k != i && !a[k][i].is_zero() {
                let factor = a[k][i].clone();
                (0..8).for_each(|j| {
                    a[k][j] =
                        a[k][j].clone() - factor.clone() * a[i][j].clone();
                });
            }
        });
    }

    Some([
        [
            a[0][4].clone(),
            a[0][5].clone(),
            a[0][6].clone(),
            a[0][7].clone(),
        ],
        [
            a[1][4].clone(),
            a[1][5].clone(),
            a[1][6].clone(),
            a[1][7].clone(),
        ],
        [
            a[2][4].clone(),
            a[2][5].clone(),
            a[2][6].clone(),
            a[2][7].clone(),
        ],
        [
            a[3][4].clone(),
            a[3][5].clone(),
            a[3][6].clone(),
            a[3][7].clone(),
        ],
    ])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::{Bezier, Bspline, CatmullRom, Hermite, Linear, Power};
    use crate::spline_segment;

    fn sample_points() -> [f64; 7] {
        [0.0, 0.125, 0.25, 0.5, 0.625, 0.875, 1.0]
    }

    fn assert_segment_equivalent<BasisIn, BasisOut>(cv_in: &[f64; 4])
    where
        BasisIn: Basis<f64>,
        BasisOut: Basis<f64>,
    {
        let cv_out = convert_segment::<BasisIn, BasisOut, f64>(cv_in);
        for t in sample_points() {
            let y_in = spline_segment::<BasisIn, _, _>(t, cv_in);
            let y_out = spline_segment::<BasisOut, _, _>(t, &cv_out);
            assert!(
                (y_in - y_out).abs() < 1e-10,
                "mismatch at t={t}: {y_in} vs {y_out}"
            );
        }
    }

    #[test]
    fn segment_identity_bspline() {
        let cv = [1.0, 2.0, 3.0, 5.0];
        assert_segment_equivalent::<Bspline, Bspline>(&cv);
    }

    #[test]
    fn segment_bspline_to_bezier_roundtrip() {
        let cv = [1.0, 2.0, 3.0, 5.0];
        let bz = convert_segment::<Bspline, Bezier, f64>(&cv);
        let bs = convert_segment::<Bezier, Bspline, f64>(&bz);
        for (a, b) in cv.iter().zip(bs.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn segment_cross_basis_equivalence() {
        let cv = [0.5, -1.25, 2.0, 0.75];
        assert_segment_equivalent::<Bspline, Bezier>(&cv);
        assert_segment_equivalent::<Bspline, CatmullRom>(&cv);
        assert_segment_equivalent::<Bspline, Hermite>(&cv);
        assert_segment_equivalent::<Bspline, Power>(&cv);

        assert_segment_equivalent::<CatmullRom, Bezier>(&cv);
        assert_segment_equivalent::<CatmullRom, Hermite>(&cv);
        assert_segment_equivalent::<CatmullRom, Bspline>(&cv);

        assert_segment_equivalent::<Bezier, Hermite>(&cv);
        assert_segment_equivalent::<Hermite, Bezier>(&cv);
        assert_segment_equivalent::<Power, Bezier>(&cv);
    }

    #[test]
    fn segment_linear_as_input_is_ok() {
        // Linear ignores cv[0] and cv[3]; the resulting curve is the
        // straight line between cv[1] and cv[2]. Any invertible basis
        // must be able to represent that exactly.
        let cv = [99.0, 1.0, 3.0, -99.0];
        assert_segment_equivalent::<Linear, Bezier>(&cv);
        assert_segment_equivalent::<Linear, Bspline>(&cv);
        assert_segment_equivalent::<Linear, Power>(&cv);
    }

    #[test]
    #[should_panic(expected = "not invertible")]
    fn segment_linear_as_output_panics() {
        let cv = [0.0, 1.0, 2.0, 3.0];
        let _ = convert_segment::<Bspline, Linear, f64>(&cv);
    }

    #[test]
    fn invert_4x4_identity() {
        let id: [[f64; 4]; 4] = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let inv = invert_4x4(id).unwrap();
        inv.iter().enumerate().for_each(|(i, row)| {
            row.iter().enumerate().for_each(|(j, &value)| {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((value - expected).abs() < 1e-12);
            });
        });
    }

    #[test]
    fn invert_4x4_linear_returns_none() {
        assert!(invert_4x4(<Linear as Basis<f64>>::MATRIX).is_none());
    }

    #[cfg(feature = "alloc")]
    mod alloc_tests {
        use super::*;
        use alloc::vec;

        fn assert_curve_equivalent<BasisIn, BasisOut>(cv_in: &[f64])
        where
            BasisIn: Basis<f64>,
            BasisOut: Basis<f64>,
        {
            let cv_out = convert_spline::<BasisIn, BasisOut, f64>(cv_in);
            assert!(
                BasisOut::is_len_ok(cv_out.len()),
                "converted length {} not valid for output basis",
                cv_out.len()
            );
            for i in 0..=20 {
                let x = i as f64 / 20.0;
                let y_in = crate::spline::<BasisIn, _, _>(x, cv_in).unwrap();
                let y_out =
                    crate::spline::<BasisOut, _, _>(x, &cv_out).unwrap();
                assert!(
                    (y_in - y_out).abs() < 1e-9,
                    "mismatch at x={x}: {y_in} vs {y_out}"
                );
            }
        }

        #[test]
        fn spline_bspline_to_bezier_multi_segment() {
            let cv = vec![0.0, 1.0, 2.0, 3.0, 2.5, 1.5, 0.5];
            assert_curve_equivalent::<Bspline, Bezier>(&cv);
        }

        #[test]
        fn spline_catmull_rom_to_hermite() {
            let cv = vec![0.0, 0.5, 1.5, 3.0, 2.0, 1.0, 0.0, -0.5];
            assert_curve_equivalent::<CatmullRom, Hermite>(&cv);
        }

        #[test]
        fn spline_bezier_to_power() {
            let cv = vec![0.0, 1.0, 2.0, 3.0, 2.5, 1.5, 0.5];
            assert_curve_equivalent::<Bezier, Power>(&cv);
        }

        #[test]
        fn spline_lengths_match_output_basis() {
            let cv = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
            let out_cr = convert_spline::<Bspline, CatmullRom, f64>(&cv);
            assert_eq!(out_cr.len(), cv.len());

            let out_bz = convert_spline::<Bspline, Bezier, f64>(&cv);
            assert_eq!(out_bz.len(), 4 + (cv.len() - 4) * 3);

            let out_her = convert_spline::<Bspline, Hermite, f64>(&cv);
            assert_eq!(out_her.len(), 4 + (cv.len() - 4) * 2);

            let out_pow = convert_spline::<Bspline, Power, f64>(&cv);
            assert_eq!(out_pow.len(), 4 + (cv.len() - 4) * 4);
        }
    }
}
