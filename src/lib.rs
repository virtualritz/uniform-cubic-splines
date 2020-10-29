#![feature(clamp)]
use core::marker::PhantomData;
use core::ops::Add;
use core::ops::Mul;
use num::{
    cast::AsPrimitive,
    traits::identities::{One, Zero}, // , ops::mul_add::MulAdd},
    Float,
    FromPrimitive,
};
use roots::{find_root_brent, FloatType, SimpleConvergency};


pub trait Basis<T: Float> {
    const STEP: usize;
    const MATRIX: [[T; 4]; 4];
}

pub struct CatmullRom;

impl Basis<f32> for CatmullRom {
    const STEP: usize = 1;
    const MATRIX: [[f32; 4]; 4] = [
        [-1.0 / 2.0, 3.0 / 2.0, -3.0 / 2.0, 1.0 / 2.0],
        [2.0 / 2.0, -5.0 / 2.0, 4.0 / 2.0, -1.0 / 2.0],
        [-1.0 / 2.0, 0.0 / 2.0, 1.0 / 2.0, 0.0 / 2.0],
        [0.0 / 2.0, 2.0 / 2.0, 0.0 / 2.0, 0.0 / 2.0],
    ];
}

pub struct Bezier;

impl Basis<f32> for Bezier {
    const STEP: usize = 3;
    const MATRIX: [[f32; 4]; 4] = [
        [-1.0, 3.0, -3.0, 1.0],
        [3.0, -6.0, 3.0, 0.0],
        [-3.0, 3.0, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
    ];
}

pub struct Interpolation<T: Float, BASIS: Basis<T>> {
    _marker: PhantomData<(T, BASIS)>,
}

impl<T, B> Interpolation<T, B>
where
    T: AsPrimitive<usize> + FloatType + Float + FromPrimitive + Ord + One + Zero,
    B: Basis<T>,
{
    pub fn evaluate<U>(x: T, knots: Vec<U>) -> U
    where
        U: Add<Output = U> + Copy + Mul<T, Output = U> + Zero,
    {
        let number_of_segments: usize = ((knots.len() - 4) / B::STEP) + 1;
        let mut x: T = x.clamp(num::Zero::zero(), num::One::one()) * T::from_usize(number_of_segments).unwrap();

        //let seg_x: T = x;
        //let segment = (seg_x.as_()).clamp(0, number_of_segments - 1);
        // TRYME: can we just x.floor() here?
        let segment = x.floor();

        // x is the position along the segment.
        x = x - segment;
        let start: usize = segment.as_() * B::STEP;

        // Get a slice for the segment.
        let cv = &knots[start..start + 4];

        let tk = (0..4)
            .map(|knot| {
                cv.iter()
                    .zip(B::MATRIX[knot].iter())
                    .fold(U::zero(), |total, (cv, basis)| total + *cv * *basis)
            })
            .collect::<Vec<U>>();

        ((tk[0] * x + tk[1]) * x + tk[2]) * x + tk[3]
    }

    /// Evaluate the inverse of a spline; i.e., solve for
    /// the x for which `evaluate(x)` == y.
    pub fn evaluate_inverse<U>(y: T, knots: Vec<T>) -> Option<T> {
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
                return Some(num::Zero::zero());
            }
            if y >= knots[high_index] {
                return Some(num::One::one());
            }
        } else {
            if y >= knots[low_index] {
                return Some(num::Zero::zero());
            }
            if y <= knots[high_index] {
                return Some(num::One::one());
            }
        }

        let spline_function = |x| Self::evaluate(x, knots.clone());

        let mut convergency = SimpleConvergency {
            eps: T::from_f64(1e-6).unwrap(),
            max_iter: 32,
        };

        match find_root_brent(y, num::Zero::zero(), &spline_function, &mut convergency) {
            Ok(x) => {
                Some(x)
            }
            Err(_) => None
        }
    }
}
