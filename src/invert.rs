use roots::{find_root_brent, FloatType, SimpleConvergency};
use lerp::Lerp;

/// Evaluate the inverse of a spline; i.e., solve for
/// the x for which [`spline`]`(x) == y`.
pub fn spline_inverse<B, T>(y: T, knots: &[T]) -> Option<T>
where
    B: Basis<T>,
    T: AsPrimitive<usize>
        + FloatType
        + Float
        + FromPrimitive
        + PartialOrd
        + One
        + Zero,
{
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
            return Some(num_traits::Zero::zero());
        }
        if y >= knots[high_index] {
            return Some(num_traits::One::one());
        }
    } else {
        if y >= knots[low_index] {
            return Some(num_traits::Zero::zero());
        }
        if y <= knots[high_index] {
            return Some(num_traits::One::one());
        }
    }

    let spline_function = |x| spline::<B, T, T>(x, knots);

    let mut convergency = SimpleConvergency {
        eps: T::from_f64(1e-15).unwrap(),
        max_iter: 32,
    };

    match find_root_brent(
        y,
        num_traits::Zero::zero(),
        &spline_function,
        &mut convergency,
    ) {
        Ok(x) => Some(x),
        Err(e) => {
            println!("{:?}", e);
            None
        }
    }
}

/// Evaluate the inverse of a spline; i.e., solve for
/// the x for which [`spline`]`(x) == y`.
pub fn spline_inverse2<B, T>(y: T, knots: &[T]) -> Option<T>
where
    B: Basis<T>,
    T: AsPrimitive<usize>
        + FloatType
        + Float
        + FromPrimitive
        + PartialOrd
        + One
        + Zero,
{
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
            return Some(num_traits::Zero::zero());
        }
        if y >= knots[high_index] {
            return Some(num_traits::One::one());
        }
    } else {
        if y >= knots[low_index] {
            return Some(num_traits::Zero::zero());
        }
        if y <= knots[high_index] {
            return Some(num_traits::One::one());
        }
    }

    let spline_function = |x| spline::<B, T, T>(x, knots);

    let number_of_segments = (knots.len() - 4) / B::STEP + 1;
    let number_of_segments_inverted = 1.0 / number_of_segments as f64;
    let mut r0 = num_traits::Zero::zero();

    // Search each interval.
    for s in 0..number_of_segments {
        let r1 =
            T::from_f64(number_of_segments_inverted * (s + 1) as f64).unwrap();
        if let Some(x) = invert(
            &spline_function,
            &y,
            &r0,
            &r1,
            32,
            &T::from_f64(1.0e-6).unwrap(),
        ) {
            return Some(x);
        }
        // Start of next interval is end of this one.
        r0 = r1;
    }
    None
}

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

fn invert<T>(
    function: &dyn Fn(T) -> T,
    y: &T,
    x_min: &T,
    x_max: &T,
    max_iterations: usize,
    epsilon: &T,
) -> Option<T>
where
    T: AsPrimitive<usize>
        + FloatType
        + Float
        + FromPrimitive
        + PartialOrd
        + One
        + Zero,
{
    // Use the Regula Falsi method, falling back to bisection if it
    // hasn't converged after 3/4 of the maximum number of iterations.
    // See, e.g., Numerical Recipes for the basic ideas behind both
    // methods.
    let mut v0 = function(*x_min);
    let mut v1 = function(*x_max);

    let mut x = *x_min;
    let increasing = v0 < v1;

    let vmin = if increasing { v0 } else { v1 };
    let vmax = if increasing { v1 } else { v0 };

    if !(*y >= vmin && *y <= vmax) {
        return None;
    }

    // Already close enough.
    if Float::abs(v0 - v1) < *epsilon {
        println!("Bye");
        return Some(x);
    }

    // How many times to try regula falsi.
    let rfiters = (3 * max_iterations) / 4;

    let mut x_min = *x_min;
    let mut x_max = *x_max;

    for iters in 0..max_iterations {
        // Interpolation factor.
        let mut t: T;
        if iters < rfiters {
            // Regula falsi.
            t = (*y - v0) / (v1 - v0);
            if t <= num_traits::Zero::zero() || t >= num_traits::One::one() {
                // RF convergence failure -- bisect instead.
                t = T::from_f64(0.5).unwrap() * num_traits::One::one();
            }
        } else {
            // Bisection.
            t = T::from_f64(0.5).unwrap() * num_traits::One::one();
        }
        x = t.lerp(x_min, x_max);
        let v = function(x);
        if (v < *y) == increasing {
            x_min = x;
            v0 = v;
        } else {
            x_max = x;
            v1 = v;
        }
        if Float::abs(x_max - x_min) < *epsilon || Float::abs(v - *y) < *epsilon
        {
            return Some(x); // converged
        }
    }
    Some(x)
}

#[test]

fn test() {
    let knots = vec![-0.4, 0.0, 0.4, 0.5, 0.9, 1.0, 1.9];

    assert!(is_len_ok::<CatmullRom>(knots.len()));
    assert!(0.4 == spline::<CatmullRom, _, _>(0.25f64, &knots));

    let knots = [0.0, 0.0, 0.5, 0.5];

    println!("{:?}", spline_inverse2::<Linear, _>(0.25f64, &knots));

    let f = |x| spline::<Linear, _, _>(x, &knots);

    let mut convergency = SimpleConvergency {
        eps: 1e-15f64,
        max_iter: 30,
    };

    println!("{:?}", find_root_brent(0.5f64, 0f64, &f, &mut convergency));
    // Returns approximately Ok(1);
}
