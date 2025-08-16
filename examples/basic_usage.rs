//! Basic usage example for uniform-cubic-splines

use uniform_cubic_splines::prelude::*;

fn main() {
    println!("=== Basic Spline Interpolation Example ===\n");

    // Example 1: Simple 1D spline interpolation
    println!("1D Catmull-Rom spline:");
    let knots = vec![0.0, 0.0, 1.0, 4.0, 3.0, 3.0, 2.0];

    for x in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let y = spline::<CatmullRom, _, _>(x, &knots).unwrap();
        println!("  spline({:.2}) = {:.4}", x, y);
    }

    // Example 2: Spline inversion (finding x for a given y)
    println!("\nSpline inversion:");
    let monotonic_knots = vec![0.0, 0.0, 1.0, 2.0, 4.0, 8.0, 16.0];
    let y_target = 3.0;

    if let Some(x) = spline_inverse::<CatmullRom, _>(y_target, &monotonic_knots)
    {
        println!("  For y = {}, x = {:.4}", y_target, x);

        // Verify by evaluating forward
        let y_check = spline::<CatmullRom, _, _>(x, &monotonic_knots).unwrap();
        println!("  Verification: spline({:.4}) = {:.4}", x, y_check);
    }

    // Example 3: Non-uniform abscissae using combination of inverse and forward
    println!("\nNon-uniform abscissae interpolation:");
    let x_values = vec![0.0, 0.0, 0.1, 0.3, 0.7, 1.0, 1.0];
    let y_values = vec![0.0, 0.0, 1.0, 4.0, 2.0, 3.0, 3.0];

    let x_query = 0.5; // We want the y value at x = 0.5

    // First find the parameter t that corresponds to x = 0.5
    let t = spline_inverse::<CatmullRom, _>(x_query, &x_values).unwrap();
    // Then evaluate the y spline at that parameter
    let y = spline::<CatmullRom, _, _>(t, &y_values).unwrap();

    println!("  At x = {}, y = {:.4}", x_query, y);

    // Example 4: Different basis functions
    println!("\nDifferent basis functions with same knots:");
    let knots = vec![0.0, 1.0, 2.0, 1.0, 0.0, -1.0];
    let x = 0.5;

    println!(
        "  Linear:      {:.4}",
        spline::<Linear, _, _>(x, &knots).unwrap()
    );
    println!(
        "  Bspline:     {:.4}",
        spline::<Bspline, _, _>(x, &knots).unwrap()
    );
    println!(
        "  CatmullRom:  {:.4}",
        spline::<CatmullRom, _, _>(x, &knots).unwrap()
    );

    // Bezier requires specific knot count (4n+3)
    let bezier_knots = vec![0.0, 1.0, 2.0, 1.0, 0.0, -1.0, -2.0];
    println!(
        "  Bezier:      {:.4}",
        spline::<Bezier, _, _>(x, &bezier_knots).unwrap()
    );

    // Hermite requires specific knot count (4n+2)
    let hermite_knots = vec![0.0, 1.0, 2.0, 1.0, 0.0, -1.0];
    println!(
        "  Hermite:     {:.4}",
        spline::<Hermite, _, _>(x, &hermite_knots).unwrap()
    );

    // Example 5: Error handling
    println!("\nError handling:");
    let too_few_knots = vec![1.0, 2.0]; // Need at least 4 knots

    match spline::<CatmullRom, _, _>(0.5, &too_few_knots) {
        Ok(y) => println!("  Result: {}", y),
        Err(e) => println!("  Error: {}", e),
    }
}
