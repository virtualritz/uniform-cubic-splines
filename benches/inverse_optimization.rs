//! Benchmark to compare linear vs binary search for spline_inverse
//! This helps identify the sweet spot where binary search becomes beneficial

use std::time::Instant;
use uniform_cubic_splines::{spline_inverse, CatmullRom};

const ITERATIONS: u32 = 100_000;

fn benchmark_inverse(
    knots: &[f64],
    y_values: &[f64],
    description: &str,
) -> f64 {
    let start = Instant::now();
    let mut sum = 0.0;

    for &y in y_values {
        for _ in 0..ITERATIONS {
            if let Some(x) = spline_inverse::<CatmullRom, _>(y, knots) {
                sum += x;
            }
        }
    }

    let elapsed = start.elapsed();
    let ns_per_op =
        elapsed.as_nanos() as f64 / (ITERATIONS as f64 * y_values.len() as f64);

    println!("{:<30} {:8.1} ns/op", description, ns_per_op);
    std::hint::black_box(sum);
    ns_per_op
}

fn main() {
    println!("=== spline_inverse Performance vs Segment Count ===");
    println!("Current implementation (linear search through segments)\n");

    // Test with different numbers of segments
    // CatmullRom: number_of_segments = (knots.len() - 4) + 1
    let segment_counts = vec![
        1,    // 4 knots
        2,    // 5 knots
        3,    // 6 knots
        5,    // 8 knots
        10,   // 13 knots
        20,   // 23 knots
        50,   // 53 knots
        100,  // 103 knots
        200,  // 203 knots
        500,  // 503 knots
        1000, // 1003 knots
    ];

    println!("{:<20} {:<15} {:<20}", "Segments", "Knots", "Time (ns)");
    println!("{:-<55}", "");

    for segments in segment_counts {
        let knot_count = segments + 3; // For CatmullRom
        let knots: Vec<f64> = (0..knot_count)
            .map(|i| i as f64 / (knot_count - 1) as f64)
            .collect();

        // Test with values at different positions
        let y_values = vec![
            knots[1] * 0.1 + knots[knot_count - 2] * 0.9, // Near start
            knots[1] * 0.5 + knots[knot_count - 2] * 0.5, // Middle
            knots[1] * 0.9 + knots[knot_count - 2] * 0.1, // Near end
        ];

        let description = format!("{:<20} {:<15}", segments, knot_count);
        benchmark_inverse(&knots, &y_values, &description);
    }

    println!("\n=== Analysis ===");
    println!("The time complexity is O(n) where n is the number of segments.");
    println!("For each segment, the algorithm:");
    println!("1. Extracts 4 control points");
    println!("2. Calls the invert function (iterative root finding)");
    println!("\nWith binary search optimization:");
    println!("- First find the segment: O(log n)");
    println!("- Then invert only that segment: O(1)");
    println!("- Total: O(log n) instead of O(n)");
    println!("\nExpected benefit: Significant for > ~10 segments");
}
