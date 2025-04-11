use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use uniform_cubic_splines::{spline, spline_inverse, CatmullRom};

fn generate_random_values(n: usize, seed: u64) -> Vec<f64> {
    let mut rng = StdRng::seed_from_u64(seed);

    // Generate n random values
    (0..n).map(|_| rng.random_range(-10.0..10.0)).collect()
}

fn generate_random_monotonic_values(n: usize, seed: u64) -> Vec<f64> {
    let mut rng = StdRng::seed_from_u64(seed);

    // Generate n random values that are monotonically increasing
    let mut current = 0.0;
    (0..n)
        .map(|_| {
            current += rng.random_range(0.1..1.0); // Ensure values are increasing
            current
        })
        .collect()
}

fn benchmark_spline(c: &mut Criterion) {
    let mut group = c.benchmark_group("Spline Functions");

    // Fixed seed for deterministic results
    const SEED: u64 = 12345;
    // Number of iterations to prepare data for
    const ITERATIONS: usize = 10000;

    // Test with different numbers of knots
    for num_points in [10, 50, 100, 500] {
        // Generate different sets of random knots for spline
        let all_spline_knots: Vec<f64> = (0..ITERATIONS)
            .flat_map(|i| generate_random_values(num_points, SEED + i as u64))
            .collect();

        // Generate different sets of monotonic knots for spline_inverse
        let all_inverse_knots: Vec<f64> = (0..ITERATIONS)
            .flat_map(|i| {
                generate_random_monotonic_values(
                    num_points,
                    SEED + 10000 + i as u64,
                )
            })
            .collect();

        // Generate x test inputs for spline
        let x_inputs: Vec<f64> = (0..ITERATIONS)
            .map(|i| {
                let x_max = (num_points - 1) as f64;
                let mut rng = StdRng::seed_from_u64(SEED + 20000 + i as u64);
                rng.random_range(0.0..x_max)
            })
            .collect();

        // Benchmark spline function with CatmullRom basis
        group.bench_function(BenchmarkId::new("spline", num_points), |b| {
            let mut local_rng = StdRng::seed_from_u64(SEED + 40000);
            b.iter_with_large_drop(|| {
                let i = local_rng.random_range(0..ITERATIONS);
                let knots =
                    &all_spline_knots[i * num_points..(i + 1) * num_points];
                let x = x_inputs[i];
                black_box(spline::<CatmullRom, _, _>(x, knots))
            });
        });

        // Generate y test inputs for spline_inverse
        let y_inputs: Vec<f64> = (0..ITERATIONS)
            .enumerate()
            .map(|(i, iter_idx)| {
                let knot_start = i * num_points;
                let y_min = all_inverse_knots[knot_start];
                let y_max = all_inverse_knots[knot_start + num_points - 1];
                let mut rng =
                    StdRng::seed_from_u64(SEED + 30000 + iter_idx as u64);
                rng.random_range(y_min..y_max)
            })
            .collect();

        // Benchmark spline_inverse function with CatmullRom basis
        group.bench_function(
            BenchmarkId::new("spline_inverse", num_points),
            |b| {
                let mut local_rng = StdRng::seed_from_u64(SEED + 50000);
                b.iter_with_large_drop(|| {
                    let i = local_rng.random_range(0..ITERATIONS);
                    let knots = &all_inverse_knots
                        [i * num_points..(i + 1) * num_points];
                    let y = y_inputs[i];
                    black_box(spline_inverse::<CatmullRom, _>(y, knots))
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, benchmark_spline);
criterion_main!(benches);
