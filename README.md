# `uniform-cubic-splines`

Uniform cubic spline interpolation & inversion.

[![Documentation](https://docs.rs/uniform-cubic-splines/badge.svg)](https://docs.rs/uniform-cubic-splines/)
[![Crate](https://img.shields.io/crates/v/uniform-cubic-splines.svg)](https://crates.io/crates/uniform-cubic-splines)

This crate supports the following types of splines:

- [B-spline](https://en.wikipedia.org/wiki/B-spline)
- [Bezier](https://en.wikipedia.org/wiki/Composite_B%C3%A9zier_curve)
- [Catmull-Rom](https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Catmull%E2%80%93Rom_spline)
- [Hermite](https://en.wikipedia.org/wiki/Cubic_Hermite_spline)
- Linear
- Power

![Curve widget with 1D Catmull-Rom spline](spline_ui.png)

_Curve widget with a 1D Catmull-Rom spline with non-uniform knot
spacing and knot multiplicity using this crate for interpolation
(drawn using `tiny-skia`)._

The crate uses generics to allow interpolation of any type for which
certain traits are defined.

I.e. you can use this crate to interpolate splines in 1D, 2D, 3D, etc.

```toml
[dependencies]
uniform-cubic-splines = { version = "0.3" }
```

## Example

Using a combination of `spline_inverse()` and `spline()` it is
possible to compute a full spline-with-nonuniform-abscissæ:

```rust
use uniform_cubic_splines::{
    basis::CatmullRom, spline_inverse, spline,
};

// We want to evaluate the spline at knot value 0.3.
let x = 0.3;

// The first and last points are never interpolated.
let knot_spacing = [0.0, 0.0, 0.1, 0.3, 1.0, 1.0];
let knots        = [0.0, 0.0, 1.3, 4.2, 3.2, 3.2];

let v = spline_inverse::<CatmullRom, _>(x, &knot_spacing).unwrap();
let y = spline::<CatmullRom, _, _>(v, &knots);

assert!(y - 4.2 < 1e-6);
```

## Features

- `monotonic_check` -- The
  [`spline_inverse()`](https://docs.rs/uniform-cubic-splines/latest/uniform_cubic_splines/fn.spline_inverse.html)/`spline_inverse_with()`
  code will check if the knot vector is monotonic (_enabled_ by default).

## `f16` & `f128` Support

This crate supports `f16` and `f128` types on a `nightly` toolchain if you use this repository as an overlay in your `Cargo.toml`.

```toml
[patch.crates-io]
uniform-cubic-splines = {
    git = "https://github.com/virtualritz/uniform-cubic-splines.git"
}
```

> This will be supported without an overlay once [this PR on `num-traits`](https://github.com/rust-num/num-traits/pull/333) is merged and published.

## Background

The code is a Rust port of the implementation found in the [Open
Shading Language](https://github.com/imageworks/OpenShadingLanguage)
C++ source.

If you come from a background of computer graphics/shading
languages used in offline rendering this crate should feel like
home.

## Speed

The code was originally a faithful but idiomatic Rust port of the C++ source. It was optimized quite a bit afterwards. The optimized code is in version `0.3.3` of the crate and after.

For perspective, the current version's `spline()` function is about twice as fast and the `spline_inverse()` up to 3.5 times as fast as the original, faithful Rust port.

### Benchmarks

| **Points** | `spline()` | `spline_inverse()` |
| ---------- | ---------- | ------------------ |
| 10         | 15.5 ns    | 95 ns              |
| 50         | 16.4 ns    | 181.2 ns           |
| 100        | 16.8 ns    | 277.7 ns           |
| 500        | 20.6 ns    | 1.275 µs           |

These were taken on Ubuntu 25.04 with `rustc 1.89.0-nightly` and `target_cpu = "native"` on an AMD Ryzen 7 6800H laptop.
