## Uniform Cubic Spline Interpolation & Inversion

This crate supports the following types of splines:
* [B-spline](https://en.wikipedia.org/wiki/B-spline)
* [Bezier](https://en.wikipedia.org/wiki/Composite_B%C3%A9zier_curve)
* [Catmull-Rom](https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Catmull%E2%80%93Rom_spline)
* [Hermite](https://en.wikipedia.org/wiki/Cubic_Hermite_spline)
* Linear
* Power

If you come from a background of shading languages used in offline
rendering this crate should feel like home.

The code is a Rust port of the resp. implementation found in the
[Open Shading Language](https://github.com/imageworks/OpenShadingLanguage)
C++ source.

# Example
Using a combination of `spline()` and `spline_inverse()` it is
possible to compute a full spline-with-nonuniform-abscissæ:
```rust
use uniform_cubic_splines::{
    spline, spline_inverse, basis::{CatmullRom, Linear}
};

// We want to evaluate the spline at knot value 0.3.
let x = 0.3;

// The first an last points are never interpolated.
let values = [0.0, 0.0, 1.3, 4.2, 3.2, 3.2];
let knots = [0.0, 0.0, 0.1, 0.3, 1.0, 1.0];

let v = spline_inverse::<Linear, _>(x, &knots).unwrap();
let y = spline::<CatmullRom, _, _>(v, &values);

assert!(y - 4.2 < 1e-6);
```
