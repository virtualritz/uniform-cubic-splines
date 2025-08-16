macro_rules! b_spline_basis {
    ($type:ident) => {
        impl Basis<$type> for Bspline {
            const NAME: &'static str = "B-spline";
            const STEP: usize = 1;
            const MATRIX: [[$type; 4]; 4] = [
                [-1.0 / 6.0, 3.0 / 6.0, -3.0 / 6.0, 1.0 / 6.0],
                [3.0 / 6.0, -6.0 / 6.0, 3.0 / 6.0, 0.0 / 6.0],
                [-3.0 / 6.0, 0.0 / 6.0, 3.0 / 6.0, 0.0 / 6.0],
                [1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0, 0.0 / 6.0],
            ];
        }
    };
}

macro_rules! bezier_basis {
    ($type:ident) => {
        impl Basis<$type> for Bezier {
            const NAME: &'static str = "Bezier";
            const STEP: usize = 3;
            const MATRIX: [[$type; 4]; 4] = [
                [-1.0, 3.0, -3.0, 1.0],
                [3.0, -6.0, 3.0, 0.0],
                [-3.0, 3.0, 0.0, 0.0],
                [1.0, 0.0, 0.0, 0.0],
            ];
        }
    };
}

macro_rules! catmull_rom_basis {
    ($type:ident) => {
        impl Basis<$type> for CatmullRom {
            const NAME: &'static str = "Catmull-Rom";
            const STEP: usize = 1;
            const MATRIX: [[$type; 4]; 4] = [
                [-1.0 / 2.0, 3.0 / 2.0, -3.0 / 2.0, 1.0 / 2.0],
                [2.0 / 2.0, -5.0 / 2.0, 4.0 / 2.0, -1.0 / 2.0],
                [-1.0 / 2.0, 0.0 / 2.0, 1.0 / 2.0, 0.0 / 2.0],
                [0.0 / 2.0, 2.0 / 2.0, 0.0 / 2.0, 0.0 / 2.0],
            ];
        }
    };
}

macro_rules! hermite_basis {
    ($type:ident) => {
        impl Basis<$type> for Hermite {
            const NAME: &'static str = "Hermite";
            const STEP: usize = 2;
            const MATRIX: [[$type; 4]; 4] = [
                [2., 1., -2., 1.],
                [-3., -2., 3., -1.],
                [0., 1., 0., 0.],
                [1., 0., 0., 0.],
            ];
        }
    };
}

macro_rules! linear_basis {
    ($type:ident) => {
        impl Basis<$type> for Linear {
            const NAME: &'static str = "Linear";
            const STEP: usize = 1;
            const MATRIX: [[$type; 4]; 4] = [
                [0., 0., 0., 0.],
                [0., 0., 0., 0.],
                [0., -1., 1., 0.],
                [0., 1., 0., 0.],
            ];
        }
    };
}

macro_rules! power_basis {
    ($type:ident) => {
        impl Basis<$type> for Power {
            const NAME: &'static str = "Power";
            const STEP: usize = 4;
            const MATRIX: [[$type; 4]; 4] = [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ];
        }
    };
}

macro_rules! _cardinal_basis {
    ($type:ident, $t:expr) => {
        impl Basis<$type> for CatmullRom {
            const NAME: &'static str = "Cardinal";
            const STEP: usize = 1;
            const MATRIX: [[$type; 4]; 4] = [
                [$t, 2.0 - $t, $t - 2, $t],
                [2 * $t, $t - 3.0, 3 - 2 * $t, -$t],
                [-$t, 0.0, $t, 0.0],
                [0.0, 1.0, 0.0, 0.0],
            ];
        }
    };
}
