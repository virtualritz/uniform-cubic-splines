//! Error types for spline operations.

use thiserror::Error;

/// Errors that can occur during spline operations.
#[derive(Error, Debug)]
pub enum SplineError {
    /// The knot vector has an invalid length for the chosen basis.
    #[error("{basis_name} curve must have at least {min_knots} knots. Found: {actual}")]
    InvalidKnotLength {
        basis_name: &'static str,
        min_knots: usize,
        actual: usize,
    },

    /// The knot vector length doesn't match the required pattern for the
    /// basis.
    #[error("{basis_name} curve requires (length - 4) to be divisible by {extra}. Found length: {actual}")]
    InvalidKnotPattern {
        basis_name: &'static str,
        extra: usize,
        actual: usize,
    },

    /// The knot vector is not monotonic (required for spline inverse).
    #[cfg(feature = "monotonic_check")]
    #[error("The knots array fed to spline_inverse() is not monotonic")]
    NonMonotonicKnots,
}

/// Result type for spline operations.
pub type SplineResult<T> = Result<T, SplineError>;
