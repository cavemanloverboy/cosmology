use std::f64::consts::PI;

use quadrature::Output;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelIterator;

use crate::utils::integration::two_dimensional_integral_const_inner_bounds;

use cdf::calculate_poisson_cdf;

const GRF_INTEGRATION_TARGET_ERROR: f64 = 1e-6;

/// These are in separate modules because they are massive and are
/// stressing VS code being open on a 48 core, 256 GB RAM machine...
mod cdf;
/// These are in separate modules because they are massive and are
/// stressing VS code being open on a 48 core, 256 GB RAM machine...
mod pdf;

/// A helper struct for calculating the k-Nearest Neighbor distributions
/// of gaussian random fields. The intended use is as follows:
///
/// 1) Initialize a [PoissonRandomField] struct with some number density
/// `nbar` and a [SpaceMode]
/// 2) Push distances to the struct (as you would a vector). Internally,
/// this will calcuate I (or its components) and dIdV (or its components).
/// 3) Use the struct's methods to calculate the values of the CDF or PDF at
/// the specified scales.
///
/// TODO: sample usage
pub struct PoissonRandomField {
    /// Average number density
    nbar: f64,

    /// Container for distances (either measurements or scales of interest)
    container: Vec<f64>,
}

impl PoissonRandomField {
    /// Initializes the [PoissonRandomField] struct using `Vec::new` for the underlying vectors.
    pub fn new(nbar: f64) -> Self {
        PoissonRandomField {
            nbar,
            container: Vec::new(),
        }
    }

    // /// Initializes the [PoissonRandomField] struct and passes in the specified capacity
    // /// to the underlying vectors for better, preemtive allocations.
    // pub fn new_with_capacity(
    //     nbar: f64,
    //     mode: SpaceMode<'b, 'c>,
    //     capacity: usize,
    // ) -> Self {
    //     PoissonRandomField {
    //         nbar,
    //         mode,
    //         distances: Vec::with_capacity(capacity)
    //     }
    // }

    /// Initializes the [PoissonRandomField] struct and passes in the specified capacity
    /// to the underlying vectors for better, preemtive allocations.
    pub fn with(mut self, rs: &[f64]) -> Self {
        // Make local copy of distances
        self.container = rs.iter().cloned().collect();

        self
    }

    pub fn get_cdf(&self, k: u8) -> Vec<f64> {
        self.container
            .iter()
            .enumerate()
            .map(|(_, r)| calculate_poisson_cdf(*r, self.nbar, k))
            .collect()
    }
}

fn check_integral(integral: &Output) {
    // TODO
}

/// This weight function is a probe for the number of 2-point configurations
/// with separation `s` in a sphere of radius `r`.
fn spherical_pair_weight_function(s: f64, r: f64) -> f64 {
    PI / 12.0 * (s - 2.0 * r).powi(2) * (s + 4.0 * r)
}

/// This is [spherical_pair_weight_function] differentiated wrt V,
/// still expressed in terms of r via dW/dr * dr/dV = dW/dr * (dV/dr)^-1.
/// TODO: unit tests
fn spherical_pair_weight_function_dv(s: f64, r: f64) -> f64 {
    (2.0 * PI * (-2.0 * r + s).powi(2) - 2.0 * PI * (-2.0 * r + s) * (4.0 * r + s))
        / (4.0 * PI * r.powi(2))
}

#[test]
fn test_spherical_weight_at_diameter() {
    // The volume formed by all two-point configurations where the two points
    // are separated by a sphere's diameter in a sphere is zero.
    assert_approx_eq::assert_approx_eq!(spherical_pair_weight_function(2.0, 1.0), 0.0, 1e-5);
}

#[test]
fn test_spherical_weight_at_radius() {
    // This volume should be the volume of two spherical caps with height = r/2 = 1/2
    // V = pi * h^2 / 3 * (3r - h)
    assert_approx_eq::assert_approx_eq!(
        spherical_pair_weight_function(1.0, 1.0),
        2.0 * PI * 0.25 / 3.0 * (3.0 - 0.5),
        1e-5
    );
}

#[test]
fn test_spherical_weight_at_zero() {
    // This volume should be the volume of the sphere.
    assert_approx_eq::assert_approx_eq!(
        spherical_pair_weight_function(0.0, 1.0),
        4.0 * PI / 3.0,
        1e-5
    );
}
