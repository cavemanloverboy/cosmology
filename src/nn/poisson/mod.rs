use quadrature::Output;
use cdf::calculate_poisson_cdf;

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
        self.container = rs.to_vec();

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

#[allow(unused)]
fn check_integral(integral: &Output) {
    // TODO
}