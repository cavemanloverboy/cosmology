use std::f64::consts::PI;

use quadrature::Output;
use rayon::prelude::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use self::{cdf::calculate_grf_cdf, pdf::calculate_grf_pdf};

// Very low to invoke max num of integrator iterations
const GRF_INTEGRATION_TARGET_ERROR: f64 = 1e-20;

/// These are in separate modules because they are massive and are
/// stressing VS code being open on a 48 core, 256 GB RAM machine...
mod cdf;
/// These are in separate modules because they are massive and are
/// stressing VS code being open on a 48 core, 256 GB RAM machine...
mod pdf;

/// A helper struct for calculating the k-Nearest Neighbor distributions
/// of gaussian random fields. The intended use is as follows:
///
/// 1) Initialize a [GaussianRandomField] struct with some number density
/// `nbar` and a [SpaceMode]
/// 2) Push distances to the struct (as you would a vector). Internally,
/// this will calcuate I (or its components) and dIdV (or its components).
/// 3) Use the struct's methods to calculate the values of the CDF or PDF at
/// the specified scales.
///
/// TODO: sample usage
pub struct GaussianRandomField<'b, 'c> {
    /// Specified [SpaceMode]. Either [SpaceMode::RealSpace] or [SpaceMode::RedshiftSpace]
    mode: SpaceMode<'b, 'c>,

    /// Container for distances (either measurements or scales of interest),
    /// the integral I, and its derivative (or its components if [SpaceMode::RedshiftSpace])
    container: Vec<(f64, TwoPointIntegral, TwoPointIntegralDerivative)>,
}

type TwoPointIntegral = f64;
type TwoPointIntegralDerivative = TwoPointIntegral;

type Radius = f64;

/// Specify what space to calculate the kNN CDFs for,
/// either [SpaceMode::RealSpace] or [SpaceMode::RedshiftSpace].
pub enum SpaceMode<'borrow, 'corr> {
    RealSpace(&'borrow (dyn Fn(Radius) -> f64 + 'corr + Send + Sync)),
    RedshiftSpace {
        real_corr: &'borrow (dyn Fn(Radius) -> f64 + 'corr + Send + Sync),
        f: f64,
    },
}

impl<'b, 'c> GaussianRandomField<'b, 'c> {
    /// Initializes the [GaussianRandomField] struct using `Vec::new` for the underlying vectors.
    pub fn new(mode: SpaceMode<'b, 'c>) -> Self {
        GaussianRandomField {
            mode,
            container: Vec::new(),
        }
    }

    // /// Initializes the [GaussianRandomField] struct and passes in the specified capacity
    // /// to the underlying vectors for better, preemtive allocations.
    // pub fn new_with_capacity(
    //     nbar: f64,
    //     mode: SpaceMode<'b, 'c>,
    //     capacity: usize,
    // ) -> Self {
    //     GaussianRandomField {
    //         nbar,
    //         mode,
    //         distances: Vec::with_capacity(capacity)
    //     }
    // }

    // /// Initializes the [GaussianRandomField] struct and passes in the specified capacity
    // /// to the underlying vectors for better, preemtive allocations.
    // #[cfg(feature = "use-mpi")]
    // pub fn with(mut self, rs: Option<Vec<f64>>, universe: Arc<Universe>) -> Self {
    //     // Initialize balancer
    //     let balancer: Balancer<Triplet> = Balancer::new(universe, false);
    //     let work = |r: &f64| Triplet((*r, self.calculate_i(*r), self.calculate_didv(*r)));
    //     if balancer.rank == 0 {
    //         let ours = balancer.distribute(Some(rs.unwrap())).unwrap();
    //         balancer.work(&ours, work);
    //         let output = balancer.collect().unwrap();
    //         // SAFETY: The triplet is a #[repr(transparent)] on the target type
    //         self.container =
    //             unsafe { std::mem::transmute::<Vec<Triplet>, Vec<(f64, f64, f64)>>(output) };
    //     } else {
    //         // Receive our work
    //         let ours: Vec<f64> = balancer.distribute(None).unwrap();
    //         // Do our work
    //         balancer.work(&ours, work);
    //         // Send back
    //         balancer.collect();
    //     }

    //     self
    // }

    /// Initializes the [GaussianRandomField] struct and passes in the specified capacity
    /// to the underlying vectors for better, preemtive allocations.
    // #[cfg(not(feature = "use-mpi"))]
    pub fn with(mut self, rs: &[f64]) -> Self {
        use rayon::prelude::IntoParallelIterator;

        self.container = rs
            .into_par_iter()
            .map(|r| (*r, self.calculate_i(*r), self.calculate_didv(*r)))
            .collect();

        self
    }

    // /// Adds a distance to the object and calculates its corresponding I and dIdV (or its components)
    // pub fn push(
    //     &mut self,
    //     r: f64,
    // ) {

    //     // Add distance to internal distances
    //     self.distances.push(r);

    //     // Calculate I (or components) and add it to the struct
    //     let i = self.calculate_i(r);
    //     self.i_container.push(i);

    //     // Calculate I (or components) and add it to the struct
    //     let didv = self.calculate_didv(r);
    //     self.didv_container.push(didv);
    // }

    /// Internal function. Calculate i for a given distance scale
    fn calculate_i(&self, r: f64) -> TwoPointIntegral {
        match self.mode {
            SpaceMode::RealSpace(corr) => {
                // Prefactor of integral (polar and azimuthal integrals)
                let prefactor = 4.0 * PI;

                // Integrand of i
                let integrand = |s: f64| {
                    // W(s) * Xi(s) * s^2
                    spherical_pair_weight_function(s, r) * corr(s) * s.powi(2)
                };

                // Carry out integral -- only need to integrate to twice the separation,
                // since the weight function is zero for s > 2r.
                let integral = quadrature::clenshaw_curtis::integrate(
                    integrand,
                    0.0,
                    2.0 * r,
                    GRF_INTEGRATION_TARGET_ERROR,
                );
                check_integral(&integral);

                let result: f64 = prefactor * integral.integral;
                result
            }

            // NOTE: this calculate I_r, the real-space integral.
            // The redshift-space integral is constructed later.
            SpaceMode::RedshiftSpace { real_corr, .. } => {
                // Prefactor of integral (polar and azimuthal integrals)
                let prefactor = 4.0 * PI;

                // Integrand of i
                let integrand = |s: f64| {
                    // W(s) * Xi(s) * s^2
                    spherical_pair_weight_function(s, r) * real_corr(s) * s.powi(2)
                };

                // Carry out integral -- only need to integrate to twice the separation,
                // since the weight function is zero for s > 2r.
                let integral = quadrature::clenshaw_curtis::integrate(
                    integrand,
                    0.0,
                    2.0 * r,
                    GRF_INTEGRATION_TARGET_ERROR,
                );
                check_integral(&integral);

                let result: f64 = prefactor * integral.integral;
                result
            }
        }
    }

    /// Internal function. Calculate didv for a distance scale. The methodology is exactly the same
    /// as for i, but in the integrand we replace W with dW/dV
    fn calculate_didv(&self, r: f64) -> TwoPointIntegralDerivative {
        match self.mode {
            SpaceMode::RealSpace(corr) => {
                // Prefactor of integral (polar and azimuthal integrals)
                let prefactor = 4.0 * PI;

                // Integrand of i
                let integrand = |s: f64| {
                    // dW/dV(s) * Xi(s) * s^2
                    spherical_pair_weight_function_dv(s, r) * corr(s) * s.powi(2)
                };

                // Carry out integral -- only need to integrate to twice the separation,
                // since the weight function is zero for s > 2r.
                let integral = quadrature::clenshaw_curtis::integrate(
                    integrand,
                    0.0,
                    2.0 * r,
                    GRF_INTEGRATION_TARGET_ERROR,
                );
                check_integral(&integral);

                let result: f64 = prefactor * integral.integral;
                result
            }

            // NOTE: this calculate dI_r/dr, the real-space integral derivative.
            // The redshift-space integral derivative is constructed later.
            SpaceMode::RedshiftSpace { real_corr, .. } => {
                // Prefactor of integral (polar and azimuthal integrals)
                let prefactor = 4.0 * PI;

                // Integrand of i
                let integrand = |s: f64| {
                    // dW/dV(s) * Xi(s) * s^2
                    spherical_pair_weight_function_dv(s, r) * real_corr(s) * s.powi(2)
                };

                // Carry out integral -- only need to integrate to twice the separation,
                // since the weight function is zero for s > 2r.
                let integral = quadrature::clenshaw_curtis::integrate(
                    integrand,
                    0.0,
                    2.0 * r,
                    GRF_INTEGRATION_TARGET_ERROR,
                );
                check_integral(&integral);

                let result: f64 = prefactor * integral.integral;
                result
            }
        }
    }

    pub fn get_cdf(&self, k: u8, nbar: f64, b: Option<f64>) -> Vec<f64> {
        self.container
            .iter()
            .enumerate()
            .map(|(_, (r, i, _))| {
                calculate_grf_cdf(*r, nbar, self.map_integral(*i, b.unwrap_or(1.0)), k)
            })
            .collect()
    }

    // #[cfg(not(feature = "use-mpi"))]
    pub fn get_pdf(&self, k: u8, nbar: f64, b: Option<f64>) -> Vec<f64> {
        self.container
            .par_iter()
            .map(|(r, i, didv)| {
                calculate_grf_pdf(
                    *r,
                    nbar,
                    self.map_integral(*i, b.unwrap_or(1.0)),
                    self.map_integral(*didv, b.unwrap_or(1.0)),
                    k,
                )
            })
            .collect()
    }

    // #[cfg(feature = "use-mpi")]
    // pub fn get_pdf(&self, k: u8, nbar: f64, b: Option<f64>, universe: Arc<Universe>) -> Vec<f64> {
    //     // Initialize balancer

    //     let balancer = Balancer::new(universe, false);
    //     let b = b.unwrap_or(1.0);

    //     let work = |(r, i, didv): &(f64, f64, f64)| {
    //         calculate_grf_pdf(
    //             *r,
    //             nbar,
    //             self.map_integral(*i, b),
    //             self.map_integral(*didv, b),
    //             k,
    //         )
    //     };
    //     if balancer.rank == 0 {
    //         // SAFETY: The triplet is a #[repr(transparent)] on the target type
    //         // NOTE: unnecessary clone but lets go with it...
    //         let distribute: Vec<Triplet> = unsafe {
    //             std::mem::transmute::<Vec<(f64, f64, f64)>, Vec<Triplet>>(self.container.clone())
    //         };
    //         let ours: Vec<(f64, f64, f64)> = unsafe {
    //             std::mem::transmute::<Vec<Triplet>, Vec<(f64, f64, f64)>>(
    //                 balancer.distribute(Some(distribute)).unwrap(),
    //             )
    //         };
    //         // Do our work
    //         balancer.work(&ours, work);
    //         // Collect all work
    //         let output = balancer.collect().unwrap();
    //         return output;
    //     } else {
    //         // Receive our work
    //         let ours: Vec<(f64, f64, f64)> = unsafe {
    //             std::mem::transmute::<Vec<Triplet>, Vec<(f64, f64, f64)>>(
    //                 balancer.distribute(None).unwrap(),
    //             )
    //         };
    //         // Do our work
    //         balancer.work(&ours, work);
    //         // Send back
    //         balancer.collect();
    //         return vec![];
    //     }
    // }

    fn map_integral(&self, i: f64, b: f64) -> f64 {
        match &self.mode {
            SpaceMode::RealSpace(_) => b * b * i,
            SpaceMode::RedshiftSpace { f, .. } => {
                i * (b * b + (2.0 * b * f / 3.0) + (f.powi(2) / 5.0))
            }
        }
    }

    pub fn get_i_didv<'a>(&'a self) -> &[(f64, TwoPointIntegral, TwoPointIntegralDerivative)] {
        &self.container
    }
}

#[allow(unused)]
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
    ((PI * (-2.0 * r + s).powi(2)) / 3. - (PI * (-2.0 * r + s) * (4.0 * r + s)) / 3.)
        / (4. * PI * r.powi(2))
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
