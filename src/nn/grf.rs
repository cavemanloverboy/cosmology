use std::f64::consts::PI;



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

    /// Average number density
    nbar: f64,

    /// Specified [SpaceMode]. Either [SpaceMode::RealSpace] or [SpaceMode::RedshiftSpace]
    mode: SpaceMode<'b, 'c>,

    /// Container for distances -- either measurements or scales of interest.
    distances: Vec<f64>,

    /// Container for I if in [SpaceMode::RealSpace], or of I if in
    /// [SpaceMode::RedshiftSpace]
    i_container: Vec<TwoPointIntegral>,

    /// Container for dIdV if in [SpaceMode::RealSpace], or of dIdV if in
    /// [SpaceMode::RedshiftSpace]. Reuses the [TwoPointIntegral] type.
    didv_container: Vec<TwoPointIntegralDerivative>,
}

pub enum TwoPointIntegral {
    RealSpace(f64),
    RedshiftSpace {
        i0: f64,
        i1: f64,
        i2: f64,
    }
}
type TwoPointIntegralDerivative = TwoPointIntegral;


type Radius = f64;
type Angle = f64;
/// Specify what space to calculate the kNN CDFs for, 
/// either [SpaceMode::RealSpace] or [SpaceMode::RedshiftSpace].
pub enum SpaceMode<'borrow, 'corr> {
    RealSpace(&'borrow (dyn Fn(Radius) -> f64 + 'corr)),
    RedshiftSpace(&'borrow (dyn Fn(Radius, Angle) -> f64 + 'corr)),
}

impl<'b, 'c> GaussianRandomField<'b,'c> {

    /// Initializes the [GaussianRandomField] struct using `Vec::new` for the underlying vectors.
    pub fn new(
        nbar: f64,
        mode: SpaceMode<'b, 'c>,
    ) -> Self {
        GaussianRandomField {
            nbar,
            mode,
            distances: Vec::new(),
            i_container: Vec::new(),
            didv_container: Vec::new()
        }
    }

    /// Initializes the [GaussianRandomField] struct and passes in the specified capacity 
    /// to the underlying vectors for better, preemtive allocations.
    pub fn new_with_capacity(
        nbar: f64,
        mode: SpaceMode<'b, 'c>,
        capacity: usize,
    ) -> Self {
        GaussianRandomField {
            nbar,
            mode,
            distances: Vec::with_capacity(capacity),
            i_container: Vec::with_capacity(capacity),
            didv_container: Vec::with_capacity(capacity)
        }
    }

    /// Adds a distance to the object and calculates its corresponding I and dIdV (or its components)
    pub fn push(
        &mut self,
        r: f64,
    ) {
        
        // Add distance to internal distances
        self.distances.push(r);

        // Calculate I (or components) and add it to the struct
        let i = self.calculate_i(r);
        self.i_container.push(i);

        // Calculate I (or components) and add it to the struct
        let didv = self.calculate_didv(r);
        self.didv_container.push(didv);
    }

    /// Internal function. Calculate i for a distance scale
    fn calculate_i(
        &self,
        r: f64,
    ) -> TwoPointIntegral {
        match self.mode {
            SpaceMode::RealSpace(corr) => {

                let result: f64 = todo!();
                TwoPointIntegral::RealSpace(result)
            },
            SpaceMode::RedshiftSpace(corr_2d) => {

                let (i0, i1, i2): (f64, f64, f64) = todo!();
                TwoPointIntegral::RedshiftSpace { i0, i1, i2 }
            }
        }
    }

    /// Internal function. Calculate didv for a distance scale
    fn calculate_didv(
        &self,
        r: f64,
    ) -> TwoPointIntegralDerivative {
        match self.mode {
            SpaceMode::RealSpace(corr) => {

                let result: f64 = todo!();
                TwoPointIntegral::RealSpace(result)
            },
            SpaceMode::RedshiftSpace(corr_2d) => {

                let (i0, i1, i2): (f64, f64, f64) = todo!();
                TwoPointIntegral::RedshiftSpace { i0, i1, i2 }
            }
        }
    }
}


/// This weight function is a probe for the number of 2-point configurations
/// with separation `s` in a sphere of radius `r`.
fn spherical_pair_weight_function(
    s: f64,
    r: f64,
) -> f64 {
    PI / 12.0 * (s - 2.0 * r).powi(2) * (s + 4.0 * r)
}



#[test]
fn test_spherical_weight_at_diameter() {

    // The volume formed by all two-point configurations where the two points
    // are separated by a sphere's diameter in a sphere is zero.
    assert_approx_eq::assert_approx_eq!(
        spherical_pair_weight_function(2.0, 1.0),
        0.0,
        1e-5
    );
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




