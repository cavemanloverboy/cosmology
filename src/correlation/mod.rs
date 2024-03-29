use crate::power::{PowerFn, PowerSpectrum};
use ouroboros::self_referencing;
use quadrature::Output;
use std::error::Error;
use std::f64::consts::PI;
use std::ops::Add;

/// Default parameter for the lower k-bound of the correlation function integral
pub const CORR_LOGK_MIN: f64 = -6.0;
/// Default parameter for the upper k-bound of the correlation function integral
pub const CORR_LOGK_MAX: f64 = 6.0;
/// Default parameter for the abs error of the integrator
pub const CORR_ABS_ERROR: f64 = 1e-10;
// /// At r = 100, this gives 1/100th of a period
// const MAX_K_STEP: f64 = 2.0 * PI / 100.0 / 100.0;

#[self_referencing]
pub struct CorrelationFunction {
    z: f64,
    params: CorrelationFunctionParameters,
    lower_logk: f64,
    upper_logk: f64,
    target_error: f64,
    #[borrows(params)]
    #[covariant]
    power_at_k: PowerFn<'this>,
}

impl CorrelationFunction {
    /// Calculates the linear theory correlation function in real-space
    pub fn correlation_function(&self, r: f64) -> f64 {
        // Define prefactor
        let prefactor = (2.0 * (PI).powi(2)).recip();

        // Linear integrand
        let integrand = |k: f64| k * self.borrow_power_at_k().power(k) * (k * r).sin() / r;

        // Carry out integral
        let intervals = r.clamp(1.0, 4.0) as usize;
        let mut result = 0.0;
        let lower_k = 1e-6 / r;
        let upper_k = 1e4 / r;
        let total_interval = upper_k - lower_k;
        let subinterval_size = total_interval / intervals as f64;
        for i in 0..intervals {
            let interval_lower = lower_k + subinterval_size * i as f64;
            let interval_upper = lower_k + subinterval_size * i.add(1) as f64;
            let cf = quadrature::double_exponential::integrate(
                integrand,
                interval_lower,
                interval_upper,
                *self.borrow_target_error(),
            );
            check_integral(&cf);
            result += cf.integral;
        }
        prefactor * result
    }

    // /// Returns the linear theory correlation function RSD monopole.
    // /// Inspired by Kaiser 1987, Hamilton 1992 RSDs (in redshift-space).
    // pub fn rsd_correlation_function_monopole(
    //     &self,
    //     r: f64,
    //     b: f64,
    // ) -> (f64, f64, f64) {

    //     // Get linear theory real-space monopole
    //     let corr_real = self.correlation_function(r);

    //     // Calculate growth rate
    //     let omega_matter = self.borrow_params().power.get_omega_matter();
    //     let f = omega_matter.powf(5.0/9.0);

    //     // Return result
    //     (corr_real, prefactor * cf_2.integral, prefactor * cf_4.integral)
    // }

    pub fn get_correlation_function(
        z: f64,
        params: CorrelationFunctionParameters,
    ) -> Result<CorrelationFunction, Box<dyn Error>> {
        // Specify domain of integration, target error
        let (lower_logk, upper_logk, target_error) = {
            if let Some(ref acc_params) = params.accuracy_params {
                // User specified bounds, error
                (
                    acc_params.lower_logk_bound,
                    acc_params.upper_logk_bound,
                    acc_params.target_error,
                )
            } else {
                // Default parameters if None
                (CORR_LOGK_MIN, CORR_LOGK_MAX, CORR_ABS_ERROR)
            }
        };

        Ok(CorrelationFunctionBuilder {
            z,
            params,
            lower_logk,
            upper_logk,
            target_error,
            power_at_k_builder: |params: &CorrelationFunctionParameters| {
                params.power.power_fn(z).unwrap()
            },
        }
        .build())
    }

    pub fn get_params<'a>(&'a self) -> &'a CorrelationFunctionParameters {
        &self.borrow_params()
    }
    pub fn get_redshift<'a>(&'a self) -> f64 {
        *self.borrow_z()
    }
}

#[allow(unused)]
fn check_integral(cf: &Output) {
    // TODO: check integral for convergence, perhaps just print warnings
    // println!("evaluation: {}", cf.num_function_evaluations);
}

/// The parameters required for calculating the 2-point correlation function.
#[derive(Clone)]
pub struct CorrelationFunctionParameters {
    /// The underlying power spectrum engine
    pub power: PowerSpectrum,

    /// Parameters controlling the accuracy of the correlation function
    pub accuracy_params: Option<CorrFuncAccuracyParameters>,
}

// Parameters controlling the accuracy of the calculated correlation function.
#[derive(Clone)]
pub struct CorrFuncAccuracyParameters {
    lower_logk_bound: f64,
    upper_logk_bound: f64,
    target_error: f64,
}

#[cfg(test)]
#[cfg(feature = "colossus-python")]
mod tests {

    use super::*;
    use crate::power::{PowerSpectrum, TransferFunction};

    macro_rules! assert_eq_tol {
        ($x:expr, $y:expr, $d:expr) => {
            // Calculate fractional delta
            let frac_delta = (($x - $y) / $y).abs();

            // Compare frac_delta
            let within = frac_delta < $d;

            if !within {
                // Construct err msg
                let msg = format!(
                    "Result {:.4e} not within {:.4e} of {:.4e}. frac_delta is {:.4e}",
                    $x, $d, $y, frac_delta,
                );

                // Panic with err msg
                panic!("{msg}");
            }
        };
    }

    macro_rules! eisenstein_corr(
        ($z:ident, $h0:ident, $om0:ident, $ob0:ident, $t0:ident) => {

          concat_idents::concat_idents!(test_name = test_eisen_corr, _, $z, _, $h0, _, $om0, _, $ob0, _, $t0, {
            #[test]
            fn test_name() {

                let z: u32 = stringify!($z)[1..].parse::<u32>().unwrap();
                let h: u32 = stringify!($h0)[1..].parse::<u32>().unwrap();
                let om0: u32 = stringify!($om0)[1..].parse::<u32>().unwrap();
                let ob0: u32 = stringify!($ob0)[1..].parse::<u32>().unwrap();
                let t0: u32 = stringify!($t0)[1..].parse::<u32>().unwrap();

                // Construct EisensteinHu model
                let power = PowerSpectrum::new(TransferFunction::EisensteinHu {
                    h: h as f64 / 100.0, // h
                    omega_matter_0: om0 as f64 / 100.0, // omega_matter_0
                    omega_baryon_0: ob0 as f64 / 100.0, // omega_baryon_0
                    temp_cmb0: t0 as f64 / 100.0, // temp_cmb_0
                    ns: 0.9665, // ns
                    sigma_8: 0.8102, // sigma8
                }).unwrap();

                // Construct correlation function
                let params = CorrelationFunctionParameters {
                    power,
                    accuracy_params: None, // default accuracy parameters
                };
                let corr = CorrelationFunction::get_correlation_function(
                    z as f64,
                    params
                ).unwrap();

              // Pick scales
              let rs = [0.1, 1.0, 10.0, 100.0];

              // Get result at redshift zero
              let result = rs.map(|r| corr.correlation_function(r));

              // Expected values, from COLOSSUS
              let expected = {
                use pyo3::prelude::*;
                use pyo3::types::*;
                Python::with_gil(|py| {

                  // Get ks into python
                  let list = PyList::new(py, &rs);
                  let locals = PyDict::new(py);
                  locals.set_item("rs", list).unwrap();

                  py.run(format!(r#"from colossus.cosmology import cosmology
import warnings
warnings.filterwarnings("ignore")
planck18 = cosmology.setCosmology("planck18")
params = {{
    "H0": {0},
    "Om0": {1},
    "Ob0": {2},
    "Tcmb0": {3},
    "ns": 0.9665,
    "sigma8": 0.8102,
}}
cosmology.addCosmology("test", params=params)
cosmo = cosmology.setCosmology("test")
x = []
for r in rs:
    x.append(cosmo.correlationFunction(r, z={4}))
                  "#, h as f64, om0 as f64 / 100.0, ob0 as f64 / 100.0,
                  t0 as f64 / 100.0, z).as_str(), None, Some(locals)).unwrap();
                  let x: Vec<_> = locals.get_item("x").unwrap().extract::<Vec<f64>>().unwrap();
                  x
                })
              };

              for i in 0..result.len() {
                let abs_diff = (result[i]-expected[i]).abs();
                let rel_diff = ((result[i]-expected[i])/expected[i]).abs();
                assert!(abs_diff < 1e-2 || rel_diff < 5e-2, "result {:.3e} has abs_diff = {abs_diff:.3e}; rel_diff = {rel_diff:.3e} wrt expectation {:.3e}", result[i], expected[i]);
              }
            }
          });
        }
      );
    dry::macro_for!($H in [h50, h60, h70, h80, h90, h100] {
        dry::macro_for!($M in [m10, m30, m50, m70, m90] {
            dry::macro_for!($B in [b1, b2, b3] {
                dry::macro_for!($T in [t270] {
                    dry::macro_for!($Z in [z0, z1, z2, z10] {
                        eisenstein_corr!($Z, $H, $M, $B, $T);
                    });
                });
            });
        });
    });

    macro_rules! eisenstein_corr_plot(
        ($z:ident, $h0:ident, $om0:ident, $ob0:ident, $t0:ident) => {

          concat_idents::concat_idents!(test_name = test_eisen_corr, _, $z, _, $h0, _, $om0, _, $ob0, _, $t0, {
            #[test]
            fn test_name() {

                let z: u32 = stringify!($z)[1..].parse::<u32>().unwrap();
                let h: u32 = stringify!($h0)[1..].parse::<u32>().unwrap();
                let om0: u32 = stringify!($om0)[1..].parse::<u32>().unwrap();
                let ob0: u32 = stringify!($ob0)[1..].parse::<u32>().unwrap();
                let t0: u32 = stringify!($t0)[1..].parse::<u32>().unwrap();

                // Construct EisensteinHu model
                let power = PowerSpectrum::new(TransferFunction::EisensteinHu {
                    h: h as f64 / 100.0, // h
                    omega_matter_0: om0 as f64 / 100.0, // omega_matter_0
                    omega_baryon_0: ob0 as f64 / 100.0, // omega_baryon_0
                    temp_cmb0: t0 as f64 / 100.0, // temp_cmb_0
                    ns: 0.9665, // ns
                    sigma_8: 0.8102, // sigma8
                }).unwrap();

                // Construct correlation function
                let params = CorrelationFunctionParameters {
                    power,
                    accuracy_params: None, // default accuracy parameters
                };
                let corr = CorrelationFunction::get_correlation_function(
                    z as f64,
                    params
                ).unwrap();

              // Pick scales
              let rs = [0.1, 1.0, 10.0, 100.0];

              // Get result at redshift zero
              let result = rs.map(|r| corr.correlation_function(r));

              // Expected values, from COLOSSUS
              let expected = {
                use pyo3::prelude::*;
                use pyo3::types::*;
                Python::with_gil(|py| {

                  // Get ks into python
                  let list = PyList::new(py, &rs);
                  let locals = PyDict::new(py);
                  locals.set_item("rs", list).unwrap();

                  py.run(format!(r#"from colossus.cosmology import cosmology
import warnings
warnings.filterwarnings("ignore")
planck18 = cosmology.setCosmology("planck18")
params = {{
    "H0": {0},
    "Om0": {1},
    "Ob0": {2},
    "Tcmb0": {3},
    "ns": 0.9665,
    "sigma8": 0.8102,
}}
cosmology.addCosmology("test", params=params)
cosmo = cosmology.setCosmology("test")
x = []
for r in rs:
    x.append(cosmo.correlationFunction(r, z={4}))
                  "#, h as f64, om0 as f64 / 100.0, ob0 as f64 / 100.0,
                  t0 as f64 / 100.0, z).as_str(), None, Some(locals)).unwrap();
                  let x: Vec<_> = locals.get_item("x").unwrap().extract::<Vec<f64>>().unwrap();
                  x
                })
              };

              for i in 0..result.len() {
                let abs_diff = (result[i]-expected[i]).abs();
                let rel_diff = ((result[i]-expected[i])/expected[i]).abs();
                assert!(abs_diff < 1e-2 || rel_diff < 5e-2, "result {:.3e} has abs_diff = {abs_diff:.3e}; rel_diff = {rel_diff:.3e} wrt expectation {:.3e}", result[i], expected[i]);
              }
            }
          });
        }
      );
}
