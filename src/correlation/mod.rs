use std::f64::consts::PI;
use std::error::Error;
use quadrature::Output;
use ouroboros::self_referencing;
use crate::{power::{PowerSpectrum, PowerFn, TransferFunction}, utils::math::spherical_jn};


/// Default parameter for the lower k-bound of the correlation function integral
pub const CORR_LOGK_MIN: f64 = -10.0;
/// Default parameter for the upper k-bound of the correlation function integral
pub const CORR_LOGK_MAX: f64 = 10.0;
// /// Default parameter for the number of points used for the correlation function integral
pub const CORR_ABS_ERROR: f64 = 1e-6;


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
    pub fn correlation_function(
        &self,
        r: f64
    ) -> f64 {

        // Define prefactor and integrand
        let prefactor = (2.0 * PI).powi(-2);
        let integrand = |k: f64| {
            k.powi(2) * self.borrow_power_at_k().power(k) * (k * r).sin() / (k * r)
        };

        // Carry out integral
        let cf = quadrature::integrate(
            integrand,
            10_f64.powf(*self.borrow_lower_logk()),
            10_f64.powf(*self.borrow_upper_logk()),
            *self.borrow_target_error()
        );
        check_integral(&cf);

        // Return result
        prefactor * cf.integral
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
        params: CorrelationFunctionParameters
    ) -> Result<CorrelationFunction, Box<dyn Error>> {

        // Specify domain of integration, target error
        let (lower_logk, upper_logk, target_error) = {
            if let Some(ref acc_params) = params.accuracy_params {
                // User specified bounds, error
                (acc_params.lower_logk_bound, acc_params.upper_logk_bound, acc_params.target_error)
            } else {
                // Default parameters if None
                (CORR_LOGK_MIN, CORR_LOGK_MAX, CORR_ABS_ERROR)
            }
        };
  
        Ok(
            CorrelationFunctionBuilder {
                z,
                params,
                lower_logk,
                upper_logk,
                target_error,
                power_at_k_builder: |params: &CorrelationFunctionParameters| params.power.power_fn(z).unwrap()
            }.build()
        )
    }
}

#[allow(unused)]
fn check_integral(cf: &Output) -> Result<(), Box<dyn Error>> {
    // TODO: check integral for convergence, perhaps just print warnings
    Ok(())
}


/// The parameters required for calculating the 2-point correlation function.
#[derive(Clone)]
pub struct CorrelationFunctionParameters {
    
    /// The underlying power spectrum engine
    pub power: PowerSpectrum,

    /// Parameters controlling the accuracy of the correlation function
    pub accuracy_params: Option<CorrFuncAccuracyParameters>
    
}

// Parameters controlling the accuracy of the calculated correlation function.
#[derive(Clone)]
pub struct CorrFuncAccuracyParameters {
    lower_logk_bound: f64,
    upper_logk_bound: f64,
    target_error: f64,
}