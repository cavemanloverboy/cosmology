use std::ops::Sub;
use std::error::Error;
use quadrature::Output;
use ouroboros::self_referencing;
use crate::power::{PowerSpectrum, PowerFn};



/// Default parameter for the lower k-bound of the correlation function integral
pub const CORR_LOGK_MIN: f64 = -10.0;
/// Default parameter for the upper k-bound of the correlation function integral
pub const CORR_LOGK_MAX: f64 = 10.0;
// /// Default parameter for the number of points used for the correlation function integral
// pub const CORR_K_POINTS: usize = 200;
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

    pub fn correlation_function(
        &self,
        r: f64
    ) -> f64 {
        let integrand = |k: f64| {
            k.powi(2) * self.borrow_power_at_k().power(k) * (k * r).sin() / (k * r)
        };
        let cf = quadrature::integrate(
            integrand,
            10_f64.powf(*self.borrow_lower_logk()),
            10_f64.powf(*self.borrow_upper_logk()),
            *self.borrow_target_error()
        );

        check_integral(&cf);
        cf.integral
    }

    pub fn get_correlation_function(
        z: f64,
        params: CorrelationFunctionParameters
    ) -> Result<CorrelationFunction, Box<dyn Error>> {

        // Get the power spectrum over the specified domain
        // To do that, first get the wavenumbers
        // let ks: Vec<f64> = {
        let lower_logk;
        let upper_logk;
        // let num_k_points;
        let target_error;
        if let Some(ref acc_params) = params.accuracy_params {
            lower_logk = acc_params.lower_logk_bound;
            upper_logk = acc_params.upper_logk_bound;
            // num_k_points = acc_params.num_k_points;
            target_error = acc_params.target_error;
        } else {
            lower_logk = CORR_LOGK_MIN;
            upper_logk = CORR_LOGK_MAX;
            // num_k_points = CORR_K_POINTS;
            target_error = CORR_ABS_ERROR;
        }
        
            // The subtraction operation here because the first point is zero steps
            // away from the lower bound.
            // let dlogk_step = (upper_logk-lower_logk)/(num_k_points.sub(1) as f64);

            // (0..CORR_K_POINTS)
            //     .map(|i| 10.0_f64.powf(lower_logk + i as f64 * dlogk_step))
            //     .collect()
        // };

        // // Then get the power over the specified domain
        // let power: Vec<f64> = params.power.calculate_power(&ks, z)?;
  
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
    // num_k_points: usize,
}