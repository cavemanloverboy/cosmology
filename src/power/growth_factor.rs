use crate::scale_factor::rk4_multi;



/// Very low to invoke max num of integrator iterations
const TARGET_ABSOLUTE_ERROR: f64 = 1e-30;
const DEFAULT_DELTA_SCALE_FACTOR: f64 = 1e-4;

/// Currently implemented only for flat universes with w = -1.
pub(crate) fn linear_growth_factor(
    omega_0: f64,
    omega_de: f64,
    z: f64
) -> Result<f64, String> {

    // Only flat universe are supported
    validate(omega_0, omega_de)?;
    let omega_r = 1.0 - omega_0 - omega_de;

    // Only universes with matter have nonzero linear growth factors
    if contains(omega_0) {
        

        if !contains(omega_de) && !contains(omega_r) {

            // If we are in a matter-only universe, then the growth mode is simply a(z)
            Ok((1.0 + z).recip())

        } else if contains(omega_de) && !contains(omega_r) {

            // If we are in a matter + de universe, then we calculate
            // H(a) * 5 Om / 2 int_0^a da / (a^3 H(a))
            let result = {

                // Friedmann equation with nonzero components
                // Prefactor of the square root, H0, cancels, so exclude it.
                let ha = |a: f64| (omega_0 * a.powi(-3) + omega_de).sqrt();

                // Construct integral parameters
                let integrand = |a: f64| a.powi(-3) * ha(a).powi(-3);
                // let integrand = |a: f64, _: f64| a.powi(-3) * ha(a).powi(-3);
                let lower_a = 1e-6;
                let upper_a = (1.0 + z).recip();

                // Prefactor
                let prefactor = ha(upper_a) * 5.0 * omega_0 / 2.0;

                // Integrate
                let result = quadrature::clenshaw_curtis::integrate(
                    integrand,
                    lower_a,
                    upper_a,
                    TARGET_ABSOLUTE_ERROR
                );
                // let result = {
                //     let mut current_a = lower_a;
                //     let mut area = lower_a;
                //     while current_a < upper_a {
                //         let da = DEFAULT_DELTA_SCALE_FACTOR.min(upper_a-current_a);
                //         area = rk4(integrand, current_a, area, da, None);
                //         current_a += da;
                //     }
                //     area
                // };

                // if result.error_estimate / result.integral > 10.0 * TARGET_ABSOLUTE_ERROR {
                //     // return Err(
                //     println!("{}",
                //         format!(
                //             "integral for linear growth factor D(a) did not converge. Error estimate: {}. Number of evaluations: {}",
                //             result.error_estimate / result.integral,
                //             result.num_function_evaluations
                //         )
                //     )
                //     // )
                // }

                // Multiply by prefactor
                // prefactor * result
                prefactor * result.integral
            };

            Ok(result)

        } else if !contains(omega_de) && contains(omega_r) {

            // If we are in a matter + radiation universe, analytic answer exists
            // Taken from Bernardeau et al 2001, which references Peebles & Groth 1975
            let x = omega_0.recip() - 1.0;
            let result = 1.0 
                + 3.0 / x 
                + 3.0 * ((1.0 + x) / x.powi(3)).sqrt() * ((1.0 + x).sqrt() - x.sqrt()).ln();

            Ok(result)

        } else if contains(omega_de) && contains(omega_r) {

            // If we are in a universe with matter, radiation, and de,
            // need explicit integration of the second order ode.
            //
            // Inspired by Linder & Jenkins 2003, we solve the coupled first-order ODEs.
            // That is, we define F = G' then we end up with the coupled first-order ODEs
            // G' = F
            // F' = -[7/2 - 3/2 w(a) / (1 + X(a))] F / a - (1-w(a))/(1+X(a)) G / a^2
            let w = |_a: f64| -1.0;
            let x = |a: f64| omega_0 * a.powi(-3) / (1.0 - omega_0 * a.powi(3));
            let g_prime = |_a: f64, g_and_f: [f64; 2]| g_and_f[1];
            let f_prime = |a: f64, g_and_f: [f64; 2]| {
                let [g, f] = g_and_f;
                -(7.0/2.0 - 3.0/2.0 * w(a) / (1.0 + x(a))) * f / a - (1.0 - w(a))/(1.0 + x(a)) * g / a.powi(2)
            };

            // Package derivative functions
            let derivatives = |a: f64, g_and_f: [f64; 2]| [g_prime(a, g_and_f), f_prime(a, g_and_f)];

            // In the radiation domination era, D ~ 0.
            // In the early matter domination era, D ~ a. Thus, lim a->0 G(a) = lim a->0 1 -> 1. 
            // Consequently, 1' = 0.
            let mut current_g_and_f = [1.0, 0.0];
            let mut current_a = 1e-5;
            let target_a = (1.0 + z).recip();
            while current_a < target_a {

                // Get step size
                let da = DEFAULT_DELTA_SCALE_FACTOR.min(target_a-current_a);

                // Update functions and independent variable
                current_g_and_f = rk4_multi(derivatives, current_a, current_g_and_f, da);
                current_a += da;
            }

            // Unpack G and transform, discard F
            let [final_g, _] = current_g_and_f;
            let growth_factor = final_g * (1.0 + z).recip();

            Ok(growth_factor)

        } else { unreachable!("all cases explicitly outlined for clarity") }
    } else {

        // No matter --> no matter overdensity
        Ok(0.0)
    }


}




fn contains(a: f64) -> bool {
    a.is_sign_positive() && a.is_normal()
}

fn validate(omega_0: f64, omega_de: f64) -> Result<(), String> {
    let is_valid = {
        omega_0 >= 0.0
        && omega_de >= 0.0
        && (omega_0 + omega_de) <= 1.0
    };

    if is_valid {
        Ok(())
    } else {
        Err(format!("The specified cosmological parameters are not valid. Only flat universe are supported"))
    }
}