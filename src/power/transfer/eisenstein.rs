use quadrature::clenshaw_curtis::integrate;
use std::f64::consts::{E as EULER, PI};

use std::error::Error;
use crate::power::growth_factor::linear_growth_factor;





#[derive(Debug)]
pub struct EisensteinHu {
    /// Present Hubble (little h, in units of km/s/Mpc)
    h: f64,

    /// Present density of matter in units of critical density
    omega_matter_0: f64,

    /// Present density of baryons in units of critical density
    omega_baryon_0: f64,

    /// Present emperature of CMB
    temp_cmb0: f64,

    /// Spectral index
    ns: f64,

    /// Power at 8 Mpc/h
    sigma_8: f64,
}

impl EisensteinHu {
    /// Provide a set of cosmological parameters:
    ///
    /// `h: f64`: Hubble constant, little h
    /// `omega_matter_0: f64`: Present density of matter in units of critical density
    /// `omega_baryon_0: f64`: Present density of baryons in units of critical density
    /// `temp_cmb0: f64`: Present emperature of CMB
    /// `ns: f64`: Spectral index
    /// `sigma_8: f64`: Power at 8 Mpc/h
    /// 
    /// and get an `EisensteinHu` struct if inputs are valid. This struct has methods to
    /// calculate the power spectrum and transfer functions.
    /// 
    /// ```rust
    /// use cosmology::power::transfer::eisenstein::*;
    /// 
    /// let h = 0.7;
    /// let Om0 = 0.3;
    /// let Omb = 0.025;
    /// let Tcmb0 = 2.7;
    /// let ns = 0.96;
    /// let s8 = 0.90;
    /// let eisen_hu = EisensteinHu::new(
    ///     h,
    ///     Om0,
    ///     Omb,
    ///     Tcmb0,
    ///     ns,
    ///     s8
    /// ).expect("These are valid cosmo params");
    /// 
    /// let ks = [0.1, 1.0];
    /// let z = 0.0;
    /// let power = eisen_hu.power_z0(&ks);
    /// let power_nb = eisen_hu.power_z0_zero_baryon(&ks);
    /// 
    /// let transfer = ks
    ///     .map(|k| eisen_hu.transfer_baryon(k));
    /// let transfer_nb = ks
    ///     .map(|k| eisen_hu.transfer_zero_baryon(k));
    /// ```
    pub fn new(
        h: f64,
        omega_matter_0: f64,
        omega_baryon_0: f64,
        temp_cmb0: f64,
        ns: f64,
        sigma_8: f64,
    ) -> Result<EisensteinHu, &'static str> {
        if omega_baryon_0 < 0.0 {
            Err("The Eisenstein & Hu 98 transfer function cannot be computed for Ob0 < 0")
        } else if temp_cmb0 <= 0.0 {
            Err("Cannot have a nonpositive CMB temperature")
        } else if omega_baryon_0 > omega_matter_0 {
            Err("Cannot have more baryons than total matter")
        } else {
            Ok(EisensteinHu {
                h,
                omega_matter_0,
                omega_baryon_0,
                temp_cmb0,
                ns,
                sigma_8,
            })
        }
    }

    /// Power spectrum for wavenumber k at z = 0 (Eisenstein & Hu 1998).
    /// The code was adapted from Benedikt Diemer's COLOSSUS code,
    /// which was adapted from Matt Becker's cosmocalc code.
    ///
    /// Given a slice of wavenumbers `&[T]`, returns the transfer function
    /// with the same units as `k`.
    ///
    /// For example usage see [`EisensteinHu`].
    pub fn power_z0(&self, ks: &[f64]) -> Vec<f64> {
        ks.iter()
            .map(|&k| {
                let sigma_8 = self._sigma8_calc_baryon();
                let current_power = (self.sigma_8 / sigma_8).powi(2)
                    * self.transfer_baryon(k).powi(2)
                    * k.powf(self.ns);

                current_power
            })
            .collect::<Vec<_>>()
    }

    /// Power spectrum for wavenumber k at z = 0 (Eisenstein & Hu 1998).
    /// The code was adapted from Benedikt Diemer's COLOSSUS code,
    /// which was adapted from Matt Becker's cosmocalc code.
    ///
    /// Given a slice of wavenumbers `&[T]`, returns the transfer function
    /// with the same units as `k`.
    ///
    /// For example usage see [`EisensteinHu`].
    pub fn power_z(&self, ks: &[f64], z: f64) -> Result<Vec<f64>, Box<dyn Error>> {
        let power_at_z0 = self.power_z0(ks);
        let growth_factor = linear_growth_factor(
            self.h,
            self.omega_matter_0,
            1.0 - self.omega_matter_0,
            z
        )?;
        let growth_factor_norm = linear_growth_factor(
            self.h,
            self.omega_matter_0,
            1.0 - self.omega_matter_0,
            0.0 // z = 0.0
        )?;
        Ok(power_at_z0
            .into_iter()
            .map(|p0| p0 * (growth_factor / growth_factor_norm).powi(2))
            .collect())
    }

    /// Power spectrum for wavenumber k at z = 0 (Eisenstein & Hu 1998).
    /// The code was adapted from Benedikt Diemer's COLOSSUS code,
    /// which was adapted from Matt Becker's cosmocalc code. This version
    /// uses the zero baryon model, which omits BAOs.
    ///
    /// Given a slice of wavenumbers `&[T]`, returns the transfer function
    /// with the same units as `k`.
    /// 
    /// API is identical to that of the `power` method. For example usage,
    /// see [`EisensteinHu`].
    pub fn power_z0_zero_baryon(&self, ks: &[f64]) -> Vec<f64> {
        ks.iter()
            .map(|&k| {
                let sigma_8 = self._sigma8_calc_zero_baryon();
                let current_power = (self.sigma_8 / sigma_8).powi(2)
                    * self.transfer_zero_baryon(k).powi(2)
                    * k.powf(self.ns);

                current_power
            })
            .collect::<Vec<_>>()
    }

    pub fn power_z_zero_baryon(&self, ks: &[f64], z: f64) -> Result<Vec<f64>, Box<dyn Error>> {
        let power_at_z0 = ks.iter()
            .map(|&k| {
                let sigma_8 = self._sigma8_calc_zero_baryon();
                let current_power = (self.sigma_8 / sigma_8).powi(2)
                    * self.transfer_zero_baryon(k).powi(2)
                    * k.powf(self.ns);

                current_power
            })
            .collect::<Vec<_>>();
        let growth_factor = linear_growth_factor(
            self.h,
            self.omega_matter_0,
            1.0 - self.omega_matter_0,
            z
        )?;
        let growth_factor_norm = linear_growth_factor(
            self.h,
            self.omega_matter_0,
            1.0 - self.omega_matter_0,
            0.0 // z = 0.0
        )?;
        Ok(power_at_z0
            .into_iter()
            .map(|p0| p0 * (growth_factor / growth_factor_norm).powi(2))
            .collect())
    }

    /// Transfer function, with baryonic effects. 
    /// For example usage, see [`EisensteinHu`].
    pub fn transfer_baryon(&self, k: f64) -> f64 {
        // Auxilary values
        let omc = self.omega_matter_0 - self.omega_baryon_0;
        let ombom0 = self.omega_baryon_0 / self.omega_matter_0;
        let omcom0 = omc / self.omega_matter_0;
        let h2 = self.h.powi(2);
        let om0h2 = self.omega_matter_0 * h2;
        let ombh2 = self.omega_baryon_0 * h2;
        let theta2p7 = self.temp_cmb0 / 2.7;
        let theta2p72 = theta2p7.powi(2);
        let theta2p74 = theta2p72.powi(2);

        // Convert kh from h/Mpc to 1/Mpc
        let kh = k * self.h;

        // Equation 2
        let zeq = 2.50e4 * om0h2 / theta2p74;

        // Equation 3
        let keq = 7.46e-2 * om0h2 / theta2p72;

        // Equation 4
        let b1d = 0.313 * om0h2.powf(-0.419) * (1.0 + 0.607 * om0h2.powf(0.674));
        let b2d = 0.238 * om0h2.powf(0.223);
        let zd = 1291.0 * om0h2.powf(0.251) / (1.0 + 0.659 * om0h2.powf(0.828))
            * (1.0 + b1d * ombh2.powf(b2d));

        // Equation 5
        let rd = 31.5 * ombh2 / theta2p74 / (zd / 1e3);
        let req = 31.5 * ombh2 / theta2p74 / (zeq / 1e3);

        // Equation 6
        let s = 2.0 / 3.0 / keq
            * sqrt(6.0 / req)
            * ln((sqrt(1.0 + rd) + sqrt(rd + req)) / (1.0 + sqrt(req)));

        // Equation 7
        let ksilk = 1.6 * ombh2.powf(0.52) * om0h2.powf(0.73) * (1.0 + (10.4 * om0h2).powf(-0.95));

        // Equation 10
        let q = kh / 13.41 / keq;

        // Equation 11
        let a1 = (46.9 * om0h2).powf(0.670) * (1.0 + (32.1 * om0h2).powf(-0.532));
        let a2 = (12.0 * om0h2).powf(0.424) * (1.0 + (45.0 * om0h2).powf(-0.582));
        let ac = a1.powf(-ombom0) * a2.powf(-ombom0.powi(3));

        // Equation 12
        let b1 = 0.944 / (1.0 + (458.0 * om0h2).powf(-0.708));
        let b2 = (0.395 * om0h2).powf(-0.0266);
        let bc = 1.0 / (1.0 + b1 * ((omcom0).powf(b2) - 1.0));

        // Equation 15
        let y = (1.0 + zeq) / (1.0 + zd);
        let gy = y
            * (-6.0 * sqrt(1.0 + y)
                + (2.0 + 3.0 * y) * ln((sqrt(1.0 + y) + 1.0) / (sqrt(1.0 + y) - 1.0)));

        // Equation 14
        let ab = 2.07 * keq * s * (1.0 + rd).powf(-3.0 / 4.0) * gy;

        // Get CDM part of transfer function

        // Equation 18
        let f = 1.0 / (1.0 + (kh * s / 5.4).powi(4));

        // Equation 20
        let c = 14.2 / ac + 386.0 / (1.0 + 69.9 * q.powf(1.08));

        // Equation 19
        let t0t = ln(EULER + 1.8 * bc * q) / (ln(EULER + 1.8 * bc * q) + c * q * q);

        // Equation 17
        let c1bc = 14.2 + 386.0 / (1.0 + 69.9 * q.powf(1.08));
        let t0t1bc = ln(EULER + 1.8 * bc * q) / (ln(EULER + 1.8 * bc * q) + c1bc * q * q);
        let transfer_cold = f * t0t1bc + (1.0 - f) * t0t;

        // Get baryon part of transfer function

        // Equation 24
        let bb = 0.5 + ombom0 + (3.0 - 2.0 * ombom0) * sqrt((17.2 * om0h2) * (17.2 * om0h2) + 1.0);

        // Equation 23
        let bnode = 8.41 * om0h2.powf(0.435);

        // Equation 22
        let st = s / (1.0 + (bnode / kh / s).powi(3)).powf(1.0 / 3.0);

        // Equation 21
        let c11 = 14.2 + 386.0 / (1.0 + 69.9 * q.powf(1.08));
        let tot11 = ln(EULER + 1.8 * q) / (ln(EULER + 1.8 * q) + c11 * q * q);
        let transfer_baryon = (tot11 / (1.0 + (kh * s / 5.2).powi(2))
            + ab / (1.0 + (bb / kh / s).powi(3)) * exp(-(kh / ksilk).powf(1.4)))
            * sin(kh * st)
            / (kh * st);

        // Total transfer function
        ombom0 * transfer_baryon + omcom0 * transfer_cold
    }

    // Ignore baryonic effects
    pub fn transfer_zero_baryon(&self, k: f64) -> f64 {
        let ombom0 = self.omega_baryon_0 / self.omega_matter_0;
        let h2 = self.h.powi(2);
        let om0h2 = self.omega_matter_0 * h2;
        let ombh2 = self.omega_baryon_0 * h2;
        let theta2p7 = self.temp_cmb0 / 2.7;

        // Convert kh from hMpc^-1 to Mpc^-1
        let kh = k * self.h;

        // Equation 26
        let s = 44.5 * ln(9.83 / om0h2) / sqrt(1.0 + 10.0 * ombh2.powf(0.75));

        // Equation 31
        let alpha_gamma =
            1.0 - 0.328 * ln(431.0 * om0h2) * ombom0 + 0.38 * ln(22.3 * om0h2) * ombom0.powi(2);

        // Equation 30
        let gamma = self.omega_matter_0
            * self.h
            * (alpha_gamma + (1.0 - alpha_gamma) / (1.0 + (0.43 * kh * s).powi(4)));

        // Equation 28
        let q = k * theta2p7 * theta2p7 / gamma;

        // Equation 29
        let c0 = 14.2 + 731.0 / (1.0 + 62.5 * q);
        let l0 = ln(2.0 * exp(1.0) + 1.8 * q);

        l0 / (l0 + c0 * q * q)
    }

    fn _sigma8_calc_baryon(&self) -> f64 {
        // Calculate sigma squared integral
        const LOGK_MIN: f64 = -20.0;
        const LOGK_MAX: f64 = 20.0;
        const SUB_INTERVALS: usize = 100;
        let sigma_2: f64 = (0..SUB_INTERVALS)
            .map(|i| {
                // Interval bounds
                let min = LOGK_MIN + i as f64 * (LOGK_MAX - LOGK_MIN) / SUB_INTERVALS as f64;
                let max = LOGK_MIN + (i + 1) as f64 * (LOGK_MAX - LOGK_MIN) / SUB_INTERVALS as f64;

                // Integrate over this sub_interval
                let integ = integrate(
                    |logk: f64| self.log_integrand_baryon(logk, 8.0),
                    min,
                    max,
                    1e-8,
                );

                // Return result of integral
                integ.integral
            })
            .sum();

        sqrt(sigma_2 / 2.0 / PI.powi(2))
    }

    fn _sigma8_calc_zero_baryon(&self) -> f64 {
        // Calculate sigma squared integral
        const LOGK_MIN: f64 = -20.0;
        const LOGK_MAX: f64 = 20.0;
        const SUB_INTERVALS: usize = 100;
        let sigma_2: f64 = (0..SUB_INTERVALS)
            .map(|i| {
                // Interval bounds
                let min = LOGK_MIN + i as f64 * (LOGK_MAX - LOGK_MIN) / SUB_INTERVALS as f64;
                let max = LOGK_MIN + (i + 1) as f64 * (LOGK_MAX - LOGK_MIN) / SUB_INTERVALS as f64;

                // Integrate over this sub_interval
                let integ = integrate(
                    |logk: f64| self.log_integrand_zero_baryon(logk, 8.0),
                    min,
                    max,
                    1e-8,
                );

                // Return result of integral
                integ.integral
            })
            .sum();

        sqrt(sigma_2 / 2.0 / PI.powi(2))
    }

    fn log_integrand_baryon(&self, logk: f64, r: f64) -> f64 {
        // k is used several times so just `exp` once.
        let k = logk.exp();

        // Calculate tophat_filter, fourier space
        let weight = self.tophat_filter(k, r);

        // Calculate power
        let power = self.transfer_baryon(k).powi(2) * k.powf(self.ns);

        // Pk * W^2 * k^2 * k dk
        power * weight.powi(2) * k.powi(3)
    }

    fn log_integrand_zero_baryon(&self, logk: f64, r: f64) -> f64 {
        // k is used several times so just `exp` once.
        let k = logk.exp();

        // Calculate tophat_filter, fourier space
        let weight = self.tophat_filter(k, r);

        // Calculate power
        let power = self.transfer_zero_baryon(k).powi(2) * k.powf(self.ns);

        // Pk * W^2 * k^2 * k dk
        power * weight.powi(2) * k.powi(3)
    }

    fn tophat_filter(&self, k: f64, r: f64) -> f64 {
        // Dimensionless product kr
        let x = k * r;

        // Numerical stablity for lim x -> 0
        if x < 1e-3 {
            1.0
        } else {
            3.0 / x.powi(3) * (sin(x) - x * cos(x))
        }
    }
}

fn sqrt(x: f64) -> f64 {
    x.sqrt()
}
fn ln(x: f64) -> f64 {
    x.ln()
}
fn exp(x: f64) -> f64 {
    x.exp()
}
fn sin(x: f64) -> f64 {
    x.sin()
}
fn cos(x: f64) -> f64 {
    x.cos()
}

#[cfg(test)]
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

#[cfg(test)]
#[cfg(not(feature = "colossus-python"))]
/// This set of unit tests use a small subset of hard coded externally
/// produced results from colossus
mod tests {
    use super::EisensteinHu;

    #[test]
    fn test_eisen_hu_transfer_30() {
        // Construct EisensteinHu model
        let eisen_hu = EisensteinHu::new(
            0.7,    // h
            0.3,    // omega_matter_0
            0.025,  // omega_baryon_0
            2.7,    // temp_cmb_0
            0.9665, // ns
            0.8102, // sigma8
        )
        .unwrap();

        // Pick wavenumbers
        let ks = [0.1, 1.0, 10.0, 100.0];

        // Expected values, from COLOSSUS
        let expected = vec![
            0.15333336269193304,
            0.005771855001319983,
            0.00011260658988099492,
            1.6961708546353218e-06,
        ];

        for i in 0..ks.len() {
            assert_eq_tol!(eisen_hu.transfer_baryon(ks[i]), expected[i], 1e-7);
        }
    }

    #[test]
    fn test_eisen_hu_transfer_100_10() {
        // Construct EisensteinHu model
        let eisen_hu = EisensteinHu::new(
            0.7,    // h
            1.0,    // omega_matter_0
            0.1,    // omega_baryon_0
            2.7,    // temp_cmb_0
            0.9665, // ns
            0.8102, // sigma8
        )
        .unwrap();

        // Pick wavenumbers
        let ks = [0.1, 1.0, 10.0, 100.0];

        // Expected values, from COLOSSUS
        let expected = vec![
            0.4338624389792491,
            0.03371743303680279,
            0.000832535760271818,
            1.3932881465826914e-05,
        ];

        for i in 0..ks.len() {
            assert_eq_tol!(eisen_hu.transfer_baryon(ks[i]), expected[i], 1e-7);
        }
    }

    #[test]
    fn test_eisen_hu_transfer_100_zerobaryon() {
        // Construct EisensteinHu model
        let eisen_hu = EisensteinHu::new(
            0.7,    // h
            1.0,    // omega_matter_0
            0.0,    // omega_baryon_0
            2.7,    // temp_cmb_0
            0.9665, // ns
            0.8102, // sigma8
        )
        .unwrap();

        // Pick wavenumbers
        let ks = [0.1, 1.0, 10.0, 100.0];

        // Expected values, from COLOSSUS
        #[cfg(not(feature = "colossus-python"))]
        let expected = vec![
            0.4924960435938662,
            0.043721843089616255,
            0.0011207655979098069,
            1.9111786222704576e-05,
        ];

        for i in 0..ks.len() {
            assert_eq_tol!(eisen_hu.transfer_zero_baryon(ks[i]), expected[i], 1e-7);
        }
    }

    #[test]
    fn test_eisen_hu_transfer_30_zerobaryon() {
        // Construct EisensteinHu model
        let eisen_hu = EisensteinHu::new(
            0.7,    // h
            0.3,    // omega_matter_0
            0.0,    // omega_baryon_0
            2.7,    // temp_cmb_0
            0.9665, // ns
            0.8102, // sigma8
        )
        .unwrap();

        // Pick wavenumbers
        let ks = [0.1, 1.0, 10.0, 100.0];

        // Expected values, from COLOSSUS
        let expected = vec![
            0.17606761051483671,
            0.006943761518101344,
            0.0001377412172820803,
            2.0957565625293516e-06,
        ];

        for i in 0..ks.len() {
            assert_eq_tol!(eisen_hu.transfer_zero_baryon(ks[i]), expected[i], 1e-7);
        }
    }

    #[test]
    fn test_current_power() {
        // Construct EisensteinHu model
        let eisen_hu = EisensteinHu::new(
            0.7,    // h
            0.3,    // omega_matter_0
            0.025,  // omega_baryon_0
            2.7,    // temp_cmb_0
            0.9665, // ns
            0.8102, // sigma8
        )
        .unwrap();

        // Pick wavenumbers
        let ks = [1.0, 10.0, 100.0];

        // Get result at redshift zero
        let result = eisen_hu.power_z0(&ks);

        // Expected values, from COLOSSUS
        let expected = vec![72.208536677773, 0.25444155136438845, 0.0005344388706918762];

        for i in 0..result.len() {
            assert_eq_tol!(result[i], expected[i], 1e-4);
        }
    }

    #[test]
    fn test_nonzero_redshift_power() {
        // Construct EisensteinHu model
        let eisen_hu = EisensteinHu::new(
            0.7,    // h
            0.3,    // omega_matter_0
            0.025,  // omega_baryon_0
            2.7,    // temp_cmb_0
            0.9665, // ns
            0.8102, // sigma8
        )
        .unwrap();

        // Pick wavenumbers
        let ks = [1.0, 10.0, 100.0];

        // Get power at z = 2.0
        let result: Vec<f64> = eisen_hu.power_z(&ks, 2.0).unwrap();

        // Expected values, from COLOSSUS
        let expected = vec![
            12.825704319955586,
            0.045193993046461566,
            9.492721010499203e-05,
        ];

        for i in 0..result.len() {
            assert_eq_tol!(result[i], expected[i], 1e-4);
        }
    }
}

#[cfg(test)]
#[cfg(feature = "colossus-python")]
/// This set of unit tests directly uses the colossus package to test
/// a larger subset of cosmological parameter space. This crate was 
/// originally tested with colossus==1.3.1
mod tests {

    macro_rules! eisenstein_transfer_baryon_test(
    ($h0:ident, $om0:ident, $ob0:ident, $t0:ident) => {

      concat_idents::concat_idents!(test_name = test_eisen_hu_transfer_, $h0, _, $om0, _, $ob0, _, $t0, {
      #[test]
      fn test_name() {

        let h: u32 = stringify!($h0)[1..].parse::<u32>().unwrap();
        let om0: u32 = stringify!($om0)[1..].parse::<u32>().unwrap();
        let ob0: u32 = stringify!($ob0)[1..].parse::<u32>().unwrap();
        let t0: u32 = stringify!($t0)[1..].parse::<u32>().unwrap();

        // Construct EisensteinHu model
        let eisen_hu = super::EisensteinHu::new(
          h as f64 / 100.0, // h
          om0 as f64 / 100.0, // omega_matter_0
          ob0 as f64 / 100.0, // omega_baryon_0
          t0 as f64 / 100.0, // temp_cmb_0
          0.9665, // ns
          0.8102, // sigma8
        ).unwrap();

        // Pick wavenumbers
        let ks = [0.1, 1.0, 10.0, 100.0];

        // Expected values, from COLOSSUS
        let expected = {
          use pyo3::prelude::*;
          use pyo3::types::*;
          Python::with_gil(|py| {

            // Get ks into python
            let list = PyList::new(py, &ks);
            let locals = PyDict::new(py);
            locals.set_item("ks", list).unwrap();

            py.run(format!(r#"from colossus.cosmology import cosmology
import warnings
warnings.filterwarnings("ignore")
x = []
for k in ks:
  x.append(cosmology.power_spectrum.transferFunction(k, {0}, {1}, {2}, {3}, model='eisenstein98'))
            "#, h as f64 / 100.0, om0 as f64 / 100.0,
             ob0 as f64 / 100.0, t0 as f64 / 100.0 ).as_str(), None, Some(locals)).unwrap();
            let x: Vec<_> = locals.get_item("x").unwrap().extract::<Vec<f64>>().unwrap();
            x
          })
        };

        for i in 0..ks.len() {
          assert_eq_tol!(
            eisen_hu.transfer_baryon(ks[i]),
            expected[i],
            1e-7
          );
        }
      }
    });
  });

  macro_rules! eisenstein_power_baryon_test(
    ($h0:ident, $om0:ident, $ob0:ident, $z0:ident) => {

      concat_idents::concat_idents!(test_name = test_eisen_hu_transfer_, $h0, _, $om0, _, $ob0, _, $t0, {
      #[test]
      fn test_name() {

        let h: u32 = stringify!($h0)[1..].parse::<u32>().unwrap();
        let om0: u32 = stringify!($om0)[1..].parse::<u32>().unwrap();
        let ob0: u32 = stringify!($ob0)[1..].parse::<u32>().unwrap();
        let z0: u32 = stringify!($z0)[1..].parse::<u32>().unwrap();

        // Construct EisensteinHu model
        let eisen_hu = super::EisensteinHu::new(
          h as f64 / 100.0, // h
          om0 as f64 / 100.0, // omega_matter_0
          ob0 as f64 / 100.0, // omega_baryon_0
          t0 as f64 / 100.0, // temp_cmb_0
          0.9665, // ns
          0.8102, // sigma8
        ).unwrap();

        // Pick wavenumbers
        let ks = [0.1, 1.0, 10.0, 100.0];

        // Expected values, from COLOSSUS
        let expected = {
          use pyo3::prelude::*;
          use pyo3::types::*;
          Python::with_gil(|py| {

            // Get ks into python
            let list = PyList::new(py, &ks);
            let locals = PyDict::new(py);
            locals.set_item("ks", list).unwrap();

            py.run(format!(r#"from colossus.cosmology import cosmology
import warnings
warnings.filterwarnings("ignore")
my_cosmo = {
    "H0": {0},    
    "Om0": {1},    
    "Ob0": {2},  
    "ns": 0.9665,
    "sigma8": 0.8102, 
}
cosmology.addCosmology("my_cosmo", my_cosmo)
cosmo = cosmology.setCosmology("my_cosmo")
x = []
for k in ks:
  x.append(cosmo.matterPowerSpectrum(k, {z}))
            "#, h as f64, om0 as f64 / 100.0,
             ob0 as f64 / 100.0).as_str(), None, Some(locals)).unwrap();
            let x: Vec<_> = locals.get_item("x").unwrap().extract::<Vec<f64>>().unwrap();
            x
          })
        };

        let power = eisen_hu.power_z(ks, z as f64);

        for i in 0..ks.len() {
          assert_eq_tol!(
            power[i],
            expected[i],
            1e-7
          );
        }
      }
    });
  });

    macro_rules! eisenstein_no_baryon_test(
    ($h0:ident, $om0:ident, $ob0:ident, $t0:ident) => {

      concat_idents::concat_idents!(test_name = test_eisen_hu_transfer_nb_, $h0, _, $om0, _, $ob0, _, $t0, {
      #[test]
      fn test_name() {

        let h: u32 = stringify!($h0)[1..].parse::<u32>().unwrap();
        let om0: u32 = stringify!($om0)[1..].parse::<u32>().unwrap();
        let ob0: u32 = stringify!($ob0)[1..].parse::<u32>().unwrap();
        let t0: u32 = stringify!($t0)[1..].parse::<u32>().unwrap();

        // Construct EisensteinHu model
        let eisen_hu = super::EisensteinHu::new(
          h as f64 / 100.0, // h
          om0 as f64 / 100.0, // omega_matter_0
          ob0 as f64 / 100.0, // omega_baryon_0
          t0 as f64 / 100.0, // temp_cmb_0
          0.9665, // ns
          0.8102, // sigma8
        ).unwrap();

        // Pick wavenumbers
        let ks = [0.1, 1.0, 10.0, 100.0];

        // Expected values, from COLOSSUS
        let expected = {
          use pyo3::prelude::*;
          use pyo3::types::*;
          Python::with_gil(|py| {

            // Get ks into python
            let list = PyList::new(py, &ks);
            let locals = PyDict::new(py);
            locals.set_item("ks", list).unwrap();

            py.run(format!(r#"from colossus.cosmology import cosmology
import warnings
warnings.filterwarnings("ignore")
x = []
for k in ks:
  x.append(cosmology.power_spectrum.transferFunction(k, {0}, {1}, {2}, {3}, model='eisenstein98_zb'))
            "#, h as f64 / 100.0, om0 as f64 / 100.0,
            ob0 as f64 / 100.0, t0 as f64 / 100.0 ).as_str(), None, Some(locals)).unwrap();
            let x: Vec<_> = locals.get_item("x").unwrap().extract::<Vec<f64>>().unwrap();
            x
          })
        };

        for i in 0..ks.len() {
          assert_eq_tol!(
            eisen_hu.transfer_zero_baryon(ks[i]),
            expected[i],
            1e-7
          );
        }
      }
    });
  });

    macro_rules! eisenstein_power(
    ($z:ident, $h0:ident, $om0:ident, $ob0:ident, $t0:ident) => {

      concat_idents::concat_idents!(test_name = test_eisen_power, _, $z, _, $h0, _, $om0, _, $ob0, _, $t0, {
        #[test]
        fn test_name() {

          let z: u32 = stringify!($z)[1..].parse::<u32>().unwrap();
          let h: u32 = stringify!($h0)[1..].parse::<u32>().unwrap();
          let om0: u32 = stringify!($om0)[1..].parse::<u32>().unwrap();
          let ob0: u32 = stringify!($ob0)[1..].parse::<u32>().unwrap();
          let t0: u32 = stringify!($t0)[1..].parse::<u32>().unwrap();

          // Construct EisensteinHu model
          let eisen_hu = super::EisensteinHu::new(
            h as f64 / 100.0, // h
            om0 as f64 / 100.0, // omega_matter_0
            ob0 as f64 / 100.0, // omega_baryon_0
            t0 as f64 / 100.0, // temp_cmb_0
            0.9665, // ns
            0.8102, // sigma8
          ).unwrap();

          // Pick wavenumbers
          let ks = [0.1, 1.0, 10.0, 100.0];

          // Get result at redshift zero
          let result = eisen_hu.power_z(&ks, z as f64).unwrap();

          // Expected values, from COLOSSUS
          let expected = {
            use pyo3::prelude::*;
            use pyo3::types::*;
            Python::with_gil(|py| {

              // Get ks into python
              let list = PyList::new(py, &ks);
              let locals = PyDict::new(py);
              locals.set_item("ks", list).unwrap();

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
for k in ks:
  x.append(cosmo.matterPowerSpectrum(k, z={4}))
              "#, h as f64, om0 as f64 / 100.0, ob0 as f64 / 100.0,
              t0 as f64 / 100.0, z).as_str(), None, Some(locals)).unwrap();
              let x: Vec<_> = locals.get_item("x").unwrap().extract::<Vec<f64>>().unwrap();
              x
            })
          };

          for i in 0..result.len() {
            assert_eq_tol!(result[i], expected[i], 1e-2);
          }
        }
      });
    }
  );

    macro_rules! eisenstein_power_no_baryon(
    ($z:ident, $h0:ident, $om0:ident, $ob0:ident, $t0:ident) => {

      concat_idents::concat_idents!(test_name = test_eisen_power_no_baryon, _, $z, _, $h0, _, $om0, _, $ob0, _, $t0, {
        #[test]
        fn test_name() {

          let z: u32 = stringify!($z)[1..].parse::<u32>().unwrap();
          let h: u32 = stringify!($h0)[1..].parse::<u32>().unwrap();
          let om0: u32 = stringify!($om0)[1..].parse::<u32>().unwrap();
          let ob0: u32 = stringify!($ob0)[1..].parse::<u32>().unwrap();
          let t0: u32 = stringify!($t0)[1..].parse::<u32>().unwrap();

          // Construct EisensteinHu model
          let eisen_hu = super::EisensteinHu::new(
            h as f64 / 100.0, // h
            om0 as f64 / 100.0, // omega_matter_0
            ob0 as f64 / 100.0, // omega_baryon_0
            t0 as f64 / 100.0, // temp_cmb_0
            0.9665, // ns
            0.8102, // sigma8
          ).unwrap();

          // Pick wavenumbers
          let ks = [0.1, 1.0, 10.0, 100.0];

          // Get result at redshift zero
          let result = eisen_hu.power_z_zero_baryon(&ks, z as f64).unwrap();

          // Expected values, from COLOSSUS
          let expected = {
            use pyo3::prelude::*;
            use pyo3::types::*;
            Python::with_gil(|py| {

              // Get ks into python
              let list = PyList::new(py, &ks);
              let locals = PyDict::new(py);
              locals.set_item("ks", list).unwrap();

              py.run(format!(r#"from colossus.cosmology import cosmology
import warnings
warnings.filterwarnings("ignore")
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
for k in ks:
  x.append(cosmo.matterPowerSpectrum(k, z={4}, model="eisenstein98_zb"))
              "#, h as f64, om0 as f64 / 100.0, ob0 as f64 / 100.0, t0 as f64 / 100.0, z).as_str(), None, Some(locals)).unwrap();
              let x: Vec<_> = locals.get_item("x").unwrap().extract::<Vec<f64>>().unwrap();
              x
            })
          };

          for i in 0..result.len() {
            assert_eq_tol!(result[i], expected[i], 1e-2);
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
                        eisenstein_power!($Z, $H, $M, $B, $T);
                        eisenstein_power_no_baryon!($Z, $H, $M, $B, $T);
                    });
                });
            });
        });
    });


    dry::macro_for!($H in [h50, h60, h70, h80, h90, h100] {
        dry::macro_for!($M in [m10, m30, m50, m70, m90, m100] {
            dry::macro_for!($B in [b1, b2, b3] {
                dry::macro_for!($T in [t268, t270, t272] {
                    eisenstein_transfer_baryon_test!($H, $M, $B, $T);
                    eisenstein_no_baryon_test!($H, $M, $B, $T);
                });
            });
        });
    });
}
