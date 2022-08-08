
use num::Float;
use special::Gamma;
use std::{f64::consts::{E as EULER, PI}, thread::current};

// SPEED OF LIGHT IN UNITS OF Mpc/s
const C_MPC_S: f64 = 299_792_458.0 / 3.08567758128e22;
const C_KM_S: f64 = 299_792.458;

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

    pub fn new(
      h: f64,
      omega_matter_0: f64,
      omega_baryon_0: f64,
      temp_cmb0: f64,
      ns: f64,
      sigma_8: f64,
    ) -> Result<EisensteinHu, &'static str> {
        if omega_baryon_0 <= std::f64::MIN {
            return Err("The Eisenstein & Hu 98 transfer function cannot be computed for Ob0 = 0")
        } else if temp_cmb0 <= 0.0 {
            return Err("Cannot have a nonpositive CMB temperature")
        } else if omega_baryon_0 > omega_matter_0 {
            return Err("Cannot have more baryons than total matter")
        } else {
          Ok(EisensteinHu {
            h,
            omega_matter_0,
            omega_baryon_0,
            temp_cmb0,
            ns,
            sigma_8
          })
        }
    }

    /// Transfer function at wavenumber k (Eisenstein & Hu 1998).
    /// The code was adapted from Benedikt Diemer's COLOSSUS code,
    /// which was adapted from Matt Becker's cosmocalc code.
    /// 
    /// Given a slice of wavenumbers `&[T]`, returns the transfer function
    /// with the same units as `k`.
    pub fn power(&self, ks: &[f64], z: f64) -> Vec<f64> {

      ks
        .iter()
        .map(|&k| {
          

          // // delta_h
          // let delta_h = {

          //   // Auxilary value: difference from 1
          //   let n_tilde = self.ns - 1.0;

          //   if self.omega_matter_0 == 1.0 {
          //     1.95e-5 * exp(-n_tilde - 0.14*n_tilde.powi(2))
          //   } else {
          //     1.94e-5 * self.omega_matter_0.powf(-0.785 - 0.05 * ln(self.omega_matter_0))
          //       * exp(-0.95*n_tilde - 0.169*n_tilde.powi(2))
          //   }
          // };
          
          // Calculate current power
          let norm = sigma2_power_law_tophat(8.0);
          let current_power = self.sigma_8.powi(2) * self.transfer(k).powi(2) * k.powf(self.ns) / norm.powi(2);
          // let current_power = 2.0 * PI.powi(2) * delta_h.powi(2) 
          //   * (C_KM_S / 100.0).powf(3.0 + self.ns) * k.powf(self.ns) * transfer.powi(2);


          // // Calculate linear growth factor
          // let dz = |z: f64| {

          //   if self.omega_matter_0 == 1.0 {
          //     (1.0 + z).recip()
          //   } else {

          //     // Auxilary values
          //     let oml = 1.0 - self.omega_matter_0;
          //     let omz = self.omega_matter_0 * (1.0 + z).powi(3)
          //       / (oml + self.omega_matter_0 * (1.0 + z).powi(3));

              
          //     (1.0 + z).recip() * 5.0 * omz / 2.0 
          //       / (omz.powf(4.0/7.0) - oml + (1.0 + omz/2.0)*(1.0 + oml/70.0))
          //   }
          // };

          current_power //* dz(z).powi(2) / dz(0.0).powi(2)
        }).collect::<Vec<_>>()
    }


    /// Transfer function. Tested against colossus 
    fn transfer(&self, k: f64) -> f64 {
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
      let s = 2.0 / 3.0 / keq * sqrt(6.0 / req) * ln((sqrt(1.0 + rd)
        + sqrt(rd + req)) / (1.0 + sqrt(req)));

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
      let gy = y * (-6.0 * sqrt(1.0 + y) + (2.0 + 3.0 * y)
        * ln((sqrt(1.0 + y) + 1.0) / (sqrt(1.0 + y) - 1.0)));

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
        * sin(kh * st) / (kh * st);

      // Total transfer function
      let transfer = ombom0 * transfer_baryon + omcom0 * transfer_cold;

      transfer
    }

    // Ignore baryonic effects
    fn transfer_zero_baryon(&self, k: f64) -> f64 {

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
      let alpha_gamma = 1.0 - 0.328 * ln(431.0 * om0h2) * ombom0 
        + 0.38 * ln(22.3 * om0h2) * ombom0.powi(2);
    
      // Equation 30
      let gamma = self.omega_matter_0 * self.h * (alpha_gamma + (1.0 - alpha_gamma) 
        / (1.0 + (0.43 * kh * s).powi(4)));
    
      // Equation 28
      let q = k * theta2p7 * theta2p7 / gamma;
    
      // Equation 29
      let c0 = 14.2 + 731.0 / (1.0 + 62.5 * q);
      let l0 = ln(2.0 * exp(1.0) + 1.8 * q);
      let transfer = l0 / (l0 + c0 * q * q);
    
      transfer
    }
}

fn sqrt(x: f64) -> f64 { x.sqrt() }
fn ln(x: f64) -> f64 { x.ln() }
fn exp(x: f64) -> f64 { x.exp() }
fn sin(x: f64) -> f64 { x.sin() }
fn cos(x: f64) -> f64 { x.cos() }

fn sigma2_power_law_tophat(R: f64) -> f64 {

      let n = 2.0;
			
			let mut ret = 9.0 * R.powf(-n - 3.0) * 2.0.powf(-n - 1.0) / (PI.powi(2) * (n - 3.0));

			if (n + 2.0).abs() < 1E-3 {
				ret *= -(1.0 + n) * PI / (2.0 * (2.0 - n).gamma() * cos(PI * n * 0.5));
      } else {
				ret *= sin(n * PI * 0.5) * (n + 2.0).gamma() / ((n - 1.0) * n);
      }
		
			return ret
}


#[allow(unused)]
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

  fn test_eisen_hu_transfer_30() {

    // Construct EisensteinHu model
    let eisen_hu = EisensteinHu::new(
      0.7, // h
      0.3, // omega_matter_0
      0.025, // omega_baryon_0
      2.7, // temp_cmb_0
      0.9665, // ns
      0.8102, // sigma8
    ).unwrap();

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
      assert_eq_tol!(
        eisen_hu.transfer(ks[i]),
        expected[i],
        1e-7
      );
    }
  }

  #[test]
  fn test_eisen_hu_transfer_100_10() {

    // Construct EisensteinHu model
    let eisen_hu = EisensteinHu::new(
      0.7, // h
      1.0, // omega_matter_0
      0.1, // omega_baryon_0
      2.7, // temp_cmb_0
      0.9665, // ns
      0.8102, // sigma8
    ).unwrap();

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
      assert_eq_tol!(
        eisen_hu.transfer(ks[i]),
        expected[i],
        1e-7
      );
    }
  }

  #[test]
  fn test_eisen_hu_transfer_100_zerobaryon() {

    // Construct EisensteinHu model
    let eisen_hu = EisensteinHu::new(
      0.7, // h
      1.0, // omega_matter_0
      0.0, // omega_baryon_0
      2.7, // temp_cmb_0
      0.9665, // ns
      0.8102, // sigma8
    ).unwrap();

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
      assert_eq_tol!(
        eisen_hu.transfer_zero_baryon(ks[i]),
        expected[i],
        1e-7
      );
    }
  }

  #[test]
  fn test_eisen_hu_transfer_30_zerobaryon() {

    // Construct EisensteinHu model
    let eisen_hu = EisensteinHu::new(
      0.7, // h
      0.3, // omega_matter_0
      0.0, // omega_baryon_0
      2.7, // temp_cmb_0
      0.9665, // ns
      0.8102, // sigma8
    ).unwrap();

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
      assert_eq_tol!(
        eisen_hu.transfer_zero_baryon(ks[i]),
        expected[i],
        1e-7
      );
    }
  }


  #[test]
  #[ignore]
  fn test_current_power() {

    // Construct EisensteinHu model
    let eisen_hu = EisensteinHu::new(
          0.7, // h
          0.3, // omega_matter_0
          0.025, // omega_baryon_0
          2.7, // temp_cmb_0
          0.9665, // ns
          0.8102, // sigma8
      ).unwrap();

      // Pick wavenumbers
      let ks = [1.0, 10.0, 100.0];

      // Get result at redshift zero
      let result = eisen_hu.power(&ks, 0.0);
      println!("{result:?}");

      // Expected values, from COLOSSUS 
      let expected = vec![
        72.208536677773,
        0.25444155136438845,
        0.0005344388706918762,
      ];

      for i in 0..result.len() {
        println!("{i}");
        assert_eq_tol!(result[i], expected[i], 1e-6);
      }
  }

  #[test]
  #[ignore]
  fn test_nonzero_redshift_power() {

    // Construct EisensteinHu model
    let eisen_hu = EisensteinHu::new(
          0.7, // h
          0.3, // omega_matter_0
          0.025, // omega_baryon_0
          2.7, // temp_cmb_0
          0.9665, // ns
          0.8102, // sigma8
      ).unwrap();

      // Pick wavenumbers
      let ks = [1.0, 10.0, 100.0];

      // Get result at redshift zero
      let result = eisen_hu.power(&ks, 2.0);
      println!("{result:?}");

      // Expected values, from COLOSSUS
      let expected = vec![
        12.825704319955586,
        0.045193993046461566,
        9.492721010499203e-05,
      ];

      for i in 0..result.len() {
        assert_eq_tol!(result[i], expected[i], 2e-1);
      }
  }
}

#[cfg(test)]
#[cfg(feature = "colossus-python")]
/// This set of unit tests directly uses the colossus package to test
/// a larger subset of cosmological parameter space
mod tests {

  macro_rules! eisenstein_baryon_test(
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
            "#, h as f64 / 100.0, om0 as f64 / 100.0, ob0 as f64 / 100.0, t0 as f64 / 100.0 ).as_str(), None, Some(locals)).unwrap();
            let x: Vec<_> = locals.get_item("x").unwrap().extract::<Vec<f64>>().unwrap();
            x
          })
        };

        for i in 0..ks.len() {
          assert_eq_tol!(
            eisen_hu.transfer(ks[i]),
            expected[i],
            1e-7
          );
        }
      }
    });
  });

  macro_rules! eisenstein_no_baryon_test(
    ($h0:ident, $om0:ident, $t0:ident) => {
      
      concat_idents::concat_idents!(test_name = test_eisen_hu_transfer_, $h0, _, $om0, _, $t0, {
      #[test]
      fn test_name() {

        let h: u32 = stringify!($h0)[1..].parse::<u32>().unwrap();
        let om0: u32 = stringify!($om0)[1..].parse::<u32>().unwrap();
        let t0: u32 = stringify!($t0)[1..].parse::<u32>().unwrap();

        // Construct EisensteinHu model
        let eisen_hu = super::EisensteinHu::new(
          h as f64 / 100.0, // h
          om0 as f64 / 100.0, // omega_matter_0
          0.0, // omega_baryon_0
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
  x.append(cosmology.power_spectrum.transferFunction(k, {0}, {1}, 0.0, {2}, model='eisenstein98_zb'))
            "#, h as f64 / 100.0, om0 as f64 / 100.0, t0 as f64 / 100.0 ).as_str(), None, Some(locals)).unwrap();
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

  eisenstein_baryon_test!(h50, m100, b10, t270);
  eisenstein_baryon_test!(h60, m100, b10, t270);
  eisenstein_baryon_test!(h70, m100, b10, t270);
  eisenstein_baryon_test!(h80, m100, b10, t270);
  eisenstein_baryon_test!(h90, m100, b10, t270);

  eisenstein_baryon_test!(h50, m75, b10, t270);
  eisenstein_baryon_test!(h60, m75, b10, t270);
  eisenstein_baryon_test!(h70, m75, b10, t270);
  eisenstein_baryon_test!(h80, m75, b10, t270);
  eisenstein_baryon_test!(h90, m75, b10, t270);

  eisenstein_baryon_test!(h50, m50, b10, t270);
  eisenstein_baryon_test!(h60, m50, b10, t270);
  eisenstein_baryon_test!(h70, m50, b10, t270);
  eisenstein_baryon_test!(h80, m50, b10, t270);
  eisenstein_baryon_test!(h90, m50, b10, t270);

  eisenstein_baryon_test!(h50, m25, b10, t270);
  eisenstein_baryon_test!(h60, m25, b10, t270);
  eisenstein_baryon_test!(h70, m25, b10, t270);
  eisenstein_baryon_test!(h80, m25, b10, t270);
  eisenstein_baryon_test!(h90, m25, b10, t270);

  eisenstein_baryon_test!(h50, m100, b2, t270);
  eisenstein_baryon_test!(h60, m100, b2, t270);
  eisenstein_baryon_test!(h70, m100, b2, t270);
  eisenstein_baryon_test!(h80, m100, b2, t270);
  eisenstein_baryon_test!(h90, m100, b2, t270);

  eisenstein_baryon_test!(h50, m75, b2, t270);
  eisenstein_baryon_test!(h60, m75, b2, t270);
  eisenstein_baryon_test!(h70, m75, b2, t270);
  eisenstein_baryon_test!(h80, m75, b2, t270);
  eisenstein_baryon_test!(h90, m75, b2, t270);

  eisenstein_baryon_test!(h50, m50, b2, t270);
  eisenstein_baryon_test!(h60, m50, b2, t270);
  eisenstein_baryon_test!(h70, m50, b2, t270);
  eisenstein_baryon_test!(h80, m50, b2, t270);
  eisenstein_baryon_test!(h90, m50, b2, t270);
  
  eisenstein_baryon_test!(h50, m25, b2, t270);
  eisenstein_baryon_test!(h60, m25, b2, t270);
  eisenstein_baryon_test!(h70, m25, b2, t270);
  eisenstein_baryon_test!(h80, m25, b2, t270);
  eisenstein_baryon_test!(h90, m25, b2, t270);

  eisenstein_baryon_test!(h50, m100, b10, t265);
  eisenstein_baryon_test!(h60, m100, b10, t265);
  eisenstein_baryon_test!(h70, m100, b10, t265);
  eisenstein_baryon_test!(h80, m100, b10, t265);
  eisenstein_baryon_test!(h90, m100, b10, t265);

  eisenstein_baryon_test!(h50, m75, b10, t265);
  eisenstein_baryon_test!(h60, m75, b10, t265);
  eisenstein_baryon_test!(h70, m75, b10, t265);
  eisenstein_baryon_test!(h80, m75, b10, t265);
  eisenstein_baryon_test!(h90, m75, b10, t265);

  eisenstein_baryon_test!(h50, m50, b10, t265);
  eisenstein_baryon_test!(h60, m50, b10, t265);
  eisenstein_baryon_test!(h70, m50, b10, t265);
  eisenstein_baryon_test!(h80, m50, b10, t265);
  eisenstein_baryon_test!(h90, m50, b10, t265);

  eisenstein_baryon_test!(h50, m25, b10, t265);
  eisenstein_baryon_test!(h60, m25, b10, t265);
  eisenstein_baryon_test!(h70, m25, b10, t265);
  eisenstein_baryon_test!(h80, m25, b10, t265);
  eisenstein_baryon_test!(h90, m25, b10, t265);

  eisenstein_baryon_test!(h50, m100, b2, t265);
  eisenstein_baryon_test!(h60, m100, b2, t265);
  eisenstein_baryon_test!(h70, m100, b2, t265);
  eisenstein_baryon_test!(h80, m100, b2, t265);
  eisenstein_baryon_test!(h90, m100, b2, t265);

  eisenstein_baryon_test!(h50, m75, b2, t265);
  eisenstein_baryon_test!(h60, m75, b2, t265);
  eisenstein_baryon_test!(h70, m75, b2, t265);
  eisenstein_baryon_test!(h80, m75, b2, t265);
  eisenstein_baryon_test!(h90, m75, b2, t265);

  eisenstein_baryon_test!(h50, m50, b2, t265);
  eisenstein_baryon_test!(h60, m50, b2, t265);
  eisenstein_baryon_test!(h70, m50, b2, t265);
  eisenstein_baryon_test!(h80, m50, b2, t265);
  eisenstein_baryon_test!(h90, m50, b2, t265);
  
  eisenstein_baryon_test!(h50, m25, b2, t265);
  eisenstein_baryon_test!(h60, m25, b2, t265);
  eisenstein_baryon_test!(h70, m25, b2, t265);
  eisenstein_baryon_test!(h80, m25, b2, t265);
  eisenstein_baryon_test!(h90, m25, b2, t265);

  eisenstein_no_baryon_test!(h50, m100, t270);
  eisenstein_no_baryon_test!(h60, m100, t270);
  eisenstein_no_baryon_test!(h70, m100, t270);
  eisenstein_no_baryon_test!(h80, m100, t270);
  eisenstein_no_baryon_test!(h90, m100, t270);

  eisenstein_no_baryon_test!(h50, m75, t270);
  eisenstein_no_baryon_test!(h60, m75, t270);
  eisenstein_no_baryon_test!(h70, m75, t270);
  eisenstein_no_baryon_test!(h80, m75, t270);
  eisenstein_no_baryon_test!(h90, m75, t270);

  eisenstein_no_baryon_test!(h50, m50, t270);
  eisenstein_no_baryon_test!(h60, m50, t270);
  eisenstein_no_baryon_test!(h70, m50, t270);
  eisenstein_no_baryon_test!(h80, m50, t270);
  eisenstein_no_baryon_test!(h90, m50, t270);

  eisenstein_no_baryon_test!(h50, m25, t270);
  eisenstein_no_baryon_test!(h60, m25, t270);
  eisenstein_no_baryon_test!(h70, m25, t270);
  eisenstein_no_baryon_test!(h80, m25, t270);
  eisenstein_no_baryon_test!(h90, m25, t270);

  eisenstein_no_baryon_test!(h50, m100, t265);
  eisenstein_no_baryon_test!(h60, m100, t265);
  eisenstein_no_baryon_test!(h70, m100, t265);
  eisenstein_no_baryon_test!(h80, m100, t265);
  eisenstein_no_baryon_test!(h90, m100, t265);

  eisenstein_no_baryon_test!(h50, m75, t265);
  eisenstein_no_baryon_test!(h60, m75, t265);
  eisenstein_no_baryon_test!(h70, m75, t265);
  eisenstein_no_baryon_test!(h80, m75, t265);
  eisenstein_no_baryon_test!(h90, m75, t265);

  eisenstein_no_baryon_test!(h50, m50, t265);
  eisenstein_no_baryon_test!(h60, m50, t265);
  eisenstein_no_baryon_test!(h70, m50, t265);
  eisenstein_no_baryon_test!(h80, m50, t265);
  eisenstein_no_baryon_test!(h90, m50, t265);

  eisenstein_no_baryon_test!(h50, m25, t265);
  eisenstein_no_baryon_test!(h60, m25, t265);
  eisenstein_no_baryon_test!(h70, m25, t265);
  eisenstein_no_baryon_test!(h80, m25, t265);
  eisenstein_no_baryon_test!(h90, m25, t265);
}