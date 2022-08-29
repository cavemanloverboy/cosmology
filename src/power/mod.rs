use transfer::eisenstein::{self, *};

use self::growth_factor::linear_growth_factor;

pub mod transfer;
pub mod growth_factor;


pub struct PowerSpectrum(TransferFunctionEngine);

pub enum TransferFunction {
    EisensteinHu {
        h: f64,
        omega_matter_0: f64,
        omega_baryon_0: f64,
        temp_cmb0: f64,
        ns: f64,
        sigma_8: f64,
    },
    EisensteinHuNoBaryon {
        h: f64,
        omega_matter_0: f64,
        omega_baryon_0: f64,
        temp_cmb0: f64,
        ns: f64,
        sigma_8: f64,
    }
}

enum TransferFunctionEngine {
    EisensteinHu(EisensteinHu),
    EisensteinHuNoBaryon(EisensteinHu),
}

impl PowerSpectrum {


    pub fn new(
        transfer: TransferFunction,
    ) -> Result<Self, &'static str> {
        match transfer {

            // Eistenstein & Hu 1998
            TransferFunction::EisensteinHu { h, omega_matter_0, omega_baryon_0, temp_cmb0, ns, sigma_8 } => {
                Ok(PowerSpectrum(
                    TransferFunctionEngine::EisensteinHu(
                        EisensteinHu::new(h, omega_matter_0, omega_baryon_0, temp_cmb0, ns, sigma_8)?
                    )
                ))
            },

            // Eistenstein & Hu 1998 (zero baryon, i.e. no BAO)
            // This uses the same underlying engine as EisensteinHu variant, but calls different methods. 
            TransferFunction::EisensteinHuNoBaryon { h, omega_matter_0, omega_baryon_0, temp_cmb0, ns, sigma_8 } => {
                Ok(PowerSpectrum(
                    TransferFunctionEngine::EisensteinHuNoBaryon(
                        EisensteinHu::new(h, omega_matter_0, omega_baryon_0, temp_cmb0, ns, sigma_8)?
                    )
                ))
            }
        }
    }

    pub fn calculate_power(
        &self,
        ks: &[f64],
        z: f64,
    ) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
        
        match &self.0 {   

            // Eistenstein & Hu 1998
            TransferFunctionEngine::EisensteinHu(e_hu_engine) => {
                e_hu_engine.power_z(ks, z)
            },

            // Eistenstein & Hu 1998 (zero baryon, i.e. no BAO)
            TransferFunctionEngine::EisensteinHuNoBaryon(e_hu_engine) => {
                Ok(e_hu_engine.power_z0_zero_baryon(ks))
            }

        }
        
    }
}






#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_create_power_spectrum() {

        let _ = PowerSpectrum::new(
            TransferFunction::EisensteinHu {
                h: 0.7,
                omega_matter_0: 0.3,
                omega_baryon_0: 0.02,
                temp_cmb0: 2.7,
                ns: 1.0,
                sigma_8: 0.06
            }
        );
    }
}