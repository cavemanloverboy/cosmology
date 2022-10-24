use std::{sync::Arc, ops::Div};

use assert_approx_eq::assert_approx_eq;


type Likelihood<const P: usize> = dyn Fn([f64; P]) -> f64 + Send + Sync;


const INITIAL_RATE: f64 = 1e-4;
const DECAY: f64 = 0.9999;
const DLOG: f64 = 1e-8;

fn rate(iter: i32) -> f64 {
    INITIAL_RATE * DECAY.powi(iter)
}
pub struct Parameters<const P: usize> {
    pub parameters: [f64; P],
    pub initial: [f64; P],
    pub likelihood: Arc<Likelihood<P>>,
    iter: i32,
}


impl<const P: usize> Parameters<P> {

    pub fn new(
        initial: [f64; P],
        likelihood: Arc<Likelihood<P>>
    ) -> Parameters<P> {
        Parameters {
            parameters: initial,
            initial,
            likelihood,
            iter: 0,
        }
    }

    /// TODO: update so that iters are additional iters, not total.
    /// or... think more about api
    pub fn update(&mut self, iters: i32) -> [f64; P] {

        while self.iter < iters {
        
            // Calculate derivative
            let derivatives: [f64; P] = self.parameters
                .iter()
                .enumerate()
                .map(|(i, &p)| {
                    
                    // Copy parameters and displace
                    let mut hi_params: [f64; P] = self.parameters;
                    hi_params[i] = if p != 0.0 {
                        // Fractional change if nonzero
                        p * (1.0 + DLOG) / 1.0 + std::f64::MIN_POSITIVE
                    } else {
                        // Displace a number if it is zero
                        1e-4
                    };
                    let dparam: f64 = hi_params[i] - self.parameters[i];

                    // Calculate likelihood
                    let hi_likelihood: f64 = (self.likelihood)(hi_params);
                    let likelihood: f64 = (self.likelihood)(self.parameters);
                    let dl = hi_likelihood - likelihood;

                    // Calculate derivative
                    dl / dparam
                })
                .collect::<Vec<f64>>()
                .try_into()
                .unwrap();

            // Calculate step size
            let step = derivatives.map(|dldp| rate(self.iter) * dldp);

            // Apply step
            let current_params = self.parameters;
            self.parameters
                .iter_mut()
                .zip(step)
                .enumerate()
                .for_each(|(i, (p, s))| {

                    // Check that the likelihood is fine if this update occurs
                    let mut new_params = current_params;
                    new_params[i] = *p + s;
                    let ll = (self.likelihood)(current_params);
                    let new_ll = (self.likelihood)(new_params);

                    if new_ll.is_normal() && new_ll > ll {
                        *p = *p + s;
                    } else { /* do nothing */ } 
                });

            // Increment
            self.iter += 1;

            // TODO DEBUG
            println!("{:?} -> {}", self.parameters, (self.likelihood)(self.parameters));
        }

        self.parameters
    }
}




#[test]
fn test_sgd() {
    use rand::Rng;

    const BATCH: usize = 32;

    // Underlying mock model parameters
    let m = 0.3;
    let b = 0.5;
    let z = 0.1;
    let mock = [m, b, z];

    let gen_batch = move || {
        let x: [f64; BATCH] = rand::random();
        let y: [f64; BATCH] = x.map(|xx| 
                m * xx + b + z * rand::thread_rng()
                    .sample::<f64, _>(rand_distr::StandardNormal)
            );
        (x, y)
    };

    // Likelihood
    let likelihood = Arc::new(
        move |[mm, bb, zz]: [f64; 3]| -> f64 {
            let (x, y) = gen_batch();

            y.into_iter()
                .zip(x)
                .fold(0.0, |acc, (yi, xi)| acc - zz.ln() - 0.5 * ((yi - (mm * xi + bb)) / zz).powi(2) )
                .div(BATCH as f64)
        }
    );

    // Construct parameters
    let initial_guess = [0.5, 0.2, 3.0];
    let mut parameters = Parameters::new(initial_guess, likelihood);

    parameters.update(100_000);

    parameters
        .parameters
        .iter()
        .enumerate()
        .for_each(|(i, p)| {
            assert_approx_eq!(p, mock[i]);
        })
}