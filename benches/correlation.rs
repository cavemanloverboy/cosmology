

use std::ops::Sub;

use cosmology::correlation::{CorrelationFunctionParameters, CorrelationFunction};
use criterion::{Criterion, criterion_group, criterion_main};
use rand_distr::num_traits::Pow;


const LOGR_MIN: f64 = 0.0;
const LOGR_MAX: f64 = 2.0;
const NUM_R_POINTS: usize = 200;

fn bench_corr(c: &mut Criterion) {

    // Initialize power engine
    let power = cosmology::power::PowerSpectrum::new(
        cosmology::power::TransferFunction::EisensteinHu {
            h: 0.7,
            omega_matter_0: 0.3,
            omega_baryon_0: 0.2, 
            temp_cmb0: 2.7,
            ns: 0.96,
            sigma_8: 0.6,
        }
    ).unwrap();
    println!("initialized power engine");

    // Initialize correlation parameters
    let params = CorrelationFunctionParameters {
        power,
        accuracy_params: None,
    };
    let z = 1.0;
    let correlation_function = CorrelationFunction::get_correlation_function(z, params).unwrap();
    println!("initialized correlation function");


    // Get scales of interest
    let rstep = (LOGR_MAX-LOGR_MIN)/(NUM_R_POINTS.sub(1) as f64);
    let rs: Vec<f64> = (0..NUM_R_POINTS)
        .map(|i| 10_f64.powf(LOGR_MIN + i as f64 * rstep))
        .collect();
    c.bench_function("correlation", move |b| {
        b.iter(|| {
            for r in &rs {
                correlation_function.correlation_function(*r);
            }
        })
    });
    
}




criterion_group!(power, bench_corr);
criterion_main!(power);