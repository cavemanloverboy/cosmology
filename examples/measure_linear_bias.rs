
use std::ops::Div;

use ndarray::Array2;
use ndarray_npy::NpzReader;

use cosmology::{
    power::{PowerSpectrum, TransferFunction},
    correlation::{CorrelationFunctionParameters, CorrelationFunction},
    nn::grf::{SpaceMode, GaussianRandomField},
    bayesian::parameter_inference_uniform_prior,
};
use rand::{thread_rng, Rng};

const REDSHIFT: f64 = 1.0;

fn main() {

    // 100,000 objects per (1000 Mpc/h)^3

    // Scales of interest
    // Get Quijote data
    const QUIJOTE_BOXSIZE: [f64; 3] = [1000.0; 3];
    const N_QUERY: usize = 100;
    let mut npz = NpzReader::new(std::fs::File::open("/home/cavey/experiments/cosmology/OneBox/OneBox.npz").unwrap()).unwrap();
    let real_data: Array2<f32> = npz.by_name("rpos.npy").unwrap();
    let real_data: Array2<f64> = real_data.map(|x| *x as f64);
    let nbar = real_data.len().div(3) as f64 / QUIJOTE_BOXSIZE[0].powi(3);
    println!("{}", real_data.len() / 3);

    // Find 1nn
    let tree = fnntw::Tree::<'_, 3>::new(
            unsafe { core::slice::from_raw_parts(real_data.as_slice().unwrap().as_ptr().cast(), real_data.len() / 3) },
            32
        ).unwrap()
        .with_boxsize(&QUIJOTE_BOXSIZE).unwrap();
    let mut rng = thread_rng();
    let mut nns: Vec<f64> = Vec::with_capacity(N_QUERY);
    for _ in 0..N_QUERY {
        let query: [f64; 3] = rng.gen::<[f64; 3]>().map(|x| x*QUIJOTE_BOXSIZE[0]);
        nns.push(tree.query_nearest(&query).unwrap().0.sqrt());
    }
    nns.sort_by(|a,b| a.partial_cmp(&b).expect("there should be no NaNs"));
    println!("Calculated NNs");
        

    // Initialize E & Hu power spectrum
    // Quijote Cosmology
    let omega_matter_0 = 0.3175;
    let omega_baryon_0 = 0.049;
    let h = 0.6711;
    let ns = 0.9624;
    let sigma_8 = 0.834;
    let power = PowerSpectrum::new(
        TransferFunction::EisensteinHu {
            h,
            omega_matter_0,
            omega_baryon_0,
            temp_cmb0: 2.7255, 
            ns,
            sigma_8,
        }
    ).unwrap();
    let params = CorrelationFunctionParameters {
        power,
        accuracy_params: None, // default accuracy parameters
    };
    let real_corr = CorrelationFunction::get_correlation_function(
        REDSHIFT,
        params
    ).unwrap();
    let real_corr_fn = move |r| real_corr.correlation_function(r);
    let mode = SpaceMode::RealSpace(&real_corr_fn);
    let grf = GaussianRandomField::new(nbar, mode)
        .with(&nns);
    println!("constructed grf. Doing bayesian inference");


    // Construct likelihood function
    let likelihood = move |_distance: &[f64; 0], parameters: &[f64; 1]| -> f64 {
        grf.get_pdf(1, Some(parameters[0]))
            .into_iter()
            .map(|x| x.ln())
            .sum()
    };
    let bounds = [[1.0, 30.0]];
    let most_likely = parameter_inference_uniform_prior::<_,0,1,8,1000,10000>(
        &[[]],
        &bounds,
        likelihood,
    );
    println!("{:?}", most_likely);

}