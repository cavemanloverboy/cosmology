
use std::ops::Sub;
use itertools::Itertools;
use plotters::prelude::*;
use ndarray::{Array2, Axis};
use ndarray_npy::NpzReader;

use cosmology::{
    power::{PowerSpectrum, TransferFunction},
    correlation::{CorrelationFunctionParameters, CorrelationFunction},
    nn::grf::{SpaceMode, GaussianRandomField},
    bayesian::parameter_inference_uniform_prior,
};
use rand::{thread_rng, Rng, seq::SliceRandom};

const REDSHIFT: f64 = 1.0;

fn main() {

    // 100,000 objects per (1000 Mpc/h)^3

    // Scales of interest
    // Get Quijote data
    const NDATA: usize = 10_000;
    const QUIJOTE_BOXSIZE: [f64; 3] = [1000.0; 3];
    const N_QUERY: usize = 40_000;
    let mut npz = NpzReader::new(std::fs::File::open("/home/cavey/experiments/cosmology/OneBox/OneBox.npz").unwrap()).unwrap();
    let real_data: Array2<f32> = npz.by_name("rpos.npy").unwrap();
    let real_data: Array2<f64> = real_data.map(|x| *x as f64);

    // subsample
    let mut real_data: Vec<[f64; 3]> = real_data.map_axis(Axis(1), |x| x.as_slice().unwrap().try_into().unwrap()).into_raw_vec();
    (&mut real_data).shuffle(&mut thread_rng());
    real_data.truncate(NDATA);
    let nbar = NDATA as f64 / QUIJOTE_BOXSIZE[0].powi(3);

    // Find 1nn
    println!("Calculating NNs");
    let tree = fnntw::Tree::<'_, 3>::new(
            &real_data,
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
    let nns = nns; // make immutable
    println!("Calculated NNs. Constructing GRF Theory");

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
    println!("Constructed GRF theory. Carrying out Bayesian inference");

    // Construct likelihood function
    let bounds = [[0.2, 10.0]];
    let loglikelihood = |_: &[f64; 0], parameters: &[f64; 1]| -> f64 {

        if parameters[0] < bounds[0][0] || parameters[0] > bounds[0][1] {
            std::f64::NEG_INFINITY
        } else {
            grf.get_pdf(1, Some(parameters[0]))
                .into_iter()
                .map(|x| x.ln())
                .sum()
        }
    };
    let most_likely = parameter_inference_uniform_prior::<_, 0, 1, 32, 1_000, 10_000>(
        &[[]],
        &bounds,
        &loglikelihood,
    );
    println!("{:?}", most_likely);


    // Plot likelihood
    let biases: Vec<f64> = (0..100)
        .map(|i| bounds[0][0] + i as f64 * (bounds[0][1]-bounds[0][0])/99.0)
        .collect();
    let lhs: Vec<f64> = biases
        .iter()
        .map(|b| loglikelihood(&[], &[*b]))
        .collect();
    println!("{lhs:?}");
    let three_fourths = *lhs.clone().select_nth_unstable_by(3*lhs.len()/4, |a, b| a.partial_cmp(&b).unwrap()).1;
    let (_min, max) = lhs.iter().minmax().into_option().unwrap();
    let root = SVGBackend::new("examples/bias_likelihood.svg", (1920, 1080))
        .into_drawing_area();
    root
        .fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        .caption("L(b|D)", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(
            bounds[0][0]..bounds[0][1],
            three_fourths..*max)
        .unwrap();

    chart
        .configure_mesh()
        .draw()
        .unwrap();

    chart
        .draw_series(LineSeries::new(
            biases.into_iter().zip(lhs),
            &BLACK,
        )).unwrap()
        .label("-log L(b|D)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw().unwrap();

    root.present().unwrap();


    // Plot most likely 1nn
    let cdf = grf.get_cdf(1, Some(most_likely[0]));
    let pcdf: Vec<(f64, f64)> = cdf
        .iter()
        .zip(&nns)
        .map(|(c,nn)| (*nn, c.min(1.0-c).clamp(f64::MIN_POSITIVE, 0.5)))
        .collect();
    let pcdf_measurements: Vec<(f64, f64)> = nns
        .iter()
        .enumerate()
        .map(|(c, nn)| {
            let c = c as f64 / nns.len().sub(1) as f64;
            (*nn, c.min(1.0-c).clamp(f64::MIN_POSITIVE, 0.5))
        })
        .collect();

    let root = SVGBackend::new("examples/most_likely_1nn.svg", (1920, 1080))
        .into_drawing_area();
    root
        .fill(&WHITE).unwrap();
    let (xmin, xmax) = nns.iter().minmax().into_option().unwrap();
    let mut chart = ChartBuilder::on(&root)
        .caption("Most likely", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(
            (*xmin..*xmax).log_scale(),
            (1e-3_f64..0.5).log_scale())
        .unwrap();

    chart
        .configure_mesh()
        .draw()
        .unwrap();

    
    chart
        .draw_series(LineSeries::new(
            pcdf,
            &BLACK,
        )).unwrap()
        .label("most likely")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));
    chart
        .draw_series(LineSeries::new(
            pcdf_measurements,
            &RED,
        )).unwrap()
        .label("measurements")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw().unwrap();

    root.present().unwrap();
}