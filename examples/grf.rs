use plotters::prelude::*;
use std::ops::Sub;

use cosmology::{
    correlation::{CorrelationFunction, CorrelationFunctionParameters},
    nn::{
        grf::{GaussianRandomField, SpaceMode},
        poisson::PoissonRandomField,
    },
    power::{PowerSpectrum, TransferFunction},
};

const REDSHIFT: f64 = 1.0;
const MIN_LOGR: f64 = 0.0;
const MAX_LOGR: f64 = 2.0;
const NUM_R_POINTS: usize = 200;

fn main() {
    // 100,000 objects per (1000 Mpc/h)^3
    let nbar = 1e4 / 1e3_f64.powi(3);

    // Scales of interest
    let scales: Vec<f64> = (0..NUM_R_POINTS)
        .map(|i| {
            10_f64.powf(MIN_LOGR + i as f64 * (MAX_LOGR - MIN_LOGR) / (NUM_R_POINTS.sub(1) as f64))
        })
        .collect();

    // Initialize E & Hu power spectrum
    let omega_matter_0 = 0.3;
    let power = PowerSpectrum::new(TransferFunction::EisensteinHu {
        h: 0.7,
        omega_matter_0,
        omega_baryon_0: 0.03,
        temp_cmb0: 2.7,
        ns: 0.9665,
        sigma_8: 0.681,
    })
    .unwrap();

    // Do real space first
    let params = CorrelationFunctionParameters {
        power,
        accuracy_params: None, // default accuracy parameters
    };
    let real_corr = CorrelationFunction::get_correlation_function(REDSHIFT, params).unwrap();
    let real_corr_fn = |r| real_corr.correlation_function(r);
    let mode = SpaceMode::RealSpace(&real_corr_fn);
    let grf = GaussianRandomField::new(nbar, mode)
        .with(&scales);
    let cdf = grf.get_cdf(1, Some(3.0));
    let pcdf: Vec<f64> = cdf
        .iter()
        .map(|c| c.min(1.0 - c).clamp(f64::MIN_POSITIVE, 0.5))
        .collect();

    // then unbiased
    let unbiased = SpaceMode::RealSpace(&real_corr_fn);
    let unbiased_grf = GaussianRandomField::new(nbar, unbiased)
        .with(&scales);
    let unbiased_cdf = unbiased_grf.get_cdf(1, None);
    let unbiased_pcdf: Vec<f64> = unbiased_cdf
        .iter()
        .map(|c| c.min(1.0 - c).clamp(f64::MIN_POSITIVE, 0.5))
        .collect();

    // // then redshift space
    // let rsd = SpaceMode::RedshiftSpace {
    //     real_corr: &real_corr_fn,
    //     f: omega_matter_0.powf(5.0 / 9.0),
    // };
    // let rsd_grf = GaussianRandomField::new(nbar, rsd).with(&scales);
    // let rsd_cdf = rsd_grf.get_cdf(1, Some(3.0));
    // let rsd_pcdf: Vec<f64> = rsd_cdf
    //     .iter()
    //     .map(|c| c.min(1.0 - c).clamp(f64::MIN_POSITIVE, 0.5))
    //     .collect();

    // Compare to poisson
    let poisson = PoissonRandomField::new(nbar).with(&scales);
    let poisson_cdf = poisson.get_cdf(1);
    let poisson_pcdf: Vec<f64> = poisson_cdf
        .iter()
        .map(|c| c.min(1.0 - c).clamp(f64::MIN_POSITIVE, 0.5))
        .collect();

    let root = BitMapBackend::new("examples/grf_vs_poisson.png", (1920 / 2, 1080)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        .caption("Gaussian vs Poisson", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(60)
        .y_label_area_size(60)
        .build_cartesian_2d(
            (scales[0]..*scales.last().unwrap()).log_scale(),
            (1e-3..1.0_f64).log_scale(),
        )
        .unwrap();

    chart
        .configure_mesh()
        .y_desc("1NN-PCDF")
        .x_desc("Distance [Mpc/h]")
        .axis_desc_style(("Cambria", 40))
        .draw()
        .unwrap();

    // Add GRF
    chart
        .draw_series(LineSeries::new(
            scales.iter().cloned().zip(pcdf.clone()),
            BLACK.stroke_width(3),
        ))
        .unwrap()
        .label("GRF (b=3)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));

    // Add unbiased
    chart
        .draw_series(LineSeries::new(
            scales.iter().cloned().zip(unbiased_pcdf.clone()),
            RED.stroke_width(3),
        ))
        .unwrap()
        .label("GRF (b=1)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Add GRF-RSD
    // chart
    //     .draw_series(LineSeries::new(
    //         scales.iter().cloned().zip(rsd_pcdf.clone()),
    //         &RED,
    //     ))
    //     .unwrap()
    //     .label("GRF-RSD")
    //     .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(LineSeries::new(
            scales.iter().cloned().zip(poisson_pcdf.clone()),
            BLUE.stroke_width(3),
        ))
        .unwrap()
        .label("Poisson")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .configure_series_labels()
        .label_font(("Calibri", 40))
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();

    root.present().unwrap();
}
