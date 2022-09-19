use std::ops::Sub;

use cosmology::{
    correlation::{CorrelationFunction, CorrelationFunctionParameters},
    nn::{
        grf::{GaussianRandomField, SpaceMode},
        poisson::PoissonRandomField,
    },
    power::{PowerSpectrum, TransferFunction},
};

use plotly::{common::{Mode, Title, Font, color::Color, Anchor, Line, DashType}, Layout, layout::{Axis, Margin, AxisType, Shape, ShapeType, Legend}, NamedColor};
use plotly::{Plot, Scatter, ImageFormat};


const REDSHIFT: f64 = 1.0;
const MIN_R: f64 = 4.0;
const MAX_R: f64 = 40.0;
const NUM_R_POINTS: usize = 1000;
const PLOT_PADDING: usize = 60*2;
const FONT_SIZE: usize = 20*2;

fn main() {
    // 100,000 objects per (1000 Mpc/h)^3
    let nbar = 1e5 / 1e3_f64.powi(3);

    // Scales of interest
    let scales: Vec<f64> = (0..NUM_R_POINTS)
        .map(|i| {
            10_f64.powf(MIN_R.log10() + i as f64 * (MAX_R.log10() - MIN_R.log10()) / NUM_R_POINTS.sub(1) as f64)
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
    let cdf = grf.get_cdf(1, Some(2.0));
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

    // then redshift space
    let rsd = SpaceMode::RedshiftSpace {
        real_corr: &real_corr_fn,
        f: omega_matter_0.powf(5.0 / 9.0),
    };
    let rsd_grf = GaussianRandomField::new(nbar, rsd).with(&scales);
    let rsd_cdf = rsd_grf.get_cdf(1, Some(2.0));
    let rsd_pcdf: Vec<f64> = rsd_cdf
        .iter()
        .map(|c| c.min(1.0 - c).clamp(f64::MIN_POSITIVE, 0.5))
        .collect();
    let unbiased_rsd = SpaceMode::RedshiftSpace {
        real_corr: &real_corr_fn,
        f: omega_matter_0.powf(5.0 / 9.0),
    };
    let unbiased_rsd_grf = GaussianRandomField::new(nbar, unbiased_rsd).with(&scales);
    let unbiased_rsd_cdf = unbiased_rsd_grf.get_cdf(1, None);
    let unbiased_rsd_pcdf: Vec<f64> = unbiased_rsd_cdf
        .iter()
        .map(|c| c.min(1.0 - c).clamp(f64::MIN_POSITIVE, 0.5))
        .collect();

    // Compare to poisson
    let poisson = PoissonRandomField::new(nbar).with(&scales);
    let poisson_cdf = poisson.get_cdf(1);
    let poisson_pcdf: Vec<f64> = poisson_cdf
        .iter()
        .map(|c| c.min(1.0 - c).clamp(f64::MIN_POSITIVE, 0.5))
        .collect();

    let poisson_line = Scatter::new(scales.clone(), poisson_pcdf)
        .name("P")
        .line(Line::new()
            .width(5.0))
        .mode(Mode::Lines);
    let g2_line = Scatter::new(scales.clone(), pcdf)
        .name("G(2)")
        .line(Line::new()
            .width(5.0))
        .mode(Mode::Lines);
    let g1_line = Scatter::new(scales.clone(), unbiased_pcdf)
        .name("G(1)")
        .line(Line::new()
            .width(5.0))
        .mode(Mode::Lines);
    let gr2_line = Scatter::new(scales.clone(), rsd_pcdf)
        .name("GR(2)")
        .line(Line::new()
            .width(5.0))
        .mode(Mode::Lines);
    let gr1_line = Scatter::new(scales.clone(), unbiased_rsd_pcdf)
        .name("GR(1)")
        .line(Line::new()
            .width(5.0))
        .mode(Mode::Lines);

    let mut plot = Plot::new();
    plot.add_trace(poisson_line);
    plot.add_trace(g1_line);
    plot.add_trace(gr1_line);
    plot.add_trace(g2_line);
    plot.add_trace(gr2_line);

    let x0 = 10_f64.powf(MIN_R.log10());
    let x1 = 10_f64.powf(MAX_R.log10());
    let y0 = 10_f64.powf(-3.0);
    let y1 = 10_f64.powf(0.0);
    let mut layout = Layout::new()
        .x_axis(Axis::new()
            .title(Title::from("Distance\t[Mpc/h]").font(Font::new().size(FONT_SIZE)).x(0.0))
            .range(vec![x0.log10(), x1.log10()])
            .type_(AxisType::Log)
            .tick_values(vec![1.0, 5.0, 10.0, 20.0, 30.0, 40.0].into_iter().filter(|&x| x > MIN_R && x < MAX_R).collect())
            .grid_color(NamedColor::DarkGray))
        .y_axis(Axis::new()
            .title(Title::from("Peaked\tCDF").font(Font::new().size(FONT_SIZE)).x(0.0))
            .range(vec![y0.log10(), y1.log10()])
            .type_(AxisType::Log)
            .tick_values(vec![1e-3, 1e-2, 1e-1, 1e0])
            .grid_color(NamedColor::DarkGray))
        .margin(Margin::new()
            .top(PLOT_PADDING/2)
            .bottom(PLOT_PADDING+5)
            .right(PLOT_PADDING-5)
            .left(PLOT_PADDING))
        .font(Font::new().size(FONT_SIZE-4))
        .legend(Legend::new()
            .x(1.0)
            .y(1.0)
            .x_anchor(Anchor::Right)
            .border_width(1));
    layout
        .add_shape(Shape::new()
            .shape_type(ShapeType::Rect)
            .x0(x0)
            .x1(x1)
            .y0(y0)
            .y1(y1));
    plot.set_layout(layout);


    // The following will save the plot in all available formats and show the plot.
    plot.save("examples/Figure2.png", ImageFormat::PNG, 1920 * 5 / 6, 1080, 1.0);
    
}
