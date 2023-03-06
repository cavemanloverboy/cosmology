use ::ndarray::{Array2, Axis};
use itertools::Itertools;
use ndarray_npy::NpzReader;
use plotly::{
    common::{Font, Title},
    *,
};
use plotters::prelude::*;
use std::ops::Sub;

use cosmology::{
    bayesian::parameter_inference_uniform_prior,
    correlation::{CorrelationFunction, CorrelationFunctionParameters},
    nn::grf::{GaussianRandomField, SpaceMode},
    power::{PowerSpectrum, TransferFunction},
};
use rand::{seq::SliceRandom, Rng, SeedableRng};

const REDSHIFT: f64 = 1.0;

#[cfg(feature = "use-mpi")]
use {
    balancer::Balancer, mpi::collective::CommunicatorCollectives, mpi::environment::Universe,
    mpi::traits::Communicator, std::sync::Arc,
};

#[cfg(feature = "use-mpi")]
lazy_static::lazy_static! {
    static ref UNIVERSE: Arc<Universe> = Arc::new(
        // mpi::initialize().unwrap()
        mpi::initialize_with_threading(mpi::Threading::Serialized).unwrap().0
    );
}

static GRF_COUNTER: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(1);

fn main() {
    // Get 1NN measurements from quijote.
    let (nns, nbar) = get_1nn_quijote_measurements();
    println!("Calculated NNs. Constructing GRF Theory");
    #[cfg(feature = "use-mpi")]
    let nns = {
        println!(
            "rank {} got nbar {nbar:.3e} and nns {:?}",
            UNIVERSE.world().rank(),
            &nns[..5]
        );
        let balancer = balancer::Balancer::<u8>::new(UNIVERSE.clone(), false);
        let our_subset = balancer.get_subset(&nns).to_vec();
        println!(
            "rank {} got nbar {nbar:.3e} and subset {:?}",
            UNIVERSE.world().rank(),
            &our_subset[..5]
        );
        our_subset
    };

    // Construct likelihood function, bounds on bias parameter's uniform prior
    let (bounds, loglikelihood) = get_bounds_and_ll(&nns, nbar);

    // Do parameter inference
    const BURN_IN: usize = 1;
    const WALKERS: usize = 8;
    const SAMPLES: usize = 4;
    const DATA_DIM: usize = 0; // the data is contained inside grf, owned by loglikelihood
    const PARAMETERS: usize = 2;

    #[cfg(not(feature = "use-mpi"))]
    {
        let (most_likely, chain) = {
            parameter_inference_uniform_prior::<_, DATA_DIM, PARAMETERS, WALKERS, BURN_IN, SAMPLES>(
                &[[]],
                &bounds,
                &loglikelihood,
            )
        };
        // Plot likelihood
        plot_likelihood(bounds, &loglikelihood);

        // Plot most likely 1nn
        plot_likely_1nn(most_likely, nbar, chain, &nns);
    }

    #[cfg(feature = "use-mpi")]
    {
        if UNIVERSE.world().rank() == 0 {
            let (most_likely, chain) = {
                parameter_inference_uniform_prior::<
                    _,
                    DATA_DIM,
                    PARAMETERS,
                    WALKERS,
                    BURN_IN,
                    SAMPLES,
                >(&[[]], &bounds, &loglikelihood)
            };
            // Plot likelihood
            // plot_likelihood(bounds, &loglikelihood);

            // Plot most likely 1nn
            plot_likely_1nn(most_likely, nbar, chain, &nns);
        } else {
            for _ in 0..WALKERS * (BURN_IN + SAMPLES) {
                loglikelihood(&[], &[2.5, 0.3]);
            }
        }
        UNIVERSE.world().barrier();
    }
}

fn get_1nn_quijote_measurements() -> (Vec<f64>, f64) {
    const NDATA: usize = 1_000;
    const QUIJOTE_BOXSIZE: [f64; 3] = [1000.0; 3];
    const N_QUERY: usize = 10;
    let mut npz = NpzReader::new(std::fs::File::open("./OneBox/OneBox.npz").unwrap()).unwrap();
    let real_data: Array2<f32> = npz.by_name("rpos.npy").unwrap();
    let real_data: Array2<f64> = real_data.map(|x| *x as f64);

    // subsample
    let mut real_data: Vec<[f64; 3]> = real_data
        .map_axis(Axis(1), |x| x.as_slice().unwrap().try_into().unwrap())
        .into_raw_vec();
    let mut seeded_rng = rand_chacha::ChaCha8Rng::from_seed([1; 32]);
    (&mut real_data).shuffle(&mut seeded_rng);
    real_data.truncate(NDATA);
    let nbar = NDATA as f64 / QUIJOTE_BOXSIZE[0].powi(3);

    // Find 1nn
    println!("Calculating NNs");
    let tree = fnntw::Tree::<'_, 3>::new(&real_data, 32)
        .unwrap()
        .with_boxsize(&QUIJOTE_BOXSIZE)
        .unwrap();
    let mut nns: Vec<f64> = Vec::with_capacity(2 * N_QUERY);
    for _ in 0..2 * N_QUERY {
        let query: [f64; 3] = seeded_rng.gen::<[f64; 3]>().map(|x| x * QUIJOTE_BOXSIZE[0]);
        nns.push(tree.query_nearest(&query).unwrap().0.sqrt());
    }
    nns.sort_by(|a, b| a.partial_cmp(&b).expect("there should be no NaNs"));
    nns.drain(..N_QUERY);
    (nns, nbar)
}

fn get_correlation_function<'c>(omega_matter_0: f64) -> Box<dyn Fn(f64) -> f64 + Send + Sync + 'c> {
    // Initialize E & Hu power spectrum
    // Quijote Cosmology
    // let omega_matter_0 = 0.3175;
    let omega_baryon_0 = 0.049;
    let h = 0.6711;
    let ns = 0.9624;
    let sigma_8 = 0.834;
    let power = PowerSpectrum::new(TransferFunction::EisensteinHu {
        h,
        omega_matter_0,
        omega_baryon_0,
        temp_cmb0: 2.7255,
        ns,
        sigma_8,
    })
    .unwrap();
    let params = CorrelationFunctionParameters {
        power,
        accuracy_params: None, // default accuracy parameters
    };
    let real_corr = CorrelationFunction::get_correlation_function(REDSHIFT, params).unwrap();
    let real_corr_fn = move |r| real_corr.correlation_function(r);
    Box::new(real_corr_fn)
}
fn construct_grf<'b, 'c>(
    real_corr_fn: &'b Box<dyn Fn(f64) -> f64 + Send + Sync + 'c>,
    nns: &Vec<f64>,
    nbar: f64,
) -> GaussianRandomField<'b, 'c>
where
    'b: 'c,
{
    let mode = SpaceMode::RealSpace(&*real_corr_fn);
    let grf = GaussianRandomField::new(mode).with(&nns);

    grf
}

fn get_bounds_and_ll<'b, 'c>(
    nns: &'b Vec<f64>,
    nbar: f64,
) -> (
    [[f64; 2]; 2],
    Box<dyn Fn(&[f64; 0], &[f64; 2]) -> f64 + Send + Sync + 'c>,
)
where
    'b: 'c,
{
    //let bounds = [[0.8, 10.0], [0.05, 0.999]];
    let bounds = [[1.75, 3.75], [0.05, 0.999]];
    let loglikelihood = move |_: &[f64; 0], parameters: &[f64; 2]| -> f64 {
        #[cfg(feature = "use-mpi")]
        let balancer: Balancer<f64> = Balancer::new(UNIVERSE.clone(), false);
        #[allow(unused_mut)] // if using mpi
        let mut bias = parameters[0];
        #[allow(unused_mut)] // if using mpi
        let mut omega_matter_0 = parameters[1];

        // Syncrhonize across all things, prior to checking for bounds!
        #[cfg(feature = "use-mpi")]
        {
            balancer.synchronize_value(&mut omega_matter_0);
            balancer.synchronize_value(&mut bias);
        }

        if bias < bounds[0][0]
            || bias > bounds[0][1]
            || omega_matter_0 < bounds[1][0]
            || omega_matter_0 > bounds[1][1]
        {
            std::f64::NEG_INFINITY
        } else {
            let start = std::time::Instant::now();

            // Construct Gaussian Random Field Theory at those scales.
            let real_corr = get_correlation_function(omega_matter_0);
            let grf = construct_grf(&real_corr, &nns, nbar);

            let local_result: f64 = grf
                .get_pdf(1, nbar, Some(bias))
                .into_iter()
                .map(|x| x.ln())
                .sum();
            #[cfg(not(feature = "use-mpi"))]
            return local_result;

            #[cfg(feature = "use-mpi")]
            {
                let mut global_result = 0.0;
                balancer.world().all_reduce_into(
                    &local_result,
                    &mut global_result,
                    mpi::collective::SystemOperation::sum(),
                );
                // println!(
                //     "rank {}: [{bias:.3}, {omega_matter_0:.3}]: local {local_result}, global {global_result}",
                //     balancer.rank
                // );
                // if balancer.rank == 0 {
                println!(
                        "finished grf #{} (in {:.2} secs), ll({bias:.2}, {omega_matter_0:.2}) = {global_result:.3e}",
                        GRF_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed), start.elapsed().as_secs_f32()
                    );
                // }
                global_result
            }
        }
    };
    (bounds, Box::new(loglikelihood))
}

fn plot_likelihood(bounds: [[f64; 2]; 2], loglikelihood: &dyn Fn(&[f64; 0], &[f64; 2]) -> f64) {
    let biases: Vec<f64> = (0..100)
        .map(|i| bounds[0][0] + i as f64 * (bounds[0][1] - bounds[0][0]) / 99.0)
        .collect();
    let omegas: Vec<f64> = (0..100)
        .map(|i| bounds[1][0] + i as f64 * (bounds[1][1] - bounds[0][0]) / 99.0)
        .collect();
    let lhs: Vec<f64> = biases
        .iter()
        .zip(&omegas)
        .map(|(b, om)| loglikelihood(&[], &[*b, *om]))
        .collect();
    println!("{lhs:?}");
    let three_fourths = *lhs
        .clone()
        .select_nth_unstable_by(3 * lhs.len() / 4, |a, b| a.partial_cmp(&b).unwrap())
        .1;
    let (_min, max) = lhs.iter().minmax().into_option().unwrap();
    let root = SVGBackend::new("examples/bias_likelihood.svg", (1920, 1080)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let mut chart = ChartBuilder::on(&root)
        .caption("L(b|D)", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(bounds[0][0]..bounds[0][1], three_fourths..*max)
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(LineSeries::new(biases.into_iter().zip(lhs), &BLACK))
        .unwrap()
        .label("-log L(b|D)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();

    root.present().unwrap();
}

fn plot_likely_1nn(most_likely: [f64; 2], nbar: f64, mut chain: Vec<[f64; 2]>, nns: &Vec<f64>) {
    let mut chain_file = std::fs::File::create("bias_omega_chain_1e4").unwrap();
    use std::io::Write;
    for pair in &chain {
        chain_file
            .write(format!("{} {}\n", pair[0], pair[1]).as_bytes())
            .unwrap();
    }
    let [most_likely_b, most_likely_omega_matter_0] = most_likely;
    let most_likely_real_corr_fn = get_correlation_function(most_likely_omega_matter_0);
    let most_likely_grf = construct_grf(&most_likely_real_corr_fn, nns, nbar);
    let cdf = most_likely_grf.get_cdf(1, nbar, Some(most_likely_b));
    let pcdf: Vec<(f64, f64)> = cdf
        .iter()
        .zip(nns)
        .map(|(c, nn)| (*nn, c.min(1.0 - c).clamp(f64::MIN_POSITIVE, 0.5)))
        .collect();
    let pcdf_measurements: Vec<(f64, f64)> = nns
        .iter()
        .enumerate()
        .map(|(c, nn)| {
            let c = c as f64 / nns.len().sub(1) as f64;
            (*nn, c.min(1.0 - c).clamp(f64::MIN_POSITIVE, 0.5))
        })
        .collect();

    let root = SVGBackend::new("examples/most_likely_1nn2.svg", (1920, 1080)).into_drawing_area();
    root.fill(&WHITE).unwrap();
    let (xmin, xmax) = nns.iter().minmax().into_option().unwrap();
    let mut chart = ChartBuilder::on(&root)
        .caption("Most likely", ("sans-serif", 50).into_font())
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d((*xmin..*xmax).log_scale(), (1e-3_f64..0.5).log_scale())
        .unwrap();

    chart.configure_mesh().draw().unwrap();

    chart
        .draw_series(LineSeries::new(pcdf, &BLACK))
        .unwrap()
        .label("most likely")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLACK));
    chart
        .draw_series(LineSeries::new(pcdf_measurements, &RED))
        .unwrap()
        .label("measurements")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()
        .unwrap();

    root.present().unwrap();

    // Extract bounds of credible interval
    let left_2p5 = chain.len() * 25 / 1000;
    let right_2p5 = chain.len() * 975 / 1000;
    chain.sort_unstable_by(|a, b| a[0].partial_cmp(&b[0]).unwrap()); // this only works for 1-parameter models. need to sort each parameter independently.
    let left_2p5_bias = chain[left_2p5][0];
    let right_2p5_bias = chain[right_2p5][0];
    println!("{right_2p5_bias}");
    let leftmost_bias = chain[0][0];
    let rightmost_bias = chain.last().unwrap()[0];
    const BINS: usize = 100;
    let bin_size_bias = (rightmost_bias - leftmost_bias) / BINS as f64;

    // Now get omega_matter
    chain.sort_unstable_by(|a, b| a[1].partial_cmp(&b[1]).unwrap()); // this only works for 1-parameter models. need to sort each parameter independently.
    let left_2p5_omega = chain[left_2p5][1];
    let right_2p5_omega = chain[right_2p5][1];
    println!("{right_2p5_omega}");
    let leftmost_omega = chain[0][1];
    let rightmost_omega = chain.last().unwrap()[1];
    let bin_size_omega = (rightmost_omega - leftmost_omega) / BINS as f64;

    // // manually calculate hist...
    // let bin_num = |x: &[f64; 1]| ((x[0] - 0.8) / (10.0 - 0.8) * 100.0) as usize;
    // let y = chain.iter().counts_by(bin_num).into_iter().map(|(bin_num, count)| );

    let hist = plotly::Histogram::new(chain.iter().map(|x| x[0]).collect::<Vec<f64>>())
        .x_bins(plotly::histogram::Bins::new(
            leftmost_bias,
            rightmost_bias,
            bin_size_bias,
        ))
        .hist_norm(plotly::histogram::HistNorm::Probability)
        .opacity(0.5);
    println!("hist: {hist:?}");
    let hist2 = plotly::Histogram::new(chain.iter().map(|x| x[1]).collect::<Vec<f64>>())
        .x_bins(plotly::histogram::Bins::new(
            leftmost_omega,
            rightmost_omega,
            bin_size_omega,
        ))
        .hist_norm(plotly::histogram::HistNorm::Probability)
        .opacity(0.5);
    println!("hist2: {hist2:?}");

    const MARGIN: usize = 120;
    const FONT_SIZE: usize = 40;
    let mut layout = plotly::Layout::new()
        .x_axis(
            plotly::layout::Axis::new()
                .title(
                    Title::from("Halo\tBias")
                        .font(Font::new().size(FONT_SIZE))
                        .x(0.0),
                )
                .grid_color(NamedColor::DarkGray),
        )
        .y_axis(
            plotly::layout::Axis::new()
                .title(
                    Title::from("Probability")
                        .font(Font::new().size(FONT_SIZE))
                        .x(0.0),
                )
                .grid_color(NamedColor::DarkGray),
        )
        .margin(
            plotly::layout::Margin::new()
                .top(MARGIN / 2)
                .bottom(MARGIN + 5)
                .right(MARGIN - 5)
                .left(MARGIN),
        )
        .font(plotly::common::Font::new().size(30));
    layout.add_shape(
        plotly::layout::Shape::new()
            .shape_type(plotly::layout::ShapeType::Line)
            .y_ref("paper")
            .x0(left_2p5_bias)
            .y0(0)
            .x1(left_2p5_bias)
            .y1(1)
            .line(
                plotly::layout::ShapeLine::new()
                    .color(plotly::NamedColor::Black)
                    .width(3.),
            ),
    );
    layout.add_shape(
        plotly::layout::Shape::new()
            .shape_type(plotly::layout::ShapeType::Line)
            .y_ref("paper")
            .x0(right_2p5_bias)
            .y0(0)
            .x1(right_2p5_bias)
            .y1(1)
            .line(
                plotly::layout::ShapeLine::new()
                    .color(plotly::NamedColor::Black)
                    .width(3.),
            ),
    );
    layout.add_shape(
        plotly::layout::Shape::new()
            .shape_type(plotly::layout::ShapeType::Rect)
            .x_ref("paper")
            .y_ref("paper")
            .x0(0)
            .x1(1)
            .y0(0)
            .y1(1),
    );

    let mut plot = plotly::Plot::new();
    plot.set_layout(layout);
    plot.add_trace(hist);
    plot.save(
        "examples/bias_credible_interval2_2.png",
        plotly::ImageFormat::PNG,
        1920,
        1080,
        1.0,
    );

    let mut layout = plotly::Layout::new()
        .x_axis(
            plotly::layout::Axis::new()
                .title(
                    Title::from("Omega")
                        .font(Font::new().size(FONT_SIZE))
                        .x(0.0),
                )
                .grid_color(NamedColor::DarkGray),
        )
        .y_axis(
            plotly::layout::Axis::new()
                .title(
                    Title::from("Probability")
                        .font(Font::new().size(FONT_SIZE))
                        .x(0.0),
                )
                .grid_color(NamedColor::DarkGray),
        )
        .margin(
            plotly::layout::Margin::new()
                .top(MARGIN / 2)
                .bottom(MARGIN + 5)
                .right(MARGIN - 5)
                .left(MARGIN),
        )
        .font(plotly::common::Font::new().size(30));
    layout.add_shape(
        plotly::layout::Shape::new()
            .shape_type(plotly::layout::ShapeType::Line)
            .y_ref("paper")
            .x0(left_2p5_omega)
            .y0(0)
            .x1(left_2p5_omega)
            .y1(1)
            .line(
                plotly::layout::ShapeLine::new()
                    .color(plotly::NamedColor::Black)
                    .width(3.),
            ),
    );
    layout.add_shape(
        plotly::layout::Shape::new()
            .shape_type(plotly::layout::ShapeType::Line)
            .y_ref("paper")
            .x0(right_2p5_omega)
            .y0(0)
            .x1(right_2p5_omega)
            .y1(1)
            .line(
                plotly::layout::ShapeLine::new()
                    .color(plotly::NamedColor::Black)
                    .width(3.),
            ),
    );
    layout.add_shape(
        plotly::layout::Shape::new()
            .shape_type(plotly::layout::ShapeType::Rect)
            .x_ref("paper")
            .y_ref("paper")
            .x0(0)
            .x1(1)
            .y0(0)
            .y1(1),
    );

    let mut plot2 = plotly::Plot::new();
    plot2.set_layout(layout);
    plot2.add_trace(hist2);
    plot2.save(
        "examples/omega_credible_interval_2.png",
        plotly::ImageFormat::PNG,
        1920,
        1080,
        1.0,
    );
}
