


use std::f64::consts::PI;

use assert_approx_eq::assert_approx_eq;
use cosmology::correlation::{CorrelationFunctionParameters, CorrelationFunction};
use criterion::{Criterion, criterion_group, criterion_main, black_box};

const LOGR_MIN: f64 = 0.0;
const LOGR_MAX: f64 = 2.0;
const NUM_R_POINTS: usize = 200;

fn bench_1d(c: &mut Criterion) {

    let mut group = c.benchmark_group("1D");

    // Integrand for half unit circle
    let half_circle = |x: f64| (1.0 - x.powi(2)).sqrt();

    group.bench_function("quadrature/clenshaw-curtis", move |b| {
        b.iter(|| {
            let a = black_box(quadrature::clenshaw_curtis::integrate(
                black_box(half_circle),
                black_box(-1.0),
                black_box(1.0),
                black_box(1e-6)
            ));
            assert_approx_eq!(a.integral, PI / 2.0, a.error_estimate);
        })
    });

    group.bench_function("quadrature/double-exponential", move |b| {
        b.iter(|| {
            let a = black_box(
                quadrature::double_exponential::integrate(
                    black_box(half_circle),
                    black_box(-1.0),
                    black_box(1.0),
                    black_box(1e-6)
                )
            );
            assert_approx_eq!(a.integral, PI / 2.0, a.error_estimate);
        })
    });

    group.bench_function("sequential", move |b| {
        b.iter(|| {
            let a = black_box(
            sequential_integration::calculate_single_integral_simpson(
                    black_box(half_circle),
                    black_box(-1.0),
                    black_box(1.0),
                    1e-6
                ).unwrap()
            );
            assert_approx_eq!(a, PI / 2.0, 1e-6);
        })
    });
    
}

fn bench_2d(c: &mut Criterion) {

    let mut group = c.benchmark_group("2D");

    // Integrand for half unit sphere
    let half_sphere = |x: f64, y: f64| (1.0 - x.powi(2) - y.powi(2)).sqrt();

    group.bench_function("quadrature/clenshaw-curtis (nested)", move |b| {
        b.iter(|| {
            let a = black_box(quadrature::clenshaw_curtis::integrate(
                |y: f64| {
                    quadrature::clenshaw_curtis::integrate(
                        |x: f64| half_sphere(black_box(x), black_box(y)),
                        black_box(-(1.0 - y.powi(2)).sqrt()),
                        black_box((1.0 - y.powi(2)).sqrt()),
                        1e-6,
                    ).integral
                },
                -1.0,
                1.0,
                1e-6
            ));
            // We are a little more lenient with this error because we are discarding
            // the error of the internal integral.
            assert_approx_eq!(a.integral, 2.0 * PI / 3.0, 1e-6);
        })
    });

    group.bench_function("quadrature/double-exponential (nested)", move |b| {
        b.iter(|| {
            let a = black_box(quadrature::double_exponential::integrate(
                |y: f64| {
                    quadrature::double_exponential::integrate(
                        |x: f64| half_sphere(black_box(x), black_box(y)),
                        black_box(-(1.0 - y.powi(2)).sqrt()),
                        black_box((1.0 - y.powi(2)).sqrt()),
                        1e-6,
                    ).integral
                },
                -1.0,
                1.0,
                1e-6
            ));
            // We are a little more lenient with this error because we are discarding
            // the error of the internal integral.
            assert_approx_eq!(a.integral, 2.0 * PI / 3.0, 1e-6);
        })
    });



    // group.bench_function("sequential", move |b| {
    //     b.iter(|| {
    //         let a = black_box(
    //             sequential_integration::calculate_double_integral_simpson(
    //                 half_sphere,
    //                 black_box(-1.0),
    //                 black_box(1.0),
    //                 1e-3,
    //                 |x: f64| -((1.0 - x.powi(2)).max(0.0)).sqrt(),
    //                 |x: f64|  ((1.0 - x.powi(2)).max(0.0)).sqrt(),
    //                 1e-3,
    //             ).unwrap()
    //         );
    //         assert_approx_eq!(a, 2.0 * PI / 3.0, 1e-6);
    //     })
    // });
}


criterion_group!(benches, bench_1d, bench_2d);
criterion_main!(benches);