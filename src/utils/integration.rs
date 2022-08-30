#![allow(unused)]

use quadrature::Output;


pub(crate) fn two_dimensional_integral(
    integrand: &dyn Fn(f64, f64) -> f64,
    outer_lower: f64,
    outer_upper: f64,
    inner_lower: &dyn Fn(f64) -> f64,
    inner_upper: &dyn Fn(f64) -> f64,
) -> f64 {

    let outer_result = quadrature::clenshaw_curtis::integrate(
        |y: f64| {
            let inner_result = quadrature::clenshaw_curtis::integrate(
                |x: f64| integrand(x, y),
                inner_lower(y),
                inner_upper(y),
                1e-6,
            );
            check_integral(&inner_result);
            inner_result.integral
        },
        outer_lower,
        outer_upper,
        1e-6
    );
    check_integral(&outer_result);
    outer_result.integral
}

pub(crate) fn two_dimensional_integral_const_inner_bounds(
    integrand: &dyn Fn(f64, f64) -> f64,
    outer_lower: f64,
    outer_upper: f64,
    inner_lower: f64,
    inner_upper: f64,
    target_error: f64,
) -> f64 {

    let outer_result = quadrature::clenshaw_curtis::integrate(
        |r: f64| {
            let inner_result = quadrature::clenshaw_curtis::integrate(
                |t: f64| integrand(r, t),
                inner_lower,
                inner_upper,
                target_error,
            );
            check_integral(&inner_result);
            inner_result.integral
        },
        outer_lower,
        outer_upper,
        target_error
    );
    check_integral(&outer_result);
    outer_result.integral
}

fn check_integral(i: &Output) {
    // TODO
}