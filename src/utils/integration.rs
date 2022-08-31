#![allow(unused)]

use quadrature::Output;

use crate::scale_factor::rk4;


pub(crate) fn rk4_integrator(
    integrand: &dyn Fn(f64) -> f64,
    mut lower_bound: f64,
    upper_bound: f64,
    mut step: f64,
) -> f64 {

    let mut result = 0.0;
    let rk4_integrand = |x: f64, _: f64| integrand(x);
    while lower_bound < upper_bound {
        step = step.min(upper_bound-lower_bound);
        result = rk4(rk4_integrand, lower_bound, result, step, None);
        lower_bound += step;
    }
    result
}

// pub(crate) fn rk4_integrator_adaptive(
//     integrand: &dyn Fn(f64) -> f64,
//     mut lower_bound: f64,
//     upper_bound: f64,
//     step_hint: f64,
//     adaptive_tol: f64,
// ) -> f64 {

//     let mut result = 0.0;
//     let rk4_integrand = |x: f64, _: f64| integrand(x);
//     let mut current_step; = step_hint;
//     while lower_bound < upper_bound {

//         // 
//         current_step = step_hint.max(current_step).min(upper_bound-lower_bound);
//         let mut temp = rk4(rk4_integrand, lower_bound, result, current_step, None);
//         let temp2 = 
//         while (1.0 - temp/rk4(rk4_integrand, lower_bound, result, step/2, None)).abs()  {
//             step
//         }
//         result = rk4(rk4_integrand, lower_bound, result, step, None);
//         lower_bound += step;
//     }
//     result
// }

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