
#![allow(unused)]


// fn legendre(

// ) -> {

// }

/// Only a few are implemented (the ones used in this library)
/// TODO: unit test
pub(crate) fn spherical_jn(
    n: u8,
    z: f64,
) -> f64 {
    match n {
        0 => z.sin()/z,
        2 => (3.0 / z.powi(2) - 1.0) * z.sin() / z - 3.0 * z.cos() / z.powi(2),
        4 => (5.0 * z * (-21.0 + 2.0 * z.powi(2)) * z.cos() + (105.0 - 45.0 * z.powi(2) + z.powi(4)) * z.sin()) / z.powi(5),
        _ => unimplemented!("only n=0,2,4 are implemented as those are the only relevant ones")
    }    
}

/// Only a few are implemented (the ones used in this library)
/// TODO: unit test
pub(crate) fn legendre_polynomial(
    n: u8,
    mu: f64,
) -> f64 {
    match n {
        0 => 1.0,
        2 => (3.0 * mu.powi(2) - 1.0) / 2.0,
        4 => (35.0 * mu.powi(4) - 30.0 * mu.powi(2) + 3.0) / 8.0,
        _ => unimplemented!("only n=0,2,4 are implemented as those are the only relevant ones")
    }    
}