#![doc = include_str!("../README.md")]
pub mod power;
pub mod scale_factor;

#[cfg(feature = "bayesian")]
pub mod bayesian;

pub mod correlation;

pub mod nn;

