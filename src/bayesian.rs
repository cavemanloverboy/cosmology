use hammer_and_sample::{sample, Model, Parallel};
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;
use std::sync::Arc;

struct Data<'a, 'b, const D: usize, const P: usize> {
    data: &'a [[f64; D]],
    // bounds: Option<(f64, f64)>,
    loglikelihood: Arc<dyn Fn(&[f64; D], &[f64; P]) -> f64 + Send + Sync + 'b>,
}

impl<const D: usize, const K: usize> Model for Data<'_, '_, D, K> {
    type Params = [f64; K];
    fn log_prob(&self, params: &Self::Params) -> f64 {
        self.data.iter().map(|d| (self.loglikelihood)(d, params)).sum()
    }
}

pub fn parameter_inference_uniform_prior<
    'a,
    'b,
    L,
    const D: usize,
    const P: usize,
    const W: usize,
    const B: usize,
    const S: usize,
>(
    data: &'a [[f64; D]],
    // TODO: bounds check
    bounds: &[[f64; 2]; P],
    loglikelihood: &L,
) -> [f64; P]
where
    L: Fn(&[f64; D], &[f64; P]) -> f64 + Send + Sync + 'b,
{
    // Check bounds
    bounds.iter().enumerate().for_each(|(i, [low, high])| {
        assert!(
            high > low,
            "Invalid bounds for parameter {i} is invalid. low = {low}, high = {high}"
        )
    });

    // Initialize helper struct
    let data = Data {
        data,
        loglikelihood: Arc::new(loglikelihood),
    };

    // Initialize walkers iter
    let walkers = (0..W as u64).map(|seed| {
        let mut rng = Pcg64::seed_from_u64(seed);

        let p = (0..P)
            .map(|i| rng.gen_range(bounds[i][0]..=bounds[i][1]))
            .collect::<Vec<_>>()
            .try_into()
            .expect("shape should be as expected");

        (p, rng)
    });

    let (chain, _accepted) = sample(&data, walkers, S, &Parallel);

    // Remove the first burn_in samples from each walker
    println!("chain[..15] = {:?}", &chain[..15]);
    println!("chain[-15..] = {:?}", &chain.iter().rev().take(15).rev().collect::<Vec<_>>());
    let chain = &chain[W * B..];

    // TODO: visualize chain

    let norm = chain.len() as f64;
    chain
        .iter()
        .fold([0.0; P], |mut acc, p| {
            acc.iter_mut().enumerate().for_each(|(i, a)| *a += p[i]);
            acc
        })
        .map(|p| p / norm)
}

#[test]
fn test_coin_flip() {
    // Coin flip example, but with floats instead of bools

    // First define bernoulli
    const P: usize = 1;
    const D: usize = 1;
    let loglikelihood = |x: &[f64; D], p: &[f64; P]| {
        if x[0] == 1.0 {
            p[0].ln()
        } else if x[0] == 0.0 {
            (1.0 - p[0]).ln()
        } else {
            unreachable!()
        }
    };

    // Then define some data
    let data = [vec![[1.0]; 1000], vec![[0.0]; 500]].concat();

    // Uniform prior (also bounds of parameters)
    let bounds = [[0.0, 1.0]];

    const WALKERS: usize = 32;
    const BURN_IN: usize = 100;
    const SAMPLES: usize = 1000;
    let params = parameter_inference_uniform_prior::<_, D, P, WALKERS, BURN_IN, SAMPLES>(
        &data, &bounds, &loglikelihood,
    );

    // Check that param is close to what we expected
    assert!(params[0] > 0.65 && params[0] < 0.67);
}

#[test]
fn test_gaussian() {
    use rand_distr::{Distribution, Normal};
    use std::f64::consts::PI;

    // First define gaussian
    const P: usize = 2;
    const D: usize = 1;
    let loglikelihood = |x: &[f64; 1], p: &[f64; P]| {
        if p[1] < 0.0 {
            return std::f64::NEG_INFINITY;
        } else {
            -0.5 * (((x[0] - p[0]) / p[1]).powf(2.0) + (2.0 * PI * p[1].powi(2)).ln())
        }
    };

    // Then define some data
    let actual_mean: f64 = 1.0;
    let actual_std: f64 = 0.3;
    let normal = Normal::new(actual_mean, actual_std).unwrap();
    let mut thread_rng = rand::thread_rng();
    let data: [[f64; 1]; 10_000] = [(); 10_000].map(|_| [normal.sample(&mut thread_rng)]);
    println!("data[..15] = {:?}", &data[..15]);

    // Initial uniform sample for walkers
    let bounds = [[-8.0, 12.0], [1e-2, 1e1]];

    const WALKERS: usize = 64;
    const BURN_IN: usize = 100;
    const SAMPLES: usize = 1000;
    let params = parameter_inference_uniform_prior::<_, D, P, WALKERS, BURN_IN, SAMPLES>(
        &data, &bounds, &loglikelihood,
    );

    // Check that param is close to what we expected
    println!("params = {params:?}");
    assert!(params[0] > 0.9 && params[0] < 1.1 && params[1] > 0.29 && params[1] < 0.31);
}

#[test]
#[should_panic]
fn test_invalid_bound() {
    // Coin flip example, but with floats instead of bools

    // First define bernoulli
    const P: usize = 1;
    const D: usize = 1;
    let loglikelihood = |x: &[f64; 1], p: &[f64; P]| {
        if x[0] == 1.0 {
            p[0].ln()
        } else if x[0] == 0.0 {
            (1.0 - p[0]).ln()
        } else {
            unreachable!("test where we only construct float arrays with 1.0 and 0.0")
        }
    };

    // Then define some data
    let data = [vec![[1.0]; 1000], vec![[0.0]; 500]].concat();

    // Uniform prior (also bounds of parameters)
    let bounds = [[1.0, 0.0]];

    const WALKERS: usize = 32;
    const BURN_IN: usize = 100;
    const SAMPLES: usize = 1000;
    let params = parameter_inference_uniform_prior::<_, D, P, WALKERS, BURN_IN, SAMPLES>(
        &data, &bounds, &loglikelihood,
    );

    // Check that param is close to what we expected
    assert!(params[0] > 0.65 && params[0] < 0.67);
}
