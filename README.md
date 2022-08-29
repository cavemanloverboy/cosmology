# cosmology: a rust crate for cosmology

An early-in-development crate intended to eventually include lots of utilities commonly used in cosmology, including

### Current Features

- [x] Solves the Friedmann equations to get the scale factor for a given cosmology over time.
  - [x] On-the-fly: supply some `dt` and you are able to `step_forward(dt)`. Useful alongside another algorithm (e.g. nbody or fluid code) where timesteps are not known in advance.
  - [x] In advance: supply some `dt` spacing and a final time or scale_factor and get a `Box<[T]>` time series for `a` and `dadt`. Useful when evaluating on a time grid known in advance.

### Planned Features

- [ ] Power spectra `P(k)` and correlation function `ξ(r)` using transfer function from (Eisenstein & Hu '98).
- [ ] Cosmological initial condition generators.
- [ ] Linear theory predictions for the kNN distributions of matter and its tracers (e.g. haloes, galaxies).

### Sample Usage

```rust
use cosmology::scale_factor::{CosmologicalParameters, ScaleFactor, LITTLE_H_TO_BIG_H};

// Specify cosmological parameters
let params = CosmologicalParameters {
    omega_m0: 1.0,
    omega_de0: 0.0,
    omega_r0: 0.0,
    omega_k0: 0.0,
    w: 1.0,
    h: 0.7,
};

// Specify initial redshift
let z0 = 9.0;
let a0 = 1.0 / (1.0 + z0);

// Specify max dloga
let max_dloga = 0.01;

// Expectation (analytic solution of the Friedmann equation)
// For a matter dominated universe, this is a(t) = (t / (2/3/H0))^(2/3)
let hubble = 0.7 * LITTLE_H_TO_BIG_H;
let age_of_universe = 2.0 / 3.0 / hubble;
let expected = |t: f64| (t / age_of_universe).powf(2.0 / 3.0);

// Initialize ScaleFactor
let mut scale_factor = ScaleFactor::new(
    params,
    z0,
    max_dloga,
    // Obtained via inverting the expected relationship
    Some(age_of_universe * a0.powf(3.0 / 2.0)),
);

// On-the-fly
// Initialize vectors which collect values
let mut a = vec![scale_factor.get_a()];
let mut dadt = vec![scale_factor.get_dadt()];
let mut t = vec![scale_factor.get_time()];
let dt = 100.0;
let mut a_expected = vec![scale_factor.get_a()];

while a.last().unwrap() < &1.0 {
    // Evolve scale factor
    scale_factor.step_forward(dt);

    // Add (t, a) to vec
    t.push(scale_factor.get_time());
    a.push(scale_factor.get_a());
    dadt.push(scale_factor.get_dadt());

    // Here you can do something that requires a, dadt, at your specified time t.
    // do_something(&a, &dadt);

    // Calculate expected value
    a_expected.push(expected(*t.last().unwrap()));
}

// Calculate the avg of L1 loss between the calculated values and the expected values
let avg_diff = a
    .iter()
    .zip(a_expected)
    .map(|(&actual, expected)| (actual - expected).abs())
    .sum::<f64>()
    / a.len() as f64;

// If --nocapture, print value
println!("avg_diff for matter-only universe {avg_diff:.2e}");

// Check that the avg of the L1 loss between the calculated values and the expected values
// is under this threshold
const ERROR_TOLERANCE: f64 = 1e-10;
assert!(avg_diff < ERROR_TOLERANCE);
```

# About the Author

I work on a lot of things including non-equilibrium fluid dynamics, cosmology (large scale structure and scalar field dark matter), machine learning, distributed systems, smart contracts, (financial) derivatives. I am a Rust zealot and was surprised to see that there is a gap to be filled in the astrophysics/cosmology crate space. This is a crate that I am looking forward to developing.

I first began to look for a crate like this one because I needed to solve the Friedmann equations for a cosmological scalar field solver I've been working on. I found no such crate. So, here I am building it. This is presently in a very early stage but I look very forward to adding features over time. Feel free to submit PRs or suggestions for new features.
