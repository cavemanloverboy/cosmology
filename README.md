# cosmology: a rust crate for cosmology

Hey y'all, it's cavemanloverboy. I work on a lot of things including non-equilibrium fluid dynamics, cosmology (large scale structure and scalar field dark matter), machine learning, distributed systems, smart contracts, (financial) derivatives. I am a Rust zealot and was surprised to see that there is a gap to be filled in the astrophysics/cosmology crate space. This is a crate that I am looking forward to developing.

I first began to look for such a crate because I needed to solve the Friedmann equations for a cosmological scalar field solver I've been working on. I found no such crate. So, here I am building it. This is presently in a very early stage but I look very forward to adding features over time. Feel free to submit PRs or suggestions for new features.

###Current Features

- [x] Solves the Friedmann equations to get the scale factor for a given cosmology over time.
  - [x] On-the-fly: supply some `dt` and you are able to `step_forward(dt)`. This is what I will be using in my quantum solver since the `dt` taken by the solver at some n-th step is not known at compile time.
  - [x] In advance: supply some `dt` spacing and a final time or scale_factor and get a `Box<[T]>` time series for `a` and `dadt`

### Planned Features

- [ ] Implement power spectra `P(k)` and correlation function `ξ(r)` using transfer function from (Eisenstein & Hu '98).
