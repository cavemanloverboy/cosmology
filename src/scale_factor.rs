type Time = f64;

/// The primary, user-facing struct containing the state of the scale factor,
/// cosmological parameters, and the integrator configuration.
pub struct ScaleFactor {
    /// Cosmological parameters for this universe
    params: CosmologicalParameters,

    /// Current state
    state: ScaleFactorState,

    /// Integrator Parameters
    config: IntegratorConfig,
}

/// Multiplicative factor used to turn h to H in units of 1/Myr.
/// (i.e. 100 km/s/Mpc converted to 1/Myr)
pub const LITTLE_H_TO_BIG_H: f64 = 1.022e-4;

/// This contains all of present day values of the cosmological parameters
/// relevant to the evolution of the scale factor, as per the Friedmann Equations.
pub struct CosmologicalParameters {
    /// Present day value of Ω_m0, the matter content of the universe
    pub omega_m0: f64,

    /// Present day value of Ω_Λ0, the dark enery content of the universe
    pub omega_de0: f64,

    /// Present day value of Ω_r0, the radiation/relativistic content of the universe
    pub omega_r0: f64,

    /// Present day value of Ω_k0, the energy associated with the universe's curvature
    pub omega_k0: f64,

    /// The EoS parameter for dark energy. TODO: implement
    #[allow(dead_code)]
    pub w: f64,

    /// The present value of the hubble parameter, as little h.
    pub h: f64,
}

/// A struct containing the state of the universe. This includes the scale factor,
/// its derivetive, and the current time.
struct ScaleFactorState {
    /// scale factor
    a: f64,
    /// Time derivative of scale factor
    dadt: f64,
    /// Time
    t: f64,
    /// Flag whether to use the positive or negative root of the Friedmann Eq
    expanding: bool,
}

/// A struct containing parameters about the integrator
struct IntegratorConfig {
    max_dloga: f64,
}

impl ScaleFactor {
    /// A constructor function forour
    pub fn new(params: CosmologicalParameters, z0: f64, max_dloga: f64, t0: Option<f64>) -> Self {
        // Construct initial state
        let state = ScaleFactorState {
            // Convert redshift to initial scale factor
            a: 1.0 / (1.0 + z0),
            // Calculate initial derivative dadt
            dadt: Self::derivative_associated(1.0 / (1.0 + z0), &params, true),
            // Initial time. Default is t = 1 age of universe at z0.
            t: t0.unwrap_or(1.0),
            // TODO: implement contracting universes
            expanding: true,
        };

        // Construct integrator config
        let config = IntegratorConfig { max_dloga };

        Self {
            params,
            state,
            config,
        }
    }

    /// Updates the time derivative of the scale factor in the `ScaleFactorState` struct.
    pub fn update_dadt(&mut self) {
        self.state.dadt =
            Self::derivative_associated(self.state.a, &self.params, self.state.expanding);
    }

    pub fn derivative(&self) -> f64 {
        Self::derivative_associated(self.state.a, &self.params, self.state.expanding)
    }

    fn derivative_associated(a: f64, params: &CosmologicalParameters, expanding: bool) -> f64 {
        params.h
            * LITTLE_H_TO_BIG_H
            * (params.omega_r0 / a.powi(2)
                + params.omega_m0 / a.powi(1)
                + params.omega_de0 * a.powi(2)
                + params.omega_k0)
                .sqrt()
            * if expanding { 1.0 } else { -1.0 }
    }

    /// For on the fly evaluation of the scale factor. This steps forward by `dt`, but breaks the interval
    /// into subintervals if dloga is above the user threshold `max_dloga` to linear order.
    pub fn step_forward(&mut self, dt: f64) {
        // Update dadt
        self.update_dadt();

        // Simple first-order check for size of dloga
        let mut steps: usize = 1;
        let mut step_dt: f64 = dt;
        while (self.state.dadt / self.state.a * step_dt) > self.config.max_dloga {
            steps = steps.checked_add(1).unwrap();
            step_dt = dt / steps as f64;
        }
        let final_time = self.state.t + dt;

        while steps > 0 {
            steps -= 1;

            // Update derivative
            self.update_dadt();
            // println!("time derivative is {}; steps is {steps}, a is {} and t is {}", self.state.dadt, self.state.a, self.state.t);

            // Derivative function in its standard form for RK4
            let f = |_tn, an| Self::derivative_associated(an, &self.params, self.state.expanding);

            // Final step
            if steps == 0 {
                // Slightly more exact (re: floating point)
                let step_dt = final_time - self.state.t;

                // Update value
                self.state.a = rk4(
                    f,
                    self.state.t,
                    self.state.a,
                    step_dt,
                    Some(self.state.dadt),
                );

                // Update time
                self.state.t = final_time;
            } else {
                // Update value
                self.state.a = rk4(
                    f,
                    self.state.t,
                    self.state.a,
                    step_dt,
                    Some(self.state.dadt),
                );

                // Update time
                self.state.t += step_dt;
            }
        }

        // Update derivative one more time before exiting to be at t_n+1
        self.update_dadt();
    }

    /// Returns the currently stored value for the scale factor
    pub fn get_a(&self) -> f64 {
        self.state.a
    }

    /// Returns the currently stored value for the time derivative of the scale factor
    pub fn get_dadt(&self) -> f64 {
        self.state.dadt
    }

    /// Returns the currently stored value for the time
    pub fn get_time(&self) -> f64 {
        self.state.t
    }

    /// Gather time series for the evolution of the scale factor with spacing `dt`.
    pub fn get_time_series(&mut self, dt: f64, terminate: Terminate) -> TimeSeries {
        // Initialize vectors which collect values
        let mut a = vec![self.get_a()];
        let mut dadt = vec![self.get_dadt()];
        let mut t = vec![self.get_time()];

        let condition = |a: &Vec<f64>, t: &Vec<Time>| match terminate {
            Terminate::ScaleFactor(ref final_a) => a.last().unwrap() < final_a,
            Terminate::Time(ref final_time) => t.last().unwrap() < final_time,
        };

        while condition(&a, &t) {
            // Evolve scale factor
            self.step_forward(dt);

            // Add (t, a) to vec
            t.push(self.get_time());
            a.push(self.get_a());
            dadt.push(self.get_dadt());
        }

        TimeSeries {
            t: t.into_boxed_slice(),
            a: a.into_boxed_slice(),
            dadt: dadt.into_boxed_slice(),
        }
    }
}

/// This returns the timeseries for a given
pub struct TimeSeries {
    pub t: Box<[Time]>,
    pub a: Box<[f64]>,
    pub dadt: Box<[f64]>,
}

pub enum Terminate {
    ScaleFactor(f64),
    Time(f64),
}

/// This function does a single step of the rk4 algorithm in 1D.
/// It is tested indirectly through the tests for solving the Friedmann Equation
/// for the evolution of the scale factor.
pub(crate) fn rk4<F>(
    // Function that evaluates derivative
    f: F,
    // Current time
    tn: Time,
    // Current function value
    yn: f64,
    // Timestep
    h: f64,
    // Option to provide already-computed derivative at tn (i.e. f(tn,yn))
    derivative: Option<f64>,
) -> f64
where
    F: Fn(Time, f64) -> f64,
{
    // Evaluate the derivative function at (tn, yn)
    let k1: f64 = derivative.unwrap_or_else(|| f(tn, yn));

    // Evaluate the derivative funciton at (tn + h/2, y + h*k1/2)
    let k2: f64 = f(tn + h / 2.0, yn + h * k1 / 2.0);

    // Evaluate the derivative funciton at (tn + h/2, y + h*k2/2)
    let k3: f64 = f(tn + h / 2.0, yn + h * k2 / 2.0);

    // Evaluate the derivative funciton at (tn + h, y + h*k3)
    let k4: f64 = f(tn + h, yn + h * k3);

    // Return new value as per RK4
    yn + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
}

/// This function does a single step of the rk4 algorithm in 1D.
/// It is tested indirectly through the tests for solving the Friedmann Equation
/// for the evolution of the scale factor.
pub(crate) fn rk4_multi<F, const N: usize>(
    // Function that evaluates derivative
    f: F,
    // Current time
    tn: Time,
    // Current function values
    yn: [f64; N],
    // Timestep
    h: f64,
) -> [f64; N]
where
    F: Fn(Time, [f64; N]) -> [f64; N],
{
    // Initialize mutable arrays (on the stack)
    let mut scratch: [f64; N] = [0.0; N];

    // Evaluate the derivative function at (tn, yn)
    let k1: [f64; N] = f(tn, yn);

    // Evaluate the derivative funciton at (tn + h/2, y + h*k1/2)
    for i in 0..N {
        // Here, scratch is yn + h * k1 / 2.0
       scratch[i] = yn[i] + h * k1[i] / 2.0;
    }
    let k2: [f64; N] = f(tn + h / 2.0, scratch);

    // Evaluate the derivative funciton at (tn + h/2, y + h*k2/2)
    for i in 0..N {
        // Here, scratch is yn + h * k2 / 2.0
       scratch[i] = yn[i] + h * k2[i] / 2.0;
    }
    let k3: [f64; N] = f(tn + h / 2.0, scratch);

    // Evaluate the derivative function at (tn + h, y + h*k3)
    for i in 0..N {
        // Here, scratch is yn + h * k3
       scratch[i] = yn[i] + h * k3[i];
    }
    let k4: [f64; N] = f(tn + h, scratch);

    // Return new value as per RK4
    for i in 0..N {
        // Here, scratch is the Rk4 result
       scratch[i] = yn[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
    }
    scratch
}


#[test]
fn test_rk4_impls_are_consistent() {

    // y = t^2 --> y' = 2t.
    let t0 = 0.0;
    let y0 = 0.0;
    let f = |t: f64, _y0: f64| 2.0 * t;

    let y0_multi = [0.0];
    let f_multi = |t: f64, _y0: [f64; 1]| [2.0 * t];

    let dt = 1.0;

    let first_result = rk4(f, t0, y0, dt, None);
    let second_result = rk4_multi(f_multi, t0, y0_multi, dt);

    assert_eq!(first_result, second_result[0]);

}

