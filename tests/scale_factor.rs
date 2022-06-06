




#[test]
fn test_matter_only() {

    use cosmology::scale_factor::{ScaleFactor, CosmologicalParameters, LITTLE_H_TO_BIG_H};
    
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
    let age_of_universe = 2.0/3.0/hubble;
    let expected = |t:f64| {
        (t / age_of_universe).powf(2.0/3.0)
    };

    // Initialize ScaleFactor
    let mut scale_factor = ScaleFactor::new(
        params,
        z0,
        max_dloga,
        // Obtained via inverting the expected relationship
        Some(age_of_universe * a0.powf(3.0/2.0)), 
    );

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

        // Calculate expected value
        a_expected.push(expected(*t.last().unwrap()));
    }

    // Calculate the avg of L1 loss between the calculated values and the expected values
    let avg_diff =  a
        .iter()
        .zip(a_expected)
        .map(|(&actual, expected)| (actual-expected).abs())
        .sum::<f64>() / a.len() as f64;
    
    // If --nocapture, print value
    println!("avg_diff for matter-only universe {avg_diff:.2e}");

    // Check that the avg of the L1 loss between the calculated values and the expected values
    // is under this threshold
    const ERROR_TOLERANCE: f64 = 1e-10;
    assert!(avg_diff < ERROR_TOLERANCE);
}

#[test]
fn test_radiation_only() {

    use cosmology::scale_factor::{ScaleFactor, CosmologicalParameters, LITTLE_H_TO_BIG_H};
    
    // Specify cosmological parameters
    let params = CosmologicalParameters {
        omega_m0: 0.0,
        omega_de0: 0.0,
        omega_r0: 1.0,
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
    // For a radiation dominated universe, this is a(t) = (t / (1/2/H0))^(1/2)
    let hubble = 0.7 * LITTLE_H_TO_BIG_H;
    let age_of_universe = 1.0/2.0/hubble;
    let expected = |t:f64| {
        (t / age_of_universe).powf(1.0/2.0)
    };

    // Initialize ScaleFactor
    let mut scale_factor = ScaleFactor::new(
        params,
        z0,
        max_dloga,
        // Obtained via inverting the expected relationship
        Some(age_of_universe * a0.powi(2)), 
    );

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

        // Calculate expected value
        a_expected.push(expected(*t.last().unwrap()));
    }

    {
        use plotters::prelude::*;

        let root = BitMapBackend::new("tests/scale-factor-radiation-only.png", (640, 480)).into_drawing_area();
        root.fill(&WHITE).unwrap();

        let mut chart = ChartBuilder::on(&root)
            .caption("Dark Energy Only", ("sans-serif", 50).into_font())
            .margin(5_u32)
            .x_label_area_size(30_u32)
            .y_label_area_size(30_u32)
            .build_cartesian_2d(0.0..age_of_universe * a0.powi(2), 0.0_f64..1.0).unwrap();

        chart.configure_mesh().draw().unwrap();


        chart
            .draw_series(LineSeries::new(
                t.iter().zip(&a).map(|(&x, &y)| (x, y)),
                &RED,
            )).unwrap()
            .label("Calculated")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

        chart
            .draw_series(
                t.iter().zip(&a).step_by(10).map(|(&x, &y)| Circle::new((x, y), 1_u32, BLUE.filled())))
            .unwrap()
            .label("Expected")
            .legend(|(x, y)| Circle::new((x, y), 3_u32, BLUE.filled()));

        chart
            .configure_series_labels()
            .position(SeriesLabelPosition::UpperLeft)
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .draw().unwrap();
    }

    // Calculate the avg of L1 loss between the calculated values and the expected values
    let avg_diff =  a
        .iter()
        .zip(a_expected)
        .map(|(&actual, expected)| (actual-expected).abs())
        .sum::<f64>() / a.len() as f64;
    
    // If --nocapture, print value
    println!("avg_diff for radiation-only universe {avg_diff:.2e}");

    // Check that the avg of the L1 loss between the calculated values and the expected values
    // is under this threshold
    const ERROR_TOLERANCE: f64 = 1e-10;
    assert!(avg_diff < ERROR_TOLERANCE);
}

#[test]
fn test_dark_energy_only() {
    
    use cosmology::scale_factor::{ScaleFactor, CosmologicalParameters, LITTLE_H_TO_BIG_H};
    
    // Specify cosmological parameters
    let params = CosmologicalParameters {
        omega_m0: 0.0,
        omega_de0: 1.0,
        omega_r0: 0.0,
        omega_k0: 0.0,
        w: 1.0,
        h: 0.7,
    };

    // Specify initial redshift
    let z0: f64 = 99.0;
    let a0: f64 = 1.0  / (1.0 + z0);

    // Specify max dloga
    let max_dloga = 0.01;

    // Expectation (analytic solution of the Friedmann equation)
    // For a DE dominated universe, infinite age but a(t) ∝ exp(H_0 t)
    let hubble = 0.7 * LITTLE_H_TO_BIG_H;
    let hubble_time = 1.0 / hubble;
    let t0 = a0.ln() * hubble_time;
    let expected = |t:f64| {
        (t / hubble_time).exp()
    };

    // Initialize ScaleFactor
    let mut scale_factor = ScaleFactor::new(
        params,
        z0,
        max_dloga,
        Some(t0), 
    );

    // Initialize vectors which collect values
    let mut a = vec![scale_factor.get_a()];
    let mut dadt = vec![scale_factor.get_dadt()];
    let mut t = vec![scale_factor.get_time()];
    let dt = 100.0;
    let mut a_expected = vec![scale_factor.get_a()];

    while t.last().unwrap() < &0.0 {
        
        // Evolve scale factor
        scale_factor.step_forward(dt);

        // Add (t, a) to vec
        t.push(scale_factor.get_time());
        a.push(scale_factor.get_a());
        dadt.push(scale_factor.get_dadt());

        // Calculate expected value
        a_expected.push(expected(*t.last().unwrap()));
    }


    {
        use plotters::prelude::*;

        let root = BitMapBackend::new("tests/scale-factor-de-only.png", (640, 480)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let mut chart = ChartBuilder::on(&root)
            .caption("Dark Energy Only", ("sans-serif", 50).into_font())
            .margin(5_u32)
            .x_label_area_size(30_u32)
            .y_label_area_size(30_u32)
            .build_cartesian_2d(t0..0.0, 0.0_f64..1.0).unwrap();

        chart.configure_mesh().draw().unwrap();

        chart
            .draw_series(LineSeries::new(
                t.iter().zip(&a).map(|(&x, &y)| (x, y)),
                &RED,
            )).unwrap()
            .label("Calculated")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

        chart
            .draw_series(
                t.iter().zip(&a).step_by(10).map(|(&x, &y)| Circle::new((x, y), 1_u32, BLUE.filled())))
            .unwrap()
            .label("Expected")
            .legend(|(x, y)| Circle::new((x, y), 3_u32, BLUE.filled()));

        chart
            .configure_series_labels()
            .position(SeriesLabelPosition::UpperLeft)
            .background_style(&WHITE.mix(0.8))
            .border_style(&BLACK)
            .draw().unwrap();
    }

    // Calculate the avg of L1 loss between the calculated values and the expected values
    let avg_diff =  a
        .iter()
        .zip(a_expected)
        .map(|(&actual, expected)| (actual-expected).abs())
        .sum::<f64>() / a.len() as f64;

    
    // If --nocapture, print value
    println!("avg_diff for dark energy-only universe {avg_diff:.2e}");

    // Check that the avg of the L1 loss between the calculated values and the expected values
    // is under this threshold
    const ERROR_TOLERANCE: f64 = 1e-10;
    assert!(avg_diff < ERROR_TOLERANCE);
}