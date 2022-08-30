

use criterion::{Criterion, criterion_group, criterion_main, black_box};


const LOGK_MIN: f64 = -2.0;
const LOGK_MAX: f64 = 2.0;
const NUM_K_POINTS: usize = 1000;

fn bench_power(c: &mut Criterion) {

    // Initialize power engine
    let power = cosmology::power::PowerSpectrum::new(
        cosmology::power::TransferFunction::EisensteinHu {
            h: 0.7,
            omega_matter_0: 0.3,
            omega_baryon_0: 0.2, 
            temp_cmb0: 2.7,
            ns: 0.96,
            sigma_8: 0.6,
        }
    ).unwrap();

    // Initialize wavenumbers of interest
    let ks: Vec<f64> = (0..NUM_K_POINTS)
        .map(|i| 10_f64.powf(LOGK_MIN + (i as f64) * (LOGK_MAX-LOGK_MIN)/(NUM_K_POINTS as f64 - 1.0)))
        .collect();


    c.bench_function("power", move |b| {
        b.iter(|| {
            black_box(power.calculate_power(black_box(&ks), black_box(1.0)).unwrap());
        })
    });
    
}




criterion_group!(power, bench_power);
criterion_main!(power);