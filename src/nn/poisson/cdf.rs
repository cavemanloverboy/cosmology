pub(super) fn calculate_poisson_cdf(r: f64, nbar: f64, k: u8) -> f64 {
    // Convert distance to corresponding spherical volume
    let v = 4.0 * std::f64::consts::PI * r.powi(3) / 3.0;

    match k {
        1 => 1. - 1. / exp(nbar * v),
        2 => (-1. + exp(nbar * v) - nbar * v) / exp(nbar * v),
        3 => {
            (0.5 * (-2. + 2. * exp(nbar * v) - 2. * nbar * v - (nbar * v).powi(2))) / exp(nbar * v)
        }
        4 => {
            (0.16666666666666666
                * (-6. + 6. * exp(nbar * v)
                    - 6. * nbar * v
                    - 3. * (nbar * v).powi(2)
                    - (nbar * v).powi(3)))
                / exp(nbar * v)
        }
        5 => {
            (0.041666666666666664
                * (-24. + 24. * exp(nbar * v)
                    - 24. * nbar * v
                    - 12. * (nbar * v).powi(2)
                    - 4. * (nbar * v).powi(3)
                    - (nbar * v).powi(4)))
                / exp(nbar * v)
        }
        6 => {
            (0.008333333333333333
                * (-120. + 120. * exp(nbar * v)
                    - 120. * nbar * v
                    - 60. * (nbar * v).powi(2)
                    - 20. * (nbar * v).powi(3)
                    - 5. * (nbar * v).powi(4)
                    - (nbar * v).powi(5)))
                / exp(nbar * v)
        }
        7 => {
            (0.001388888888888889
                * (-720. + 720. * exp(nbar * v)
                    - 720. * nbar * v
                    - 360. * (nbar * v).powi(2)
                    - 120. * (nbar * v).powi(3)
                    - 30. * (nbar * v).powi(4)
                    - 6. * (nbar * v).powi(5)
                    - (nbar * v).powi(6)))
                / exp(nbar * v)
        }
        8 => {
            (0.0001984126984126984
                * (-5040. + 5040. * exp(nbar * v)
                    - 5040. * nbar * v
                    - 2520. * (nbar * v).powi(2)
                    - 840. * (nbar * v).powi(3)
                    - 210. * (nbar * v).powi(4)
                    - 42. * (nbar * v).powi(5)
                    - 7. * (nbar * v).powi(6)
                    - (nbar * v).powi(7)))
                / exp(nbar * v)
        }
        9 => {
            (0.0000248015873015873
                * (-40320. + 40320. * exp(nbar * v)
                    - 40320. * nbar * v
                    - 20160. * (nbar * v).powi(2)
                    - 6720. * (nbar * v).powi(3)
                    - 1680. * (nbar * v).powi(4)
                    - 336. * (nbar * v).powi(5)
                    - 56. * (nbar * v).powi(6)
                    - 8. * (nbar * v).powi(7)
                    - (nbar * v).powi(8)))
                / exp(nbar * v)
        }
        10 => {
            (2.7557319223985893e-6
                * (-362880. + 362880. * exp(nbar * v)
                    - 362880. * nbar * v
                    - 181440. * (nbar * v).powi(2)
                    - 60480. * (nbar * v).powi(3)
                    - 15120. * (nbar * v).powi(4)
                    - 3024. * (nbar * v).powi(5)
                    - 504. * (nbar * v).powi(6)
                    - 72. * (nbar * v).powi(7)
                    - 9. * (nbar * v).powi(8)
                    - (nbar * v).powi(9)))
                / exp(nbar * v)
        }
        11 => {
            (2.755731922398589e-7
                * (-3.6288e6 + 3.6288e6 * exp(nbar * v)
                    - 3.6288e6 * nbar * v
                    - 1.8144e6 * (nbar * v).powi(2)
                    - 604800. * (nbar * v).powi(3)
                    - 151200. * (nbar * v).powi(4)
                    - 30240. * (nbar * v).powi(5)
                    - 5040. * (nbar * v).powi(6)
                    - 720. * (nbar * v).powi(7)
                    - 90. * (nbar * v).powi(8)
                    - 10. * (nbar * v).powi(9)
                    - (nbar * v).powi(10)))
                / exp(nbar * v)
        }
        12 => {
            (2.505210838544172e-8
                * (-3.99168e7 + 3.99168e7 * exp(nbar * v)
                    - 3.99168e7 * nbar * v
                    - 1.99584e7 * (nbar * v).powi(2)
                    - 6.6528e6 * (nbar * v).powi(3)
                    - 1.6632e6 * (nbar * v).powi(4)
                    - 332640. * (nbar * v).powi(5)
                    - 55440. * (nbar * v).powi(6)
                    - 7920. * (nbar * v).powi(7)
                    - 990. * (nbar * v).powi(8)
                    - 110. * (nbar * v).powi(9)
                    - 11. * (nbar * v).powi(10)
                    - (nbar * v).powi(11)))
                / exp(nbar * v)
        }
        13 => {
            (2.08767569878681e-9
                * (-4.790016e8 + 4.790016e8 * exp(nbar * v)
                    - 4.790016e8 * nbar * v
                    - 2.395008e8 * (nbar * v).powi(2)
                    - 7.98336e7 * (nbar * v).powi(3)
                    - 1.99584e7 * (nbar * v).powi(4)
                    - 3.99168e6 * (nbar * v).powi(5)
                    - 665280. * (nbar * v).powi(6)
                    - 95040. * (nbar * v).powi(7)
                    - 11880. * (nbar * v).powi(8)
                    - 1320. * (nbar * v).powi(9)
                    - 132. * (nbar * v).powi(10)
                    - 12. * (nbar * v).powi(11)
                    - (nbar * v).powi(12)))
                / exp(nbar * v)
        }
        14 => {
            (1.6059043836821613e-10
                * (-6.2270208e9 + 6.2270208e9 * exp(nbar * v)
                    - 6.2270208e9 * nbar * v
                    - 3.1135104e9 * (nbar * v).powi(2)
                    - 1.0378368e9 * (nbar * v).powi(3)
                    - 2.594592e8 * (nbar * v).powi(4)
                    - 5.189184e7 * (nbar * v).powi(5)
                    - 8.64864e6 * (nbar * v).powi(6)
                    - 1.23552e6 * (nbar * v).powi(7)
                    - 154440. * (nbar * v).powi(8)
                    - 17160. * (nbar * v).powi(9)
                    - 1716. * (nbar * v).powi(10)
                    - 156. * (nbar * v).powi(11)
                    - 13. * (nbar * v).powi(12)
                    - (nbar * v).powi(13)))
                / exp(nbar * v)
        }
        15 => {
            (1.1470745597729725e-11
                * (-8.71782912e10 + 8.71782912e10 * exp(nbar * v)
                    - 8.71782912e10 * nbar * v
                    - 4.35891456e10 * (nbar * v).powi(2)
                    - 1.45297152e10 * (nbar * v).powi(3)
                    - 3.6324288e9 * (nbar * v).powi(4)
                    - 7.2648576e8 * (nbar * v).powi(5)
                    - 1.2108096e8 * (nbar * v).powi(6)
                    - 1.729728e7 * (nbar * v).powi(7)
                    - 2.16216e6 * (nbar * v).powi(8)
                    - 240240. * (nbar * v).powi(9)
                    - 24024. * (nbar * v).powi(10)
                    - 2184. * (nbar * v).powi(11)
                    - 182. * (nbar * v).powi(12)
                    - 14. * (nbar * v).powi(13)
                    - (nbar * v).powi(14)))
                / exp(nbar * v)
        }
        16 => {
            (7.647163731819816e-13
                * (-1.307674368e12 + 1.307674368e12 * exp(nbar * v)
                    - 1.307674368e12 * nbar * v
                    - 6.53837184e11 * (nbar * v).powi(2)
                    - 2.17945728e11 * (nbar * v).powi(3)
                    - 5.4486432e10 * (nbar * v).powi(4)
                    - 1.08972864e10 * (nbar * v).powi(5)
                    - 1.8162144e9 * (nbar * v).powi(6)
                    - 2.594592e8 * (nbar * v).powi(7)
                    - 3.24324e7 * (nbar * v).powi(8)
                    - 3.6036e6 * (nbar * v).powi(9)
                    - 360360. * (nbar * v).powi(10)
                    - 32760. * (nbar * v).powi(11)
                    - 2730. * (nbar * v).powi(12)
                    - 210. * (nbar * v).powi(13)
                    - 15. * (nbar * v).powi(14)
                    - (nbar * v).powi(15)))
                / exp(nbar * v)
        }

        _ => todo!("only k=1 through k=16 NN CDFs are implemented for the Gaussian Random Fields"),
    }
}
/// This is just here to simplify the conversion from mathematica output
/// as simply as possible, i.e. without any regex (or similar) manipulation.
#[inline(always)]
fn exp(x: f64) -> f64 {
    x.exp()
}
