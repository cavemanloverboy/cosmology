use std::f64::consts::PI;

#[allow(unused)]
fn calculate_poisson_pdf(r: f64, nbar: f64, k: u8) -> f64 {
    // Convert distance to corresponding spherical volume
    let v = 4.0 * std::f64::consts::PI * r.powi(3) / 3.0;

    // int P(V) dV = 1
    // int P(V(R)) * dR dV/dR = 1
    // ----> P(R) = dV/dR * P(V)
    // ----> P(R) = 4 * PI * R^2 * dCDF(V)/dV
    let jacobian = 4.0 * PI * r.powi(2);

    jacobian
        * match k {
            1 => nbar / exp(nbar * v),
            2 => {
                (-1. * nbar + exp(nbar * v) * nbar) / exp(nbar * v)
                    - (nbar * (-1. + exp(nbar * v) - 1. * nbar * v)) / exp(nbar * v)
            }
            3 => {
                (0.5 * (-2. * nbar + 2. * exp(nbar * v) * nbar - 2. * nbar.powi(2) * v))
                    / exp(nbar * v)
                    - (0.5
                        * nbar
                        * (-2. + 2. * exp(nbar * v)
                            - 2. * nbar * v
                            - 1. * nbar.powi(2) * v.powi(2)))
                        / exp(nbar * v)
            }
            4 => {
                (0.16666666666666666
                    * (-6. * nbar + 6. * exp(nbar * v) * nbar
                        - 6. * nbar.powi(2) * v
                        - 3. * nbar * (nbar * v).powi(2)))
                    / exp(nbar * v)
                    - (0.16666666666666666
                        * nbar
                        * (-6. + 6. * exp(nbar * v)
                            - 6. * nbar * v
                            - 3. * nbar.powi(2) * v.powi(2)
                            - 1. * nbar.powi(3) * v.powi(3)))
                        / exp(nbar * v)
            }
            5 => {
                (0.041666666666666664
                    * (-24. * nbar + 24. * exp(nbar * v) * nbar
                        - 24. * nbar.powi(2) * v
                        - 12. * nbar * (nbar * v).powi(2)
                        - 4. * nbar * (nbar * v).powi(3)))
                    / exp(nbar * v)
                    - (0.041666666666666664
                        * nbar
                        * (-24. + 24. * exp(nbar * v)
                            - 24. * nbar * v
                            - 12. * nbar.powi(2) * v.powi(2)
                            - 4. * nbar.powi(3) * v.powi(3)
                            - 1. * nbar.powi(4) * v.powi(4)))
                        / exp(nbar * v)
            }
            6 => {
                (0.008333333333333333
                    * (-120. * nbar + 120. * exp(nbar * v) * nbar
                        - 120. * nbar.powi(2) * v
                        - 60. * nbar * (nbar * v).powi(2)
                        - 20. * nbar * (nbar * v).powi(3)
                        - 5. * nbar * (nbar * v).powi(4)))
                    / exp(nbar * v)
                    - (0.008333333333333333
                        * nbar
                        * (-120. + 120. * exp(nbar * v)
                            - 120. * nbar * v
                            - 60. * nbar.powi(2) * v.powi(2)
                            - 20. * nbar.powi(3) * v.powi(3)
                            - 5. * nbar.powi(4) * v.powi(4)
                            - 1. * nbar.powi(5) * v.powi(5)))
                        / exp(nbar * v)
            }
            7 => {
                (0.001388888888888889
                    * (-720. * nbar + 720. * exp(nbar * v) * nbar
                        - 720. * nbar.powi(2) * v
                        - 360. * nbar * (nbar * v).powi(2)
                        - 120. * nbar * (nbar * v).powi(3)
                        - 30. * nbar * (nbar * v).powi(4)
                        - 6. * nbar * (nbar * v).powi(5)))
                    / exp(nbar * v)
                    - (0.001388888888888889
                        * nbar
                        * (-720. + 720. * exp(nbar * v)
                            - 720. * nbar * v
                            - 360. * nbar.powi(2) * v.powi(2)
                            - 120. * nbar.powi(3) * v.powi(3)
                            - 30. * nbar.powi(4) * v.powi(4)
                            - 6. * nbar.powi(5) * v.powi(5)
                            - 1. * nbar.powi(6) * v.powi(6)))
                        / exp(nbar * v)
            }
            8 => {
                (0.0001984126984126984
                    * (-5040. * nbar + 5040. * exp(nbar * v) * nbar
                        - 5040. * nbar.powi(2) * v
                        - 2520. * nbar * (nbar * v).powi(2)
                        - 840. * nbar * (nbar * v).powi(3)
                        - 210. * nbar * (nbar * v).powi(4)
                        - 42. * nbar * (nbar * v).powi(5)
                        - 7. * nbar * (nbar * v).powi(6)))
                    / exp(nbar * v)
                    - (0.0001984126984126984
                        * nbar
                        * (-5040. + 5040. * exp(nbar * v)
                            - 5040. * nbar * v
                            - 2520. * nbar.powi(2) * v.powi(2)
                            - 840. * nbar.powi(3) * v.powi(3)
                            - 210. * nbar.powi(4) * v.powi(4)
                            - 42. * nbar.powi(5) * v.powi(5)
                            - 7. * nbar.powi(6) * v.powi(6)
                            - 1. * nbar.powi(7) * v.powi(7)))
                        / exp(nbar * v)
            }
            9 => {
                (0.0000248015873015873
                    * (-40320. * nbar + 40320. * exp(nbar * v) * nbar
                        - 40320. * nbar.powi(2) * v
                        - 20160. * nbar * (nbar * v).powi(2)
                        - 6720. * nbar * (nbar * v).powi(3)
                        - 1680. * nbar * (nbar * v).powi(4)
                        - 336. * nbar * (nbar * v).powi(5)
                        - 56. * nbar * (nbar * v).powi(6)
                        - 8. * nbar * (nbar * v).powi(7)))
                    / exp(nbar * v)
                    - (0.0000248015873015873
                        * nbar
                        * (-40320. + 40320. * exp(nbar * v)
                            - 40320. * nbar * v
                            - 20160. * nbar.powi(2) * v.powi(2)
                            - 6720. * nbar.powi(3) * v.powi(3)
                            - 1680. * nbar.powi(4) * v.powi(4)
                            - 336. * nbar.powi(5) * v.powi(5)
                            - 56. * nbar.powi(6) * v.powi(6)
                            - 8. * nbar.powi(7) * v.powi(7)
                            - 1. * nbar.powi(8) * v.powi(8)))
                        / exp(nbar * v)
            }
            10 => {
                (2.7557319223985893e-6
                    * (-362880. * nbar + 362880. * exp(nbar * v) * nbar
                        - 362880. * nbar.powi(2) * v
                        - 181440. * nbar * (nbar * v).powi(2)
                        - 60480. * nbar * (nbar * v).powi(3)
                        - 15120. * nbar * (nbar * v).powi(4)
                        - 3024. * nbar * (nbar * v).powi(5)
                        - 504. * nbar * (nbar * v).powi(6)
                        - 72. * nbar * (nbar * v).powi(7)
                        - 9. * nbar * (nbar * v).powi(8)))
                    / exp(nbar * v)
                    - (2.7557319223985893e-6
                        * nbar
                        * (-362880. + 362880. * exp(nbar * v)
                            - 362880. * nbar * v
                            - 181440. * nbar.powi(2) * v.powi(2)
                            - 60480. * nbar.powi(3) * v.powi(3)
                            - 15120. * nbar.powi(4) * v.powi(4)
                            - 3024. * nbar.powi(5) * v.powi(5)
                            - 504. * nbar.powi(6) * v.powi(6)
                            - 72. * nbar.powi(7) * v.powi(7)
                            - 9. * nbar.powi(8) * v.powi(8)
                            - 1. * nbar.powi(9) * v.powi(9)))
                        / exp(nbar * v)
            }
            11 => {
                (2.755731922398589e-7
                    * (-3.6288e6 * nbar + 3.6288e6 * exp(nbar * v) * nbar
                        - 3.6288e6 * nbar.powi(2) * v
                        - 1.8144e6 * nbar * (nbar * v).powi(2)
                        - 604800. * nbar * (nbar * v).powi(3)
                        - 151200. * nbar * (nbar * v).powi(4)
                        - 30240. * nbar * (nbar * v).powi(5)
                        - 5040. * nbar * (nbar * v).powi(6)
                        - 720. * nbar * (nbar * v).powi(7)
                        - 90. * nbar * (nbar * v).powi(8)
                        - 10. * nbar * (nbar * v).powi(9)))
                    / exp(nbar * v)
                    - (2.755731922398589e-7
                        * nbar
                        * (-3.6288e6 + 3.6288e6 * exp(nbar * v)
                            - 3.6288e6 * nbar * v
                            - 1.8144e6 * nbar.powi(2) * v.powi(2)
                            - 604800. * nbar.powi(3) * v.powi(3)
                            - 151200. * nbar.powi(4) * v.powi(4)
                            - 30240. * nbar.powi(5) * v.powi(5)
                            - 5040. * nbar.powi(6) * v.powi(6)
                            - 720. * nbar.powi(7) * v.powi(7)
                            - 90. * nbar.powi(8) * v.powi(8)
                            - 10. * nbar.powi(9) * v.powi(9)
                            - 1. * nbar.powi(10) * v.powi(10)))
                        / exp(nbar * v)
            }
            12 => {
                (2.505210838544172e-8
                    * (-3.99168e7 * nbar + 3.99168e7 * exp(nbar * v) * nbar
                        - 3.99168e7 * nbar.powi(2) * v
                        - 1.99584e7 * nbar * (nbar * v).powi(2)
                        - 6.6528e6 * nbar * (nbar * v).powi(3)
                        - 1.6632e6 * nbar * (nbar * v).powi(4)
                        - 332640. * nbar * (nbar * v).powi(5)
                        - 55440. * nbar * (nbar * v).powi(6)
                        - 7920. * nbar * (nbar * v).powi(7)
                        - 990. * nbar * (nbar * v).powi(8)
                        - 110. * nbar * (nbar * v).powi(9)
                        - 11. * nbar * (nbar * v).powi(10)))
                    / exp(nbar * v)
                    - (2.505210838544172e-8
                        * nbar
                        * (-3.99168e7 + 3.99168e7 * exp(nbar * v)
                            - 3.99168e7 * nbar * v
                            - 1.99584e7 * nbar.powi(2) * v.powi(2)
                            - 6.6528e6 * nbar.powi(3) * v.powi(3)
                            - 1.6632e6 * nbar.powi(4) * v.powi(4)
                            - 332640. * nbar.powi(5) * v.powi(5)
                            - 55440. * nbar.powi(6) * v.powi(6)
                            - 7920. * nbar.powi(7) * v.powi(7)
                            - 990. * nbar.powi(8) * v.powi(8)
                            - 110. * nbar.powi(9) * v.powi(9)
                            - 11. * nbar.powi(10) * v.powi(10)
                            - 1. * nbar.powi(11) * v.powi(11)))
                        / exp(nbar * v)
            }
            13 => {
                (2.08767569878681e-9
                    * (-4.790016e8 * nbar + 4.790016e8 * exp(nbar * v) * nbar
                        - 4.790016e8 * nbar.powi(2) * v
                        - 2.395008e8 * nbar * (nbar * v).powi(2)
                        - 7.98336e7 * nbar * (nbar * v).powi(3)
                        - 1.99584e7 * nbar * (nbar * v).powi(4)
                        - 3.99168e6 * nbar * (nbar * v).powi(5)
                        - 665280. * nbar * (nbar * v).powi(6)
                        - 95040. * nbar * (nbar * v).powi(7)
                        - 11880. * nbar * (nbar * v).powi(8)
                        - 1320. * nbar * (nbar * v).powi(9)
                        - 132. * nbar * (nbar * v).powi(10)
                        - 12. * nbar * (nbar * v).powi(11)))
                    / exp(nbar * v)
                    - (2.08767569878681e-9
                        * nbar
                        * (-4.790016e8 + 4.790016e8 * exp(nbar * v)
                            - 4.790016e8 * nbar * v
                            - 2.395008e8 * nbar.powi(2) * v.powi(2)
                            - 7.98336e7 * nbar.powi(3) * v.powi(3)
                            - 1.99584e7 * nbar.powi(4) * v.powi(4)
                            - 3.99168e6 * nbar.powi(5) * v.powi(5)
                            - 665280. * nbar.powi(6) * v.powi(6)
                            - 95040. * nbar.powi(7) * v.powi(7)
                            - 11880. * nbar.powi(8) * v.powi(8)
                            - 1320. * nbar.powi(9) * v.powi(9)
                            - 132. * nbar.powi(10) * v.powi(10)
                            - 12. * nbar.powi(11) * v.powi(11)
                            - 1. * nbar.powi(12) * v.powi(12)))
                        / exp(nbar * v)
            }
            14 => {
                (1.6059043836821613e-10
                    * (-6.2270208e9 * nbar + 6.2270208e9 * exp(nbar * v) * nbar
                        - 6.2270208e9 * nbar.powi(2) * v
                        - 3.1135104e9 * nbar * (nbar * v).powi(2)
                        - 1.0378368e9 * nbar * (nbar * v).powi(3)
                        - 2.594592e8 * nbar * (nbar * v).powi(4)
                        - 5.189184e7 * nbar * (nbar * v).powi(5)
                        - 8.64864e6 * nbar * (nbar * v).powi(6)
                        - 1.23552e6 * nbar * (nbar * v).powi(7)
                        - 154440. * nbar * (nbar * v).powi(8)
                        - 17160. * nbar * (nbar * v).powi(9)
                        - 1716. * nbar * (nbar * v).powi(10)
                        - 156. * nbar * (nbar * v).powi(11)
                        - 13. * nbar * (nbar * v).powi(12)))
                    / exp(nbar * v)
                    - (1.6059043836821613e-10
                        * nbar
                        * (-6.2270208e9 + 6.2270208e9 * exp(nbar * v)
                            - 6.2270208e9 * nbar * v
                            - 3.1135104e9 * nbar.powi(2) * v.powi(2)
                            - 1.0378368e9 * nbar.powi(3) * v.powi(3)
                            - 2.594592e8 * nbar.powi(4) * v.powi(4)
                            - 5.189184e7 * nbar.powi(5) * v.powi(5)
                            - 8.64864e6 * nbar.powi(6) * v.powi(6)
                            - 1.23552e6 * nbar.powi(7) * v.powi(7)
                            - 154440. * nbar.powi(8) * v.powi(8)
                            - 17160. * nbar.powi(9) * v.powi(9)
                            - 1716. * nbar.powi(10) * v.powi(10)
                            - 156. * nbar.powi(11) * v.powi(11)
                            - 13. * nbar.powi(12) * v.powi(12)
                            - 1. * nbar.powi(13) * v.powi(13)))
                        / exp(nbar * v)
            }
            15 => {
                (1.1470745597729725e-11
                    * (-8.71782912e10 * nbar + 8.71782912e10 * exp(nbar * v) * nbar
                        - 8.71782912e10 * nbar.powi(2) * v
                        - 4.35891456e10 * nbar * (nbar * v).powi(2)
                        - 1.45297152e10 * nbar * (nbar * v).powi(3)
                        - 3.6324288e9 * nbar * (nbar * v).powi(4)
                        - 7.2648576e8 * nbar * (nbar * v).powi(5)
                        - 1.2108096e8 * nbar * (nbar * v).powi(6)
                        - 1.729728e7 * nbar * (nbar * v).powi(7)
                        - 2.16216e6 * nbar * (nbar * v).powi(8)
                        - 240240. * nbar * (nbar * v).powi(9)
                        - 24024. * nbar * (nbar * v).powi(10)
                        - 2184. * nbar * (nbar * v).powi(11)
                        - 182. * nbar * (nbar * v).powi(12)
                        - 14. * nbar * (nbar * v).powi(13)))
                    / exp(nbar * v)
                    - (1.1470745597729725e-11
                        * nbar
                        * (-8.71782912e10 + 8.71782912e10 * exp(nbar * v)
                            - 8.71782912e10 * nbar * v
                            - 4.35891456e10 * nbar.powi(2) * v.powi(2)
                            - 1.45297152e10 * nbar.powi(3) * v.powi(3)
                            - 3.6324288e9 * nbar.powi(4) * v.powi(4)
                            - 7.2648576e8 * nbar.powi(5) * v.powi(5)
                            - 1.2108096e8 * nbar.powi(6) * v.powi(6)
                            - 1.729728e7 * nbar.powi(7) * v.powi(7)
                            - 2.16216e6 * nbar.powi(8) * v.powi(8)
                            - 240240. * nbar.powi(9) * v.powi(9)
                            - 24024. * nbar.powi(10) * v.powi(10)
                            - 2184. * nbar.powi(11) * v.powi(11)
                            - 182. * nbar.powi(12) * v.powi(12)
                            - 14. * nbar.powi(13) * v.powi(13)
                            - 1. * nbar.powi(14) * v.powi(14)))
                        / exp(nbar * v)
            }
            16 => {
                (7.647163731819816e-13
                    * (-1.307674368e12 * nbar + 1.307674368e12 * exp(nbar * v) * nbar
                        - 1.307674368e12 * nbar.powi(2) * v
                        - 6.53837184e11 * nbar * (nbar * v).powi(2)
                        - 2.17945728e11 * nbar * (nbar * v).powi(3)
                        - 5.4486432e10 * nbar * (nbar * v).powi(4)
                        - 1.08972864e10 * nbar * (nbar * v).powi(5)
                        - 1.8162144e9 * nbar * (nbar * v).powi(6)
                        - 2.594592e8 * nbar * (nbar * v).powi(7)
                        - 3.24324e7 * nbar * (nbar * v).powi(8)
                        - 3.6036e6 * nbar * (nbar * v).powi(9)
                        - 360360. * nbar * (nbar * v).powi(10)
                        - 32760. * nbar * (nbar * v).powi(11)
                        - 2730. * nbar * (nbar * v).powi(12)
                        - 210. * nbar * (nbar * v).powi(13)
                        - 15. * nbar * (nbar * v).powi(14)))
                    / exp(nbar * v)
                    - (7.647163731819816e-13
                        * nbar
                        * (-1.307674368e12 + 1.307674368e12 * exp(nbar * v)
                            - 1.307674368e12 * nbar * v
                            - 6.53837184e11 * nbar.powi(2) * v.powi(2)
                            - 2.17945728e11 * nbar.powi(3) * v.powi(3)
                            - 5.4486432e10 * nbar.powi(4) * v.powi(4)
                            - 1.08972864e10 * nbar.powi(5) * v.powi(5)
                            - 1.8162144e9 * nbar.powi(6) * v.powi(6)
                            - 2.594592e8 * nbar.powi(7) * v.powi(7)
                            - 3.24324e7 * nbar.powi(8) * v.powi(8)
                            - 3.6036e6 * nbar.powi(9) * v.powi(9)
                            - 360360. * nbar.powi(10) * v.powi(10)
                            - 32760. * nbar.powi(11) * v.powi(11)
                            - 2730. * nbar.powi(12) * v.powi(12)
                            - 210. * nbar.powi(13) * v.powi(13)
                            - 15. * nbar.powi(14) * v.powi(14)
                            - 1. * nbar.powi(15) * v.powi(15)))
                        / exp(nbar * v)
            }

            _ => todo!(
                "only k=1 through k=16 NN PDFs are implemented for the Gaussian Random Fields"
            ),
        }
}

/// This is just here to simplify the conversion from mathematica output
/// as simply as possible, i.e. without any regex (or similar) manipulation.
#[inline(always)]
#[allow(unused)]
fn exp(x: f64) -> f64 {
    x.exp()
}
