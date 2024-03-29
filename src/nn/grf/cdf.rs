pub(super) fn calculate_grf_cdf(r: f64, nbar: f64, i: f64, k: u8) -> f64 {
    // Convert distance to corresponding spherical volume
    let v = 4.0 * std::f64::consts::PI * r.powi(3) / 3.0;

    match k {
        1 => 1. - exp(0.5 * i * nbar.powi(2) - nbar * v),
        2 => {
            1. - exp(0.5 * i * nbar.powi(2) - nbar * v)
                - exp(0.5 * i * nbar.powi(2) - nbar * v) * (-i * nbar.powi(2) + nbar * v)
        }
        3 => {
            0.5 * (2. * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                - exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v) * (-i * nbar.powi(2) + nbar * v)
                - exp(0.5 * i * nbar.powi(2) - nbar * v) * (-i * nbar.powi(2) + nbar * v).powi(2))
        }
        4 => {
            0.16666666666666666
                * (6. * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v) * (-i * nbar.powi(2) + nbar * v)
                    - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(3)
                    + 3. * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)))
        }
        5 => {
            0.041666666666666664
                * (24. * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(2) * nbar.powi(4)
                    - 24. * exp(0.5 * i * nbar.powi(2) - nbar * v) * (-i * nbar.powi(2) + nbar * v)
                    - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(2)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(4)
                    + 12.
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 4. * (-3.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)))
        }
        6 => {
            0.008333333333333333
                * (120. * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 120.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 15.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 10.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(3)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(5)
                    + 60.
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 20.
                        * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3))
                    + 5. * (-3.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)))
        }
        7 => {
            0.001388888888888889
                * (720. * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 15. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(3) * nbar.powi(6)
                    - 720.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 45.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v).powi(2)
                    - 15.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(4)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(6)
                    + 360.
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 120.
                        * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3))
                    + 30.
                        * (-3. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(2) * nbar.powi(4)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4))
                    + 6. * (-15.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v)
                        - 10.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)))
        }
        8 => {
            0.0001984126984126984
                * (5040. * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 5040.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 105.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(6)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 105.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v).powi(3)
                    - 21.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(5)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(7)
                    + 2520.
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 840.
                        * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3))
                    + 210.
                        * (-3. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(2) * nbar.powi(4)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4))
                    + 42.
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 10.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(5))
                    + 7. * (-15.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(6)
                        - 45.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                        - 15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)))
        }
        9 => {
            0.0000248015873015873
                * (40320. * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 105. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(4) * nbar.powi(8)
                    - 40320.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 420.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(6)
                        * (-i * nbar.powi(2) + nbar * v).powi(2)
                    - 210.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v).powi(4)
                    - 28.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(6)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(8)
                    + 20160.
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 6720.
                        * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3))
                    + 1680.
                        * (-3. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(2) * nbar.powi(4)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4))
                    + 336.
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 10.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(5))
                    + 56.
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(6))
                    + 8. * (-105.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(6)
                        * (-i * nbar.powi(2) + nbar * v)
                        - 105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                        - 21.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)))
        }
        10 => {
            2.7557319223985893e-6
                * (362880. * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 362880.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 945.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(4)
                        * nbar.powi(8)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 1260.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(6)
                        * (-i * nbar.powi(2) + nbar * v).powi(3)
                    - 378.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v).powi(5)
                    - 36.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(7)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(9)
                    + 181440.
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 60480.
                        * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3))
                    + 15120.
                        * (-3. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(2) * nbar.powi(4)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4))
                    + 3024.
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 10.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(5))
                    + 504.
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(6))
                    + 72.
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 21.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(7))
                    + 9. * (-105.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(4)
                        * nbar.powi(8)
                        - 420.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                        - 210.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                        - 28.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)))
        }
        11 => {
            2.755731922398589e-7
                * (3.6288e6 * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 945. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(5) * nbar.powi(10)
                    - 3.6288e6
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 4725.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(4)
                        * nbar.powi(8)
                        * (-i * nbar.powi(2) + nbar * v).powi(2)
                    - 3150.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(6)
                        * (-i * nbar.powi(2) + nbar * v).powi(4)
                    - 630.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v).powi(6)
                    - 45.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(8)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(10)
                    + 1.8144e6
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 604800.
                        * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3))
                    + 151200.
                        * (-3. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(2) * nbar.powi(4)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4))
                    + 30240.
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 10.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(5))
                    + 5040.
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(6))
                    + 720.
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 21.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(7))
                    + 90.
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            - 420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 210.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 28.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(8))
                    + 10.
                        * (-945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 1260.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 378.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 36.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(9)))
        }
        12 => {
            2.505210838544172e-8
                * (3.99168e7 * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 3.99168e7
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 10395.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(5)
                        * nbar.powi(10)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 17325.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(4)
                        * nbar.powi(8)
                        * (-i * nbar.powi(2) + nbar * v).powi(3)
                    - 6930.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(6)
                        * (-i * nbar.powi(2) + nbar * v).powi(5)
                    - 990.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v).powi(7)
                    - 55.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(9)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(11)
                    + 1.99584e7
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 6.6528e6
                        * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3))
                    + 1.6632e6
                        * (-3. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(2) * nbar.powi(4)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4))
                    + 332640.
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 10.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(5))
                    + 55440.
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(6))
                    + 7920.
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 21.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(7))
                    + 990.
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            - 420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 210.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 28.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(8))
                    + 110.
                        * (-945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 1260.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 378.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 36.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(9))
                    + 11.
                        * (-945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            - 4725.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 3150.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 630.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(8)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(10)))
        }
        13 => {
            2.08767569878681e-9
                * (4.790016e8 * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 10395. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(6) * nbar.powi(12)
                    - 4.790016e8
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 62370.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(5)
                        * nbar.powi(10)
                        * (-i * nbar.powi(2) + nbar * v).powi(2)
                    - 51975.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(4)
                        * nbar.powi(8)
                        * (-i * nbar.powi(2) + nbar * v).powi(4)
                    - 13860.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(6)
                        * (-i * nbar.powi(2) + nbar * v).powi(6)
                    - 1485.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v).powi(8)
                    - 66.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(10)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(12)
                    + 2.395008e8
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 7.98336e7
                        * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3))
                    + 1.99584e7
                        * (-3. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(2) * nbar.powi(4)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4))
                    + 3.99168e6
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 10.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(5))
                    + 665280.
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(6))
                    + 95040.
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 21.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(7))
                    + 11880.
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            - 420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 210.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 28.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(8))
                    + 1320.
                        * (-945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 1260.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 378.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 36.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(9))
                    + 132.
                        * (-945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            - 4725.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 3150.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 630.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(8)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(10))
                    + 12.
                        * (-10395.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 17325.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 6930.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 990.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - 55.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(9)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(11)))
        }
        14 => {
            1.6059043836821613e-10
                * (6.2270208e9 * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 6.2270208e9
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 135135.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(6)
                        * nbar.powi(12)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 270270.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(5)
                        * nbar.powi(10)
                        * (-i * nbar.powi(2) + nbar * v).powi(3)
                    - 135135.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(4)
                        * nbar.powi(8)
                        * (-i * nbar.powi(2) + nbar * v).powi(5)
                    - 25740.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(6)
                        * (-i * nbar.powi(2) + nbar * v).powi(7)
                    - 2145.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v).powi(9)
                    - 78.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(11)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(13)
                    + 3.1135104e9
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 1.0378368e9
                        * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3))
                    + 2.594592e8
                        * (-3. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(2) * nbar.powi(4)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4))
                    + 5.189184e7
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 10.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(5))
                    + 8.64864e6
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(6))
                    + 1.23552e6
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 21.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(7))
                    + 154440.
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            - 420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 210.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 28.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(8))
                    + 17160.
                        * (-945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 1260.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 378.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 36.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(9))
                    + 1716.
                        * (-945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            - 4725.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 3150.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 630.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(8)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(10))
                    + 156.
                        * (-10395.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 17325.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 6930.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 990.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - 55.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(9)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(11))
                    + 13.
                        * (-10395.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            - 62370.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(5)
                                * nbar.powi(10)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 51975.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 13860.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - 1485.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(8)
                            - 66.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(10)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(12)))
        }
        15 => {
            1.1470745597729725e-11
                * (8.71782912e10 * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 135135. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(7) * nbar.powi(14)
                    - 8.71782912e10
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 945945.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(6)
                        * nbar.powi(12)
                        * (-i * nbar.powi(2) + nbar * v).powi(2)
                    - 945945.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(5)
                        * nbar.powi(10)
                        * (-i * nbar.powi(2) + nbar * v).powi(4)
                    - 315315.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(4)
                        * nbar.powi(8)
                        * (-i * nbar.powi(2) + nbar * v).powi(6)
                    - 45045.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(6)
                        * (-i * nbar.powi(2) + nbar * v).powi(8)
                    - 3003.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v).powi(10)
                    - 91.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(12)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(14)
                    + 4.35891456e10
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 1.45297152e10
                        * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3))
                    + 3.6324288e9
                        * (-3. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(2) * nbar.powi(4)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4))
                    + 7.2648576e8
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 10.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(5))
                    + 1.2108096e8
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(6))
                    + 1.729728e7
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 21.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(7))
                    + 2.16216e6
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            - 420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 210.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 28.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(8))
                    + 240240.
                        * (-945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 1260.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 378.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 36.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(9))
                    + 24024.
                        * (-945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            - 4725.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 3150.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 630.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(8)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(10))
                    + 2184.
                        * (-10395.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 17325.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 6930.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 990.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - 55.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(9)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(11))
                    + 182.
                        * (-10395.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            - 62370.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(5)
                                * nbar.powi(10)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 51975.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 13860.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - 1485.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(8)
                            - 66.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(10)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(12))
                    + 14.
                        * (-135135.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 270270.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(5)
                                * nbar.powi(10)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 135135.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 25740.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - 2145.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(9)
                            - 78.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(11)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(13)))
        }
        16 => {
            7.647163731819816e-13
                * (1.307674368e12 * (1. - exp(0.5 * i * nbar.powi(2) - nbar * v))
                    - 1.307674368e12
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 2.027025e6
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(7)
                        * nbar.powi(14)
                        * (-i * nbar.powi(2) + nbar * v)
                    - 4.729725e6
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(6)
                        * nbar.powi(12)
                        * (-i * nbar.powi(2) + nbar * v).powi(3)
                    - 2.837835e6
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(5)
                        * nbar.powi(10)
                        * (-i * nbar.powi(2) + nbar * v).powi(5)
                    - 675675.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(4)
                        * nbar.powi(8)
                        * (-i * nbar.powi(2) + nbar * v).powi(7)
                    - 75075.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(6)
                        * (-i * nbar.powi(2) + nbar * v).powi(9)
                    - 4095.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v).powi(11)
                    - 105.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v).powi(13)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(15)
                    + 6.53837184e11
                        * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2))
                    + 2.17945728e11
                        * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3))
                    + 5.4486432e10
                        * (-3. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i.powi(2) * nbar.powi(4)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4))
                    + 1.08972864e10
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 10.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(5))
                    + 1.8162144e9
                        * (-15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(6))
                    + 2.594592e8
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 21.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(7))
                    + 3.24324e7
                        * (-105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            - 420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 210.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 28.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(8))
                    + 3.6036e6
                        * (-945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 1260.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 378.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 36.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(9))
                    + 360360.
                        * (-945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            - 4725.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 3150.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 630.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(8)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(10))
                    + 32760.
                        * (-10395.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 17325.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 6930.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 990.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - 55.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(9)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(11))
                    + 2730.
                        * (-10395.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            - 62370.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(5)
                                * nbar.powi(10)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 51975.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 13860.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - 1485.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(8)
                            - 66.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(10)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(12))
                    + 210.
                        * (-135135.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            * (-i * nbar.powi(2) + nbar * v)
                            - 270270.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(5)
                                * nbar.powi(10)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                            - 135135.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                            - 25740.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                            - 2145.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(9)
                            - 78.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(11)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(13))
                    + 15.
                        * (-135135.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(7)
                            * nbar.powi(14)
                            - 945945.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(6)
                                * nbar.powi(12)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                            - 945945.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(5)
                                * nbar.powi(10)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                            - 315315.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                            - 45045.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(8)
                            - 3003.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(10)
                            - 91.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(12)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(14)))
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
