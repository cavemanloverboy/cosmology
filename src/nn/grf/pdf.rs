use std::f64::consts::PI;

pub(super) fn calculate_grf_pdf(r: f64, nbar: f64, i: f64, didv: f64, k: u8) -> f64 {
    // Convert distance to corresponding spherical volume
    let v = 4.0 * std::f64::consts::PI * r.powi(3) / 3.0;

    // int P(V) dV = 1
    // int P(V(R)) * dR dV/dR = 1
    // ----> P(R) = dV/dR * P(V)
    // ----> P(R) = 4 * PI * R^2 * dCDF(V)/dV
    let jacobian = 4.0 * PI * r.powi(2);

    jacobian
        * match k {
            1 => -exp(0.5 * i * nbar.powi(2) - nbar * v) * (-nbar + 0.5 * nbar.powi(2) * didv),
            2 => {
                -exp(0.5 * i * nbar.powi(2) - nbar * v) * (nbar - nbar.powi(2) * didv)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v) * (-nbar + 0.5 * nbar.powi(2) * didv)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                        * (-nbar + 0.5 * nbar.powi(2) * didv)
            }
            3 => {
                0.5 * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                    - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v) * (nbar - nbar.powi(2) * didv)
                    - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                        * (nbar - nbar.powi(2) * didv)
                    - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-nbar + 0.5 * nbar.powi(2) * didv)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(2)
                        * (-nbar + 0.5 * nbar.powi(2) * didv)
                    - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v)
                        * (-nbar + 0.5 * nbar.powi(2) * didv)
                    - exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * (-i * nbar.powi(2) + nbar * v).powi(2)
                        * (-nbar + 0.5 * nbar.powi(2) * didv))
            }
            4 => {
                0.16666666666666666
                    * (-3.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * nbar.powi(2)
                        * (-i * nbar.powi(2) + nbar * v)
                        * didv
                        - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (nbar - nbar.powi(2) * didv)
                        - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (nbar - nbar.powi(2) * didv)
                        - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 3. * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                            - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v)
                                * (nbar - nbar.powi(2) * didv)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            5 => {
                0.041666666666666664
                    * (-6. * exp(0.5 * i * nbar.powi(2) - nbar * v) * i * nbar.powi(4) * didv
                        - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * didv
                        - 24.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 12.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (nbar - nbar.powi(2) * didv)
                        - 24.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 24.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 12.
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 4. * (-3.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v)
                            * didv
                            - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (nbar - nbar.powi(2) * didv)
                            - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                                * (nbar - nbar.powi(2) * didv)
                            - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            6 => {
                0.008333333333333333
                    * (-30.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i
                        * nbar.powi(4)
                        * (-i * nbar.powi(2) + nbar * v)
                        * didv
                        - 10.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * didv
                        - 120.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (nbar - nbar.powi(2) * didv)
                        - 30.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (nbar - nbar.powi(2) * didv)
                        - 5. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (nbar - nbar.powi(2) * didv)
                        - 120.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 120.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 10.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 60.
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 20.
                            * (-3.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 5. * (-6.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * didv
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                                * didv
                            - 12.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * (nbar - nbar.powi(2) * didv)
                            - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                                * (nbar - nbar.powi(2) * didv)
                            - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            7 => {
                0.001388888888888889
                    * (-45.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(6)
                        * didv
                        - 90.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * didv
                        - 15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * didv
                        - 720.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 90.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 60.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (nbar - nbar.powi(2) * didv)
                        - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * (nbar - nbar.powi(2) * didv)
                        - 720.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 720.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 45.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 360.
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 120.
                            * (-3.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 30.
                            * (-6.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * didv
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 6. * (-30.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v)
                            * didv
                            - 10.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                                * didv
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (nbar - nbar.powi(2) * didv)
                            - 30.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                                * (nbar - nbar.powi(2) * didv)
                            - 5. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                                * (nbar - nbar.powi(2) * didv)
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - 10.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            8 => {
                0.0001984126984126984
                    * (-315.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(2)
                        * nbar.powi(6)
                        * (-i * nbar.powi(2) + nbar * v)
                        * didv
                        - 210.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * didv
                        - 21.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * didv
                        - 5040.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (nbar - nbar.powi(2) * didv)
                        - 315.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (nbar - nbar.powi(2) * didv)
                        - 105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (nbar - nbar.powi(2) * didv)
                        - 7. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * (nbar - nbar.powi(2) * didv)
                        - 5040.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 5040.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 21.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 2520.
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 840.
                            * (-3.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 210.
                            * (-6.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * didv
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 42.
                            * (-30.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 30.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 5. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 7. * (-45.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(6)
                            * didv
                            - 90.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                                * didv
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                                * didv
                            - 90.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v)
                                * (nbar - nbar.powi(2) * didv)
                            - 60.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                                * (nbar - nbar.powi(2) * didv)
                            - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                                * (nbar - nbar.powi(2) * didv)
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - 45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - 15.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            9 => {
                0.0000248015873015873
                    * (-420.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(8)
                        * didv
                        - 1260.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * didv
                        - 420.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * didv
                        - 28.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * didv
                        - 40320.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 840.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 840.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (nbar - nbar.powi(2) * didv)
                        - 168.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * (nbar - nbar.powi(2) * didv)
                        - 8. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * (nbar - nbar.powi(2) * didv)
                        - 40320.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 40320.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 420.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 210.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 28.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 20160.
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 6720.
                            * (-3.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1680.
                            * (-6.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * didv
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 336.
                            * (-30.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 30.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 5. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 56.
                            * (-45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 60.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 8. * (-315.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v)
                            * didv
                            - 210.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                                * didv
                            - 21.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                                * didv
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (nbar - nbar.powi(2) * didv)
                            - 315.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                                * (nbar - nbar.powi(2) * didv)
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                                * (nbar - nbar.powi(2) * didv)
                            - 7. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                                * (nbar - nbar.powi(2) * didv)
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - 21.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            10 => {
                2.7557319223985893e-6
                    * (-3780.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(3)
                        * nbar.powi(8)
                        * (-i * nbar.powi(2) + nbar * v)
                        * didv
                        - 3780.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * didv
                        - 756.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * didv
                        - 36.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * didv
                        - 362880.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (nbar - nbar.powi(2) * didv)
                        - 3780.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (nbar - nbar.powi(2) * didv)
                        - 1890.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (nbar - nbar.powi(2) * didv)
                        - 252.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * (nbar - nbar.powi(2) * didv)
                        - 9. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)
                            * (nbar - nbar.powi(2) * didv)
                        - 362880.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 362880.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 1260.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 378.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 36.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(9)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 181440.
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 60480.
                            * (-3.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 15120.
                            * (-6.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * didv
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 3024.
                            * (-30.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 30.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 5. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 504.
                            * (-45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 60.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 72.
                            * (-315.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 315.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 7. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 9. * (-420.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(8)
                            * didv
                            - 1260.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                                * didv
                            - 420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                                * didv
                            - 28.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                                * didv
                            - 840.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v)
                                * (nbar - nbar.powi(2) * didv)
                            - 840.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(3)
                                * (nbar - nbar.powi(2) * didv)
                            - 168.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(5)
                                * (nbar - nbar.powi(2) * didv)
                            - 8. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(7)
                                * (nbar - nbar.powi(2) * didv)
                            - 105.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(8)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - 420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v).powi(2)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - 210.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v).powi(4)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - 28.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v).powi(6)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)
                            - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * (-i * nbar.powi(2) + nbar * v).powi(8)
                                * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            11 => {
                2.755731922398589e-7
                    * (-4725.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(4)
                        * nbar.powi(10)
                        * didv
                        - 18900.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * didv
                        - 9450.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * didv
                        - 1260.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * didv
                        - 45.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)
                            * didv
                        - 3.6288e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 9450.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 12600.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (nbar - nbar.powi(2) * didv)
                        - 3780.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * (nbar - nbar.powi(2) * didv)
                        - 360.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * (nbar - nbar.powi(2) * didv)
                        - 10.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(9)
                            * (nbar - nbar.powi(2) * didv)
                        - 3.6288e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 3.6288e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 4725.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 3150.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 630.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 45.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(10)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 1.8144e6
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 604800.
                            * (-3.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 151200.
                            * (-6.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * didv
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 30240.
                            * (-30.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 30.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 5. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 5040.
                            * (-45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 60.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 720.
                            * (-315.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 315.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 7. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 90.
                            * (-420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * didv
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 168.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 8. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 10.
                            * (-3780.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 756.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 1890.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 252.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 9. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 378.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            12 => {
                2.505210838544172e-8
                    * (-51975.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(4)
                        * nbar.powi(10)
                        * (-i * nbar.powi(2) + nbar * v)
                        * didv
                        - 69300.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * didv
                        - 20790.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * didv
                        - 1980.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * didv
                        - 55.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(9)
                            * didv
                        - 3.99168e7
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 10395.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (nbar - nbar.powi(2) * didv)
                        - 51975.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (nbar - nbar.powi(2) * didv)
                        - 34650.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (nbar - nbar.powi(2) * didv)
                        - 6930.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * (nbar - nbar.powi(2) * didv)
                        - 495.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)
                            * (nbar - nbar.powi(2) * didv)
                        - 11.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(10)
                            * (nbar - nbar.powi(2) * didv)
                        - 3.99168e7
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 3.99168e7
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 10395.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 17325.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 6930.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 990.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 55.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(9)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(11)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 1.99584e7
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 6.6528e6
                            * (-3.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1.6632e6
                            * (-6.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * didv
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 332640.
                            * (-30.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 30.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 5. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 55440.
                            * (-45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 60.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 7920.
                            * (-315.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 315.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 7. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 990.
                            * (-420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * didv
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 168.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 8. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 110.
                            * (-3780.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 756.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 1890.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 252.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 9. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 378.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 11.
                            * (-4725.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(10)
                                * didv
                                - 18900.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 9450.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * didv
                                - 9450.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 12600.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 360.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (nbar - nbar.powi(2) * didv)
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 4725.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 3150.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 630.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            13 => {
                2.08767569878681e-9
                    * (-62370.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(5)
                        * nbar.powi(12)
                        * didv
                        - 311850.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * didv
                        - 207900.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * didv
                        - 41580.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * didv
                        - 2970.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)
                            * didv
                        - 66.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(10)
                            * didv
                        - 4.790016e8
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 124740.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 207900.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (nbar - nbar.powi(2) * didv)
                        - 83160.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * (nbar - nbar.powi(2) * didv)
                        - 11880.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * (nbar - nbar.powi(2) * didv)
                        - 660.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(9)
                            * (nbar - nbar.powi(2) * didv)
                        - 12.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(11)
                            * (nbar - nbar.powi(2) * didv)
                        - 4.790016e8
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 10395.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 4.790016e8
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 62370.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 51975.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 13860.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 1485.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 66.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(10)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(12)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 2.395008e8
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 7.98336e7
                            * (-3.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1.99584e7
                            * (-6.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * didv
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 3.99168e6
                            * (-30.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 30.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 5. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 665280.
                            * (-45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 60.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 95040.
                            * (-315.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 315.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 7. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 11880.
                            * (-420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * didv
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 168.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 8. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1320.
                            * (-3780.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 756.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 1890.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 252.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 9. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 378.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 132.
                            * (-4725.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(10)
                                * didv
                                - 18900.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 9450.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * didv
                                - 9450.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 12600.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 360.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (nbar - nbar.powi(2) * didv)
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 4725.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 3150.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 630.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 12.
                            * (-51975.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(10)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 69300.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 20790.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 1980.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 55.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * didv
                                - 10395.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (nbar - nbar.powi(2) * didv)
                                - 51975.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 34650.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6930.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 495.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 11.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10395.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 17325.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6930.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 990.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 55.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            14 => {
                1.6059043836821613e-10
                    * (-810810.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(5)
                        * nbar.powi(12)
                        * (-i * nbar.powi(2) + nbar * v)
                        * didv
                        - 1.35135e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * didv
                        - 540540.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * didv
                        - 77220.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * didv
                        - 4290.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(9)
                            * didv
                        - 78.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(11)
                            * didv
                        - 6.2270208e9
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 135135.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            * (nbar - nbar.powi(2) * didv)
                        - 810810.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (nbar - nbar.powi(2) * didv)
                        - 675675.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (nbar - nbar.powi(2) * didv)
                        - 180180.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * (nbar - nbar.powi(2) * didv)
                        - 19305.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)
                            * (nbar - nbar.powi(2) * didv)
                        - 858.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(10)
                            * (nbar - nbar.powi(2) * didv)
                        - 13.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(12)
                            * (nbar - nbar.powi(2) * didv)
                        - 6.2270208e9
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 6.2270208e9
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 135135.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 270270.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 135135.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 25740.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 2145.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(9)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 78.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(11)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(13)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 3.1135104e9
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1.0378368e9
                            * (-3.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 2.594592e8
                            * (-6.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * didv
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 5.189184e7
                            * (-30.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 30.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 5. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 8.64864e6
                            * (-45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 60.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1.23552e6
                            * (-315.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 315.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 7. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 154440.
                            * (-420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * didv
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 168.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 8. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 17160.
                            * (-3780.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 756.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 1890.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 252.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 9. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 378.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1716.
                            * (-4725.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(10)
                                * didv
                                - 18900.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 9450.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * didv
                                - 9450.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 12600.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 360.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (nbar - nbar.powi(2) * didv)
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 4725.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 3150.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 630.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 156.
                            * (-51975.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(10)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 69300.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 20790.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 1980.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 55.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * didv
                                - 10395.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (nbar - nbar.powi(2) * didv)
                                - 51975.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 34650.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6930.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 495.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 11.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10395.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 17325.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6930.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 990.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 55.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 13.
                            * (-62370.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(5)
                                * nbar.powi(12)
                                * didv
                                - 311850.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 207900.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 41580.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 2970.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * didv
                                - 66.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * didv
                                - 124740.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 207900.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 83160.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 11880.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 660.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (nbar - nbar.powi(2) * didv)
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10395.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(6)
                                    * nbar.powi(12)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 62370.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 51975.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 13860.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 1485.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 66.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(12)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            15 => {
                1.1470745597729725e-11
                    * (-945945.
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(6)
                        * nbar.powi(14)
                        * didv
                        - 5.67567e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(12)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * didv
                        - 4.729725e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * didv
                        - 1.26126e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * didv
                        - 135135.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)
                            * didv
                        - 6006.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(10)
                            * didv
                        - 91.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(12)
                            * didv
                        - 8.71782912e10
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 1.89189e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 3.78378e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (nbar - nbar.powi(2) * didv)
                        - 1.89189e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * (nbar - nbar.powi(2) * didv)
                        - 360360.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * (nbar - nbar.powi(2) * didv)
                        - 30030.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(9)
                            * (nbar - nbar.powi(2) * didv)
                        - 1092.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(11)
                            * (nbar - nbar.powi(2) * didv)
                        - 14.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(13)
                            * (nbar - nbar.powi(2) * didv)
                        - 8.71782912e10
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 135135.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(7)
                            * nbar.powi(14)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 8.71782912e10
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 945945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 945945.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 315315.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 45045.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 3003.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(10)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 91.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(12)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(14)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 4.35891456e10
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1.45297152e10
                            * (-3.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 3.6324288e9
                            * (-6.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * didv
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 7.2648576e8
                            * (-30.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 30.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 5. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1.2108096e8
                            * (-45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 60.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1.729728e7
                            * (-315.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 315.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 7. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 2.16216e6
                            * (-420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * didv
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 168.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 8. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 240240.
                            * (-3780.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 756.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 1890.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 252.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 9. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 378.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 24024.
                            * (-4725.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(10)
                                * didv
                                - 18900.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 9450.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * didv
                                - 9450.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 12600.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 360.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (nbar - nbar.powi(2) * didv)
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 4725.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 3150.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 630.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 2184.
                            * (-51975.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(10)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 69300.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 20790.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 1980.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 55.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * didv
                                - 10395.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (nbar - nbar.powi(2) * didv)
                                - 51975.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 34650.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6930.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 495.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 11.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10395.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 17325.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6930.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 990.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 55.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 182.
                            * (-62370.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(5)
                                * nbar.powi(12)
                                * didv
                                - 311850.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 207900.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 41580.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 2970.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * didv
                                - 66.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * didv
                                - 124740.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 207900.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 83160.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 11880.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 660.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (nbar - nbar.powi(2) * didv)
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10395.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(6)
                                    * nbar.powi(12)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 62370.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 51975.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 13860.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 1485.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 66.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(12)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 14.
                            * (-810810.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(5)
                                * nbar.powi(12)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 1.35135e6
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 540540.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 77220.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 4290.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * didv
                                - 78.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * didv
                                - 135135.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(6)
                                    * nbar.powi(12)
                                    * (nbar - nbar.powi(2) * didv)
                                - 810810.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 675675.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 180180.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 19305.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 858.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (nbar - nbar.powi(2) * didv)
                                - 13.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(12)
                                    * (nbar - nbar.powi(2) * didv)
                                - 135135.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(6)
                                    * nbar.powi(12)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 270270.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 135135.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 25740.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 2145.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 78.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(13)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }
            16 => {
                7.647163731819816e-13
                    * (-1.4189175e7
                        * exp(0.5 * i * nbar.powi(2) - nbar * v)
                        * i.powi(6)
                        * nbar.powi(14)
                        * (-i * nbar.powi(2) + nbar * v)
                        * didv
                        - 2.837835e7
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(12)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * didv
                        - 1.4189175e7
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * didv
                        - 2.7027e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * didv
                        - 225225.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(9)
                            * didv
                        - 8190.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(11)
                            * didv
                        - 105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(13)
                            * didv
                        - 1.307674368e12
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (nbar - nbar.powi(2) * didv)
                        - 2.027025e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(7)
                            * nbar.powi(14)
                            * (nbar - nbar.powi(2) * didv)
                        - 1.4189175e7
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            * (-i * nbar.powi(2) + nbar * v).powi(2)
                            * (nbar - nbar.powi(2) * didv)
                        - 1.4189175e7
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v).powi(4)
                            * (nbar - nbar.powi(2) * didv)
                        - 4.729725e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(6)
                            * (nbar - nbar.powi(2) * didv)
                        - 675675.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(8)
                            * (nbar - nbar.powi(2) * didv)
                        - 45045.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(10)
                            * (nbar - nbar.powi(2) * didv)
                        - 1365.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(12)
                            * (nbar - nbar.powi(2) * didv)
                        - 15.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(14)
                            * (nbar - nbar.powi(2) * didv)
                        - 1.307674368e12
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 1.307674368e12
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 2.027025e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(7)
                            * nbar.powi(14)
                            * (-i * nbar.powi(2) + nbar * v)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 4.729725e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(6)
                            * nbar.powi(12)
                            * (-i * nbar.powi(2) + nbar * v).powi(3)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 2.837835e6
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(5)
                            * nbar.powi(10)
                            * (-i * nbar.powi(2) + nbar * v).powi(5)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 675675.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(4)
                            * nbar.powi(8)
                            * (-i * nbar.powi(2) + nbar * v).powi(7)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 75075.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(3)
                            * nbar.powi(6)
                            * (-i * nbar.powi(2) + nbar * v).powi(9)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 4095.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i.powi(2)
                            * nbar.powi(4)
                            * (-i * nbar.powi(2) + nbar * v).powi(11)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - 105.
                            * exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * i
                            * nbar.powi(2)
                            * (-i * nbar.powi(2) + nbar * v).powi(13)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        - exp(0.5 * i * nbar.powi(2) - nbar * v)
                            * (-i * nbar.powi(2) + nbar * v).powi(15)
                            * (-nbar + 0.5 * nbar.powi(2) * didv)
                        + 6.53837184e11
                            * (-exp(0.5 * i * nbar.powi(2) - nbar * v) * nbar.powi(2) * didv
                                - 2. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 2.17945728e11
                            * (-3.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * nbar.powi(2)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 5.4486432e10
                            * (-6.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * didv
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 4. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1.08972864e10
                            * (-30.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i
                                * nbar.powi(4)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 30.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 5. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 1.8162144e9
                            * (-45.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 90.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 60.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 15.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 2.594592e8
                            * (-315.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(2)
                                * nbar.powi(6)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 315.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 7. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 21.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 3.24324e7
                            * (-420.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * didv
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 840.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 168.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 8. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 105.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 420.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 210.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 28.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 3.6036e6
                            * (-3780.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(3)
                                * nbar.powi(8)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 756.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 1890.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 252.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 9. * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 378.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 36.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 360360.
                            * (-4725.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(10)
                                * didv
                                - 18900.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 9450.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 1260.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * didv
                                - 9450.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 12600.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3780.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 360.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (nbar - nbar.powi(2) * didv)
                                - 945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 4725.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 3150.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 630.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 32760.
                            * (-51975.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(4)
                                * nbar.powi(10)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 69300.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 20790.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 1980.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 55.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * didv
                                - 10395.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (nbar - nbar.powi(2) * didv)
                                - 51975.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 34650.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 6930.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 495.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 11.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10395.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 17325.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 6930.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 990.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 55.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 2730.
                            * (-62370.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(5)
                                * nbar.powi(12)
                                * didv
                                - 311850.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 207900.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 41580.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 2970.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * didv
                                - 66.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * didv
                                - 124740.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 207900.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 83160.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 11880.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 660.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (nbar - nbar.powi(2) * didv)
                                - 12.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * (nbar - nbar.powi(2) * didv)
                                - 10395.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(6)
                                    * nbar.powi(12)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 62370.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 51975.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 13860.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 1485.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 66.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(12)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 210.
                            * (-810810.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(5)
                                * nbar.powi(12)
                                * (-i * nbar.powi(2) + nbar * v)
                                * didv
                                - 1.35135e6
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * didv
                                - 540540.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * didv
                                - 77220.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * didv
                                - 4290.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * didv
                                - 78.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * didv
                                - 135135.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(6)
                                    * nbar.powi(12)
                                    * (nbar - nbar.powi(2) * didv)
                                - 810810.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (nbar - nbar.powi(2) * didv)
                                - 675675.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (nbar - nbar.powi(2) * didv)
                                - 180180.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (nbar - nbar.powi(2) * didv)
                                - 19305.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (nbar - nbar.powi(2) * didv)
                                - 858.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (nbar - nbar.powi(2) * didv)
                                - 13.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(12)
                                    * (nbar - nbar.powi(2) * didv)
                                - 135135.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(6)
                                    * nbar.powi(12)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 270270.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 135135.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 25740.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 2145.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 78.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(13)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv))
                        + 15.
                            * (-945945.
                                * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                * i.powi(6)
                                * nbar.powi(14)
                                * didv
                                - 5.67567e6
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(12)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * didv
                                - 4.729725e6
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * didv
                                - 1.26126e6
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * didv
                                - 135135.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * didv
                                - 6006.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * didv
                                - 91.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(12)
                                    * didv
                                - 1.89189e6
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(6)
                                    * nbar.powi(12)
                                    * (-i * nbar.powi(2) + nbar * v)
                                    * (nbar - nbar.powi(2) * didv)
                                - 3.78378e6
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(3)
                                    * (nbar - nbar.powi(2) * didv)
                                - 1.89189e6
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(5)
                                    * (nbar - nbar.powi(2) * didv)
                                - 360360.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(7)
                                    * (nbar - nbar.powi(2) * didv)
                                - 30030.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(9)
                                    * (nbar - nbar.powi(2) * didv)
                                - 1092.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(11)
                                    * (nbar - nbar.powi(2) * didv)
                                - 14.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(13)
                                    * (nbar - nbar.powi(2) * didv)
                                - 135135.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(7)
                                    * nbar.powi(14)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 945945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(6)
                                    * nbar.powi(12)
                                    * (-i * nbar.powi(2) + nbar * v).powi(2)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 945945.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(5)
                                    * nbar.powi(10)
                                    * (-i * nbar.powi(2) + nbar * v).powi(4)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 315315.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(4)
                                    * nbar.powi(8)
                                    * (-i * nbar.powi(2) + nbar * v).powi(6)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 45045.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(3)
                                    * nbar.powi(6)
                                    * (-i * nbar.powi(2) + nbar * v).powi(8)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 3003.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i.powi(2)
                                    * nbar.powi(4)
                                    * (-i * nbar.powi(2) + nbar * v).powi(10)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - 91.
                                    * exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * i
                                    * nbar.powi(2)
                                    * (-i * nbar.powi(2) + nbar * v).powi(12)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)
                                - exp(0.5 * i * nbar.powi(2) - nbar * v)
                                    * (-i * nbar.powi(2) + nbar * v).powi(14)
                                    * (-nbar + 0.5 * nbar.powi(2) * didv)))
            }

            _ => todo!(
                "only k=1 through k=16 NN PDFs are implemented for the Gaussian Random Fields"
            ),
        }
        .max(f64::MIN_POSITIVE)
}

/// This is just here to simplify the conversion from mathematica output
/// as simply as possible, i.e. without any regex (or similar) manipulation.
#[inline(always)]
fn exp(x: f64) -> f64 {
    x.exp()
}
