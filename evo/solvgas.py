# ------------------------------------------------------------------------
# solvgas.py

import constants as cnst
from scipy.optimize import newton
from math import exp
import numpy as np


# ------------------------------------------------------------------------
# FUGACITY COEFFICIENTS
# ------------------------------------------------------------------------


def find_Y(P, T, species_list):
    """
    Calculates the fugacity coefficients for all relevant species.

    Parameters
    ----------
    P : float
        Pressure (bar)
    T : float
        Temperature (K)
    species_list : [strings]
        A list of each molecule in the system given as capitalised strings.

    Returns
    -------
    tuple
        a tuple of all the fugacity coefficients, as floats.

    References
    ----------
    Shaw, H.R., & Wones, D.R. (1964). Fugacity coefficients for hydrogen
    gas between 0 degrees and 1000C, for pressures to 3000 atm.
    American Journal of Science.

    Holland, T., & Powell, R. (1991). A Compensated-Redlich-Kwong (CORK)
    equation for volumes and fugacities of CO2 and H2O in the range 1 bar
    to 50 kbar and 100-1600C. Contributions to Mineralogy and Petrology.

    Shi, P., & Saxena, S. K. (1992). Thermodynamic modeling of the C-H-O-S
    fluid system. American Mineralogist.
    """

    h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, ocs_y, n2_y = (
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    )

    for species in species_list:
        if species == "H2":
            # Shaw and Wones 1964
            logyH2 = (
                P * np.exp((-3.8402 * (T ** (1 / 8))) + 0.541)
                - P**2 * np.exp(-0.1263 * (T**0.5) - 15.98)
                + 300 * np.exp(-0.11901 * T - 5.941) * np.exp(-(P / 300) - 1)
            )

            h2_y = np.exp(logyH2)

        elif species == "CO2":
            # Holland and Powell 1991

            P0 = 5.00  # kbar
            R = cnst.R / 1000  # kJ/mol.K
            a = 741.2 + (-0.10891 * T) + (-3.4203e-4 * T**2)
            b = 3.057  # kJ/kbar/mol
            c = -2.26924e-1 + 7.73793e-5 * T
            d = 1.33790e-2 + (-1.01740e-5 * T)
            P_kb = P / 1000.0  # kbar

            lnf = (
                R * T * np.log(1000.0 * P_kb)
                + b * P_kb
                + (a / (b * (T) ** 0.5))
                * (np.log(R * T + b * P_kb) - np.log(R * T + 2.0 * b * P_kb))
                + (2 / 3) * c * P_kb**1.5
                + (d / 2.0) * P_kb**2.0
            ) / (R * T)

            co2_y = np.exp(lnf) / (P_kb * 1000.0)

        elif species == "N2":
            # Holland and Powell 1991 extension to other species

            Tc = cnst.PTcrit[species][0]
            Pc = cnst.PTcrit[species][1]

            R = cnst.R / 1000  # kJ/mol.K
            a = 5.45963e-5 * (Tc ** (5 / 2) / Pc) - 8.6392e-6 * (Tc ** (3 / 2) / Pc) * T
            b = 9.18301e-4 * (Tc / Pc)  # kJ/kbar/mol
            c = -3.30558e-5 * (Tc / Pc ** (3 / 2)) + T * (2.30524e-6 / Pc ** (3 / 2))
            d = 6.93054e-7 * (Tc / Pc**2) + (-8.38293e-8 / Pc**2) * T
            P_kb = P / 1000.0  # kbar

            lnf = (
                R * T * np.log(1000.0 * P_kb)
                + b * P_kb
                + (a / (b * np.sqrt(T)))
                * (np.log(R * T + b * P_kb) - np.log(R * T + 2.0 * b * P_kb))
                + (2 / 3) * c * P_kb**1.5
                + (d / 2.0) * P_kb**2.0
            ) / (R * T)

            n2_y = np.exp(lnf) / (P_kb * 1000.0)

        elif species == "H2O":
            # Holland and Powell 1991 - above 673 K ONLY

            P0 = 2.0
            P_kb = P / 1000.0  # kbar
            R = cnst.R / 1000  # KJ/ mol.K
            a = (
                1113.4
                + (-0.22291 * (T - 673))
                + (-3.8022e-4 * (T - 673) ** 2)
                + 1.7791e-7 * (T - 673) ** 3
            )
            b = 1.465  # kJ/kbar/mol
            c = -3.025650e-2 + (-5.343144e-6 * T)
            d = -3.2297554e-3 + (2.2215221e-6 * T)

            def find_V(T, P_kb, guess):
                def volMRK(
                    T, P_kb, guess
                ):  # Temp, Pressure, H2O/CO2, a or agas, initial guess for NR
                    # PL: having trimmed down, there may just be one real root for this,
                    # negating need to newton algoritham, just use a solve.
                    def f(x):
                        return (
                            P_kb * (x**3)
                            - R * T * x**2
                            - (b * R * T + (b**2) * P_kb - a / T**0.5) * x
                            - (a * b / T**0.5)
                        )

                    # where Vmrk is designated by x

                    def df(x):  # Defines the differentiated eq of Vmrk
                        return (
                            3 * P_kb * x**2
                            - 2 * R * T * x
                            - (b * R * T + b**2 * P_kb - a / T**0.5)
                        )

                    # Carries out NR on eq f, with initial guess 'guess'
                    # and to an allowable error of tol
                    return newton(f, guess, fprime=df, tol=1e20)

                assert volMRK(T, P_kb, guess) > 0, "MRK Volume is negative"

                if P_kb > P0:
                    V = (
                        volMRK(T, P_kb, guess)
                        + c * ((P_kb - P0) ** 0.5)
                        + d * (P_kb - P0)
                    )
                else:
                    V = volMRK(T, P_kb, guess)

                assert V > 0, "Total volume is negative"
                return V

            guess = (R * T / P_kb) + b
            Z = (P_kb * find_V(T, P_kb, guess)) / (R * T)

            if P_kb > P0:
                lnYvirial = (1 / (R * T)) * (
                    ((2 / 3) * c * (P_kb - P0) ** 1.5) + ((d / 2) * (P_kb - P0) ** 2)
                )  # EQ. A.3
            else:
                lnYvirial = 0

            MRKa = a / (b * R * T**1.5)
            MRKb = (b * P_kb) / (R * T)

            h2o_y = np.exp(
                Z - 1 - np.log(Z - MRKb) - MRKa * np.log(1 + (MRKb / Z)) + lnYvirial
            )  # EQ. A.2

        elif species == "SO2":
            # Shi and Saxena 1992

            Pr = P / cnst.PTcrit[species][1]
            P0r = 1.0 / cnst.PTcrit[species][1]
            Tr = T / cnst.PTcrit[species][0]
            A = (
                0.92854
                + 0.43269e-1 * Tr
                + -0.24671 * Tr**-1
                + 0.24999 * Tr**-2
                + -0.53182 * Tr**-3
                + -0.16461e-1 * np.log(Tr)
            )
            B = (
                0.84866e-3
                + -0.18379e-2 * Tr
                + 0.66787e-1 * Tr**-1
                + -0.29427e-1 * Tr**-2
                + 0.29003e-1 * Tr**-3
                + 0.54808e-2 * np.log(Tr)
            )
            C = (
                -0.35456e-3
                + 0.23316e-4 * Tr
                + 0.94159e-3 * Tr**-1
                + -0.81653e-3 * Tr**-2
                + 0.23154e-3 * Tr**-3
                + 0.55542e-4 * np.log(Tr)
            )
            D = 0.0
            integral = (
                A * np.log(Pr / P0r)
                + B * (Pr - P0r)
                + (C / 2.0) * (Pr**2 - P0r**2)
                + (D / 3.0) * (Pr**3 - P0r**3)
            )

            so2_y = (np.exp(integral) + 0) / P

        else:  # For O2, CO, CH4, H2S, S2 and OCS
            # Shi and Saxena 1992

            def q_less(a, b, c, d, e, Tr):
                # when P < 1000, except H2S
                return a + b * Tr**-1 + c * Tr**-1.5 + d * Tr**-3 + e * Tr**-4

            def q_geq(a, b, c, d, e, f, g, h, Tr):
                # when P>= 1000, except H2S
                return (
                    a
                    + b * Tr
                    + c * Tr**-1
                    + d * Tr**2
                    + e * Tr**-2
                    + f * Tr**3
                    + g * Tr**-3
                    + h * np.log(Tr)
                )

            def Z(a, b, c, d, pr, p0r):
                """Returns the SS(1991) integral Z(p, t)/P, eq. 11"""
                return (
                    a * np.log(pr / p0r)
                    + b * (pr - p0r)
                    + (c / 2.0) * (pr**2 - p0r**2)
                    + (d / 3.0) * (pr**3 - p0r**3)
                )

            def Z1000(pcrit, Tr):
                # z at 1000 bar, except H2S
                Pr = 1000 / pcrit
                P0r = 1 / pcrit

                A = 1
                B = (0.09827 / Tr) + (-0.2709 / Tr**3)
                C = (-0.00103 / Tr**1.5) + (0.01427 / Tr**4)
                D = 0.0

                return Z(A, B, C, D, Pr, P0r)

            def Z5000(pcrit, Tr):
                # z at 5000 bar, except H2S
                Pr = 5000 / pcrit
                P0r = 1000 / pcrit

                A = 1 + (-0.5917 / Tr**2)
                B = 0.09122 / Tr
                C = (-0.0001416 / Tr**2) + (-0.000002835 * np.log(Tr))
                D = 0.0

                return Z(A, B, C, D, Pr, P0r)

            def H2S500(pcrit, Tr):
                # z for H2S at 500 bar
                Pr = 500.0 / pcrit
                P0r = 1.0 / pcrit

                A = q_geq(1.4721, 1.1177, 3.9657, 0, -10.028, 0, 4.5484, -3.8200, Tr)
                B = q_geq(
                    0.16066, 0.10887, 0.29014, 0, -0.99593, 0, -0.18627, -0.45515, Tr
                )
                C = q_geq(
                    -0.28933,
                    -7.0522e-02,
                    0.39828,
                    0,
                    -5.0533e-02,
                    0,
                    0.1176,
                    0.33972,
                    Tr,
                )
                D = 0.0

                return Z(A, B, C, D, Pr, P0r)

            def find_Z0(Pr, Tr, species):
                # Returns the A, B, C and D parameters for Z, and Z0

                if species != "H2S":
                    if P < 1000.0:
                        P0r = 1.0 / cnst.PTcrit[species][1]

                        A = q_less(1, 0, 0, 0, 0, Tr)
                        B = q_less(0, 0.09827, 0, -0.2709, 0, Tr)
                        C = q_less(0, 0, -0.00103, 0, 0.01427, Tr)
                        D = 0
                        Z0 = np.log(1.0)

                    elif P == 1000.0:
                        P0r = 1.0 / cnst.PTcrit[species][1]

                        A = 0
                        B = 0
                        C = 0
                        D = 0
                        Z0 = Z1000(cnst.PTcrit[species][1], Tr)

                    elif P > 1000.0 and P < 5000.0:
                        P0r = 1000.0 / cnst.PTcrit[species][1]

                        A = q_geq(1, 0, 0, 0, -0.5917, 0, 0, 0, Tr)
                        B = q_geq(0, 0, 0.09122, 0, 0, 0, 0, 0, Tr)
                        C = q_geq(0, 0, 0, 0, -0.0001416, 0, 0, -2.835e-06, Tr)
                        D = 0
                        Z0 = Z1000(cnst.PTcrit[species][1], Tr)

                    elif P == 5000.0:
                        P0r = 1000.0 / cnst.PTcrit[species][1]

                        A = 0
                        B = 0
                        C = 0
                        D = 0
                        Z0 = Z5000(cnst.PTcrit[species][1], Tr) + Z1000(
                            cnst.PTcrit[species][1], Tr
                        )

                    else:  # P > 5000 bar
                        P0r = 5000.0 / cnst.PTcrit[species][1]

                        A = q_geq(2.0614, 0, 0, 0, -2.235, 0, 0, -0.3941, Tr)
                        B = q_geq(0, 0, 0.05513, 0, 0.03934, 0, 0, 0, Tr)
                        C = q_geq(0, 0, -1.894e-06, 0, -1.109e-05, 0, -2.189e-05, 0, Tr)
                        D = q_geq(0, 0, 5.053e-11, 0, 0, -6.303e-21, 0, 0, Tr)
                        Z0 = Z5000(cnst.PTcrit[species][1], Tr) + Z1000(
                            cnst.PTcrit[species][1], Tr
                        )

                elif species == "H2S":
                    if P < 500.0:
                        P0r = 1.0 / cnst.PTcrit[species][1]

                        A = q_geq(
                            1.4721, 1.1177, 3.9657, 0, -10.028, 0, 4.5484, -3.8200, Tr
                        )
                        B = q_geq(
                            0.16066,
                            0.10887,
                            0.29014,
                            0,
                            -0.99593,
                            0,
                            -0.18627,
                            -0.45515,
                            Tr,
                        )
                        C = q_geq(
                            -0.28933,
                            -7.0522e-02,
                            0.39828,
                            0,
                            -5.0533e-02,
                            0,
                            0.1176,
                            0.33972,
                            Tr,
                        )
                        D = 0.0
                        Z0 = np.log(1.0)

                    elif P == 500.0:
                        P0r = 1.0 / cnst.PTcrit[species][1]

                        A = 0.0
                        B = 0.0
                        C = 0.0
                        D = 0.0
                        Z0 = H2S500(cnst.PTcrit["H2S"][1], Tr)

                    elif P > 500.0:
                        P0r = 500.0 / cnst.PTcrit["H2S"][1]

                        A = q_geq(
                            0.59941,
                            -1.557e-03,
                            0.04525,
                            0,
                            0.36687,
                            0,
                            -0.79248,
                            0.26058,
                            Tr,
                        )
                        B = q_geq(
                            2.2545e-02,
                            1.7473e-03,
                            0.048253,
                            0,
                            -0.01989,
                            0,
                            0.032794,
                            -0.010985,
                            Tr,
                        )
                        C = q_geq(
                            5.7375e-04,
                            -2.0944e-06,
                            -0.0011894,
                            0,
                            0.0014661,
                            0,
                            -0.00075605,
                            -0.00027985,
                            Tr,
                        )
                        D = 0.0
                        Z0 = H2S500(cnst.PTcrit["H2S"][1], Tr)

                return A, B, C, D, P0r, Z0

            def y(P, T, species):
                """
                Main equation for calculating fugacity coefficients for a
                single `species` using Shi & Saxena 1992.

                Parameters
                ----------
                P : float
                    Pressure (bar)
                T : float
                    Temperature (K)
                species : string
                    The species name

                Returns
                -------
                float
                    fugacity coefficient for `species`.
                """

                # only calibrated to 1 bar, they start increasing gamma again < 1 bar
                # which should happen - gases become more ideal at low pressure.
                if P < 1.0:
                    return 1

                Tr = T / cnst.PTcrit[species][0]
                Pr = P / cnst.PTcrit[species][1]

                A, B, C, D, P0r, Z0 = find_Z0(Pr, Tr, species)

                intZ = Z(A, B, C, D, Pr, P0r)

                return np.exp(intZ + Z0) / P

            if species == "O2":
                o2_y = y(P, T, "O2")
            elif species == "CO":
                co_y = y(P, T, "CO")
            elif species == "CH4":
                ch4_y = y(P, T, "CH4")
            elif species == "S2":
                s2_y = y(P, T, "S2")
            elif species == "H2S":
                h2s_y = y(P, T, "H2S")
            elif species == "OCS":
                ocs_y = y(P, T, "OCS")
            # elif species == 'N2':
            #     n2_y = y(P, T, 'N2')

    return h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, n2_y, ocs_y


def set_Y(sys, molecules):
    """
    Sets the fugacity coefficient for each Molecule instance in the system

    Parameters
    ----------
    sys : ThermoSystem class
        Active instance of the ThermoSystem class
    molecules : list of Molecule classes
        A list of the active Molecule class instances. Order-sensitive.
    """
    if sys.run.GAS_SYS == "OH":
        H2O, O2, H2 = molecules
        H2O.Y, O2.Y, H2.Y = find_Y(sys.P, sys.T, sys.SC)[:3]

    elif sys.run.GAS_SYS == "COH":
        H2O, O2, H2, CO, CO2, CH4 = molecules
        H2O.Y, O2.Y, H2.Y, CO.Y, CO2.Y, CH4.Y = find_Y(sys.P, sys.T, sys.SC)[:6]

    elif sys.run.GAS_SYS == "SOH":
        H2O, O2, H2, S2, SO2, H2S = molecules
        H2O.Y, O2.Y, H2.Y, S2.Y, SO2.Y, H2S.Y = (
            find_Y(sys.P, sys.T, sys.SC)[:3] + find_Y(sys.P, sys.T, sys.SC)[6:9]
        )

    elif sys.run.GAS_SYS == "COHS" and sys.OCS is False:
        H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S = molecules
        H2O.Y, O2.Y, H2.Y, CO.Y, CO2.Y, CH4.Y, S2.Y, SO2.Y, H2S.Y = find_Y(
            sys.P, sys.T, sys.SC
        )[:9]

    elif sys.run.GAS_SYS == "COHS" and sys.OCS is True:
        H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, OCS = molecules
        H2O.Y, O2.Y, H2.Y, CO.Y, CO2.Y, CH4.Y, S2.Y, SO2.Y, H2S.Y, OCS.Y = (
            find_Y(sys.P, sys.T, sys.SC)[:9] + find_Y(sys.P, sys.T, sys.SC)[-1:]
        )

    elif sys.run.GAS_SYS == "COHSN":
        H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2 = molecules
        H2O.Y, O2.Y, H2.Y, CO.Y, CO2.Y, CH4.Y, S2.Y, SO2.Y, H2S.Y, N2.Y = find_Y(
            sys.P, sys.T, sys.SC
        )[:10]


# ------------------------------------------------------------------------
# EQUILIBRIUM CONSTANTS
# ------------------------------------------------------------------------


def get_K(sys, molecules):
    """
    Calculates equilibrium constants using Gibbs Free Energy minimisation.

    Pressure independent, calculated on startup for a fixed temperature,
    based on the methods used in Unterborn (2016).

    Parameters
    ----------
    sys : ThermoSystem class
        Active instance of the ThermoSystem class
    molecules : list of Molecule classes
        A list of the active Molecule class instances. Order-sensitive.

    Returns
    -------
    K : dict
        A dictionary of equilibrium constants.

    References
    ----------
    Unterborn, C.T., Dismukes, E.E. & Panero, W.R. (2016) Scaling the Earth:
    A sensitivity analysis of terrestrial exoplanetary interior models. ApJ.
    """

    K = {}  # Temporary store for K values
    del_Gf_dict = {}  # Stores dG of formation values with molecules as key value pairs
    R = cnst.R

    for mol in molecules:
        # Stores the Gibbs energy in a dictionary with names not object identifiers.
        del_Gf_dict[mol.Mol] = mol.get_G(sys.T)

    # K1: H2 + 1/2 O2 -> H2O
    K1_Gf = (
        del_Gf_dict["H2O"] - del_Gf_dict["H2"] - (0.5 * (del_Gf_dict["O2"]))
    ) * 1000  # so in J rather than kJ
    K["K1"] = exp(-K1_Gf / (R * sys.T))
    # K['K1'] = 10**(12510/T - 0.979*log10(T) + 0.483)

    if "CO2" in sys.SC:
        # K2: CO + 1/2 O2 -> CO2
        K2_Gf = (
            (del_Gf_dict["CO2"] - 0.5 * del_Gf_dict["O2"]) - del_Gf_dict["CO"]
        ) * 1000
        K["K2"] = np.exp(-K2_Gf / (R * sys.T))

        # K3: CH4 + 2 O2 -> CO2 + 2 H2O
        K3_Gf = (
            del_Gf_dict["CO2"]
            + 2 * del_Gf_dict["H2O"]
            - del_Gf_dict["CH4"]
            - 2 * del_Gf_dict["O2"]
        ) * 1000
        K["K3"] = np.exp(-K3_Gf / (R * sys.T))

    if "S2" in sys.SC:
        # K4: 0.5 S2 + H2O -> H2S + 0.5 O2
        K4_Gf = (
            del_Gf_dict["H2S"]
            + 0.5 * del_Gf_dict["O2"]
            - 0.5 * del_Gf_dict["S2"]
            - del_Gf_dict["H2O"]
        ) * 1000
        K["K4"] = exp(-K4_Gf / (R * sys.T))

        # K['K4'] = 10**(-8117/T + 0.188*np.log10(T) - 0.352)

        # K5: 0.5 S2 + O2 -> SO2
        K5_Gf = (
            del_Gf_dict["SO2"] - 0.5 * del_Gf_dict["S2"] - del_Gf_dict["O2"]
        ) * 1000
        K["K5"] = np.exp(-K5_Gf / (R * sys.T))

    if sys.OCS is True:
        # K6: OCS + H2O -> CO2 + H2S
        K6_Gf = (
            del_Gf_dict["CO2"]
            + del_Gf_dict["H2S"]
            - del_Gf_dict["OCS"]
            - del_Gf_dict["H2O"]
        ) * 1000
        K["K6"] = np.exp(-K6_Gf / (R * sys.T))

        # k6 = np.exp(0.482 + 16.166e2*(1/T) +
        #  0.081e-3*T - 5.715e3*(1/T**2) - 2.224e-1*np.log(T))

    return K
