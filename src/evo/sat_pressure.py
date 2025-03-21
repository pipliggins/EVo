"""
Calculates the pressure of volatile saturation for use either alone,
or for initializing an EVo run.
Using system fO2, and melt H2O, CO2 and S contents, gas partial pressures are
calculated for the initial guess pressure,
then the pressure is iterated until the sum of the partial pressures equals the
overburden pressure.
Once the saturation pressure is calculated, WgT is set to 1e-8 and the atomic masses
are calculated. P is then rounded down to the nearest 1 bar,
and the solver is allowed to find the starting conditions at that pressure.

Using this method, the first and second values in mol fraction lists etc WILL NOT MATCH
(recalculate) as they are for different pressures.
"""

import re
import warnings

import numpy as np
from scipy.optimize import fsolve

import evo.constants as cnst
import evo.conversions as cnvs
import evo.messages as msgs
import evo.solubility_laws as sl
import evo.solvgas as sg


def get_molfrac(P, fugacities, gamma):
    """
    Converts fugacities into mole fractions.

    Parameters
    ----------
    P : float
        Pressure (bar)
    fugacities : tuple of floats
        tuple of all 10 fugacities
    gamma : tuple of floats
        tuple of all 10 fugacity coefficients

    Returns
    -------
    tuple of floats
        tuple of all the relevant mole fractions.
    """

    fh2o, fo2, fh2, fco, fco2, fch4, fs2, fso2, fh2s, fn2 = fugacities
    h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, n2_y = gamma

    if fco2 == 0 and fs2 == 0:
        return fh2o / (h2o_y * P), fo2 / (o2_y * P), fh2 / (h2_y * P)
    elif fco2 == 0:
        return (
            fh2o / (h2o_y * P),
            fo2 / (o2_y * P),
            fh2 / (h2_y * P),
            fs2 / (s2_y * P),
            fso2 / (so2_y * P),
            fh2s / (h2s_y * P),
        )
    elif fs2 == 0:
        return (
            fh2o / (h2o_y * P),
            fo2 / (o2_y * P),
            fh2 / (h2_y * P),
            fco / (co_y * P),
            fco2 / (co2_y * P),
            fch4 / (ch4_y * P),
        )
    elif fn2 == 0:
        return (
            fh2o / (h2o_y * P),
            fo2 / (o2_y * P),
            fh2 / (h2_y * P),
            fco / (co_y * P),
            fco2 / (co2_y * P),
            fch4 / (ch4_y * P),
            fs2 / (s2_y * P),
            fso2 / (so2_y * P),
            fh2s / (h2s_y * P),
        )
    else:
        return (
            fh2o / (h2o_y * P),
            fo2 / (o2_y * P),
            fh2 / (h2_y * P),
            fco / (co_y * P),
            fco2 / (co2_y * P),
            fch4 / (ch4_y * P),
            fs2 / (s2_y * P),
            fso2 / (so2_y * P),
            fh2s / (h2s_y * P),
            fn2 / (n2_y * P),
        )


def p_tot(P, fugacities, gamma):
    """
    Returns the difference between the sum of the gas phase partial pressures,
    and the total system pressure.

    Parameters
    ----------
    P : float
        Total system pressure (bar)
    fugacities : tuple of floats
        tuple of all 10 fugacities
    gamma : tuple of floats
        tuple of all 10 fugacity coefficients

    Returns
    -------
    float
        The difference between the sum of the gas phase partial pressures,
        and the total system pressure.
    """

    fh2o, fo2, fh2, fco, fco2, fch4, fs2, fso2, fh2s, fn2 = fugacities
    h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, n2_y = gamma

    if fco2 == 0 and fs2 == 0:
        return fh2o / h2o_y + fo2 / o2_y + fh2 / h2_y - P
    elif fco2 == 0:
        return (
            fh2o / h2o_y
            + fo2 / o2_y
            + fh2 / h2_y
            + fs2 / s2_y
            + fso2 / so2_y
            + fh2s / h2s_y
            - P
        )
    elif fs2 == 0:
        return (
            fh2o / h2o_y
            + fo2 / o2_y
            + fh2 / h2_y
            + fco / co_y
            + fco2 / co2_y
            + fch4 / ch4_y
            - P
        )
    elif fn2 == 0:
        return (
            fh2o / h2o_y
            + fo2 / o2_y
            + fh2 / h2_y
            + fco / co_y
            + fco2 / co2_y
            + fch4 / ch4_y
            + fs2 / s2_y
            + fso2 / so2_y
            + fh2s / h2s_y
            - P
        )
    else:
        return (
            fh2o / h2o_y
            + fo2 / o2_y
            + fh2 / h2_y
            + fco / co_y
            + fco2 / co2_y
            + fch4 / ch4_y
            + fs2 / s2_y
            + fso2 / so2_y
            + fh2s / h2s_y
            + fn2 / n2_y
            - P
        )


def sat_pressure(run, sys, gas, melt, mols):
    """
    Calculates the volatile saturation pressure.

    The volatile saturation pressure is calculated for a given melt
    volatile content. Atomic masses for the system are then calculated
    and fixed, and the pressure is rounded down for the first pressure
    step using `decompress()`.

    Parameters
    ----------
    run : RunDef class
        The active instance of the RunDef class
    sys : ThermoSystem class
        The active instance of the ThermoSystem class
    gas : Gas class
        The active instance of the Gas class
    melt : Melt class
        The active instance of the Melt class
    mols : {[H2O, O2, H2],
            [H2O, O2, H2, CO, CO2, CH4],
            [H2O, O2, H2, S2, SO2, H2S],
            [H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S],
            [H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2]}
        A list of the active Molecule class instances.

    Returns
    -------
    P_sat : float
        The volatile saturation pressure (bar)
    values : list of floats
        The composition of the gas phase at `P_sat`, as mole fractions
    gamma : list of floats
        a list of all calculated fugacity coefficients at `P_sat`, as floats.
    mols : {[H2O, O2, H2],
            [H2O, O2, H2, CO, CO2, CH4],
            [H2O, O2, H2, S2, SO2, H2S],
            [H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S],
            [H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2]}
        An updated list of the active Molecule class instances.
    melt.graphite_sat : bool
        If the melt is graphite saturated at `P_sat`, returns True.
    """

    def get_f(P, sys, melt, gamma):
        """
        Returns gas fugacities at the current pressure according to the melt

        Parameters
        ----------
        P : float
            Current pressure
        sys : ThermoSystem class
            Active instance of the ThermoSystem class
        melt : Melt class
            Active instance of the Melt class
        gamma : list of floats
            list of fugacity coefficients

        Returns
        -------
        fh2o, fo2, fh2, fco, fco2, fch4, fs2, fso2, fh2s, fn2 : tuple of floats
            tuple of all the fugacity coefficients
        """

        O2.Y = gamma[1]

        if run.FO2_buffer_SET:
            fo2 = np.exp(
                cnvs.generate_fo2(sys, sys.run.FO2_buffer_START, sys.run.FO2_buffer, P)
            )

        else:
            fo2 = cnvs.c2fo2(melt.Cm(), sys.T, P * 1e5, sys.run.FO2_MODEL)

        fh2o = sl.h2o_fugacity(run.WTH2O_START, H2O, name=run.H2O_MODEL)
        fh2 = fh2o / (sys.K["K1"] * fo2**0.5)

        if run.WTCO2_SET:
            fco2 = sl.co2_fugacity(
                run.WTCO2_START, CO2, fo2, sys.T, P, melt, name=run.C_MODEL
            )  # This returns graph_fCO2 if graphite saturated.
            fco = fco2 / (sys.K["K2"] * fo2**0.5)
            fch4 = (fco2 * fh2o**2) / (sys.K["K3"] * fo2**2)

        elif run.GRAPHITE_SATURATED:
            melt.graphite_sat = True
            fco2 = sl.graphite_fco2(sys.T, P, fo2)
            fco = fco2 / (sys.K["K2"] * fo2**0.5)
            fch4 = (fco2 * fh2o**2) / (sys.K["K3"] * fo2**2)

        else:
            fco2 = 0
            fco = 0
            fch4 = 0

        if run.SULFUR_SET:
            fs2 = sl.S2_fugacity(
                run.SULFUR_START,
                fo2,
                P,
                sys.T,
                melt,
                run,
                sulfidename=run.SULFIDE_CAPACITY,
                sulfatename=run.SULFATE_CAPACITY,
            )
            fso2 = sys.K["K5"] * fo2 * fs2**0.5
            fh2s = (sys.K["K4"] * fh2o * fs2**0.5) / (fo2**0.5)

        else:
            fs2 = 0
            fso2 = 0
            fh2s = 0

        if run.NITROGEN_SET:
            n2_y = gamma[-1]
            fn2 = sl.n2_fugacity(run.NITROGEN_START, n2_y, fo2, P, name=run.N_MODEL)
        else:
            fn2 = 0

        return fh2o, fo2, fh2, fco, fco2, fch4, fs2, fso2, fh2s, fn2

    # start of function ---------------------------------------------------------------

    if run.GAS_SYS == "OH":
        H2O, O2, H2 = mols
    elif run.GAS_SYS == "COH":
        H2O, O2, H2, CO, CO2, CH4 = mols
    elif run.GAS_SYS == "SOH":
        H2O, O2, H2, S2, SO2, H2S = mols
    elif run.GAS_SYS == "COHS":
        H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S = mols
    elif run.GAS_SYS == "COHSN":
        H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2 = mols

    def find_p(P, sys, melt):
        """
        Function to solve iteratively to find P_sat.

        Parameters
        ----------
        P : float
            Pressure (bar)
        sys : ThermoSystem class
            Active instance of the ThermoSystem class
        melt : Melt class
            Active instance of the Melt class

        Returns
        -------
        dP : float
            The difference between the sum of the gas phase partial
            pressures, and the total system pressure.
        """
        gamma = sg.find_Y(P, sys.T, sys.SC)[:10]
        fugacity = get_f(P, sys, melt, gamma)

        dP = p_tot(P, fugacity, gamma)

        return dP

    with warnings.catch_warnings():
        warnings.filterwarnings("error")
        try:
            P_sat = fsolve(find_p, 1.0, args=(sys, melt))
        except RuntimeWarning as err:
            raise RuntimeError(
                "Failed to find saturation pressure; solver not converging."
            ) from err

    # store the saturation pressure
    P_sat = float(re.sub(r"[\[\]]", "", np.array_str(P_sat)))
    sys.P = P_sat
    gamma = sg.find_Y(P_sat, sys.T, sys.SC)[:10]
    fugacities = get_f(P_sat, sys, melt, gamma)

    # Reset the system fO2 to match that at the saturation pressure
    if run.FO2_buffer_SET:
        sys.FO2 = cnvs.generate_fo2(
            sys, sys.run.FO2_buffer_START, sys.run.FO2_buffer, sys.P
        )

    else:
        sys.FO2 = np.log(cnvs.c2fo2(melt.Cm(), sys.T, sys.Ppa, sys.run.FO2_MODEL))

    # Updates the fugacity coefficients stored in the molecule classes
    sg.set_Y(sys, mols)

    # set the gas fraction at saturation to 1e-6 wt%
    sys.WgT[0] = 1e-8

    # Find the atomic masses at saturation, enforcing sum(mol fracs) = 1

    if run.GAS_SYS == "OH":
        mO2 = get_molfrac(sys.P, fugacities, gamma)[1]

        mH2 = ((1.0 - mO2) * H2O.Y) / (
            sys.K["K1"] * H2.Y * ((O2.Y * mO2 * sys.P) ** 0.5) + H2O.Y
        )

        mH2O = 1 - mO2 - mH2

        lst = {"o2": mO2, "h2": mH2, "h2o": mH2O}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

        # WtO
        sys.atomicM["o"] = cnst.m["o"] * (
            (sys.WgT[0] / sum(mjMj)) * (mH2O + 2 * mO2)
            + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
        )

        # WtH(total)
        sys.atomicM["h"] = (
            2
            * cnst.m["h"]
            * (
                (sys.WgT[0] / sum(mjMj)) * (mH2O + mH2)
                + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
                + sl.h2_melt(mH2, H2, sys.P, melt, name=run.H2_MODEL)
            )
        )

        sys.atomicM["c"], sys.atomicM["s"], sys.atomicM["n"] = 0.0, 0.0, 0.0

        # add gas speciation to lists to act as initial guess for first 'proper' P step
        lists = [gas.mH2O, gas.mO2, gas.mH2]
        empty_list = [gas.mCO, gas.mCO2, gas.mCH4, gas.mS2, gas.mSO2, gas.mH2S, gas.mN2]
        values = [mH2O, mO2, mH2]

    elif run.GAS_SYS == "COH":
        mH2O, mO2, mH2 = get_molfrac(sys.P, fugacities, gamma)[:3]

        mCH4 = (1 - mH2 - mO2 - mH2O) / (
            1
            + (
                (sys.K["K3"] * CH4.Y * (O2.Y * mO2) ** 2)
                / (CO2.Y * (H2O.Y * mH2O) ** 2)
            )
            + (
                (sys.K["K3"] * CH4.Y * (O2.Y * mO2) ** 2)
                / (
                    sys.K["K2"]
                    * CO.Y
                    * (O2.Y * mO2 * sys.P) ** 0.5
                    * (H2O.Y * mH2O) ** 2
                )
            )
        )

        mCO2 = (
            sys.K["K2"]
            * (O2.Y * mO2 * sys.P) ** 0.5
            * CO.Y
            * (1 - mH2 - mO2 - mH2O - mCH4)
        ) / (CO2.Y + sys.K["K2"] * (O2.Y * mO2 * sys.P) ** 0.5 * CO.Y)

        mCO = 1 - mH2O - mH2 - mO2 - mCH4 - mCO2

        if run.GRAPHITE_SATURATED is True:
            melt.graph_current = run.GRAPHITE_START / cnst.m["c"]

        elif run.GRAPHITE_SATURATED is False and melt.graphite_sat is True:
            msgs.graphite_warn_saturation(
                melt, (CO2.Y * mCO2 * sys.P), (O2.Y * mO2 * sys.P), CO2
            )

        lst = {"o2": mO2, "h2": mH2, "h2o": mH2O, "co": mCO, "co2": mCO2, "ch4": mCH4}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

        # WtO
        sys.atomicM["o"] = cnst.m["o"] * (
            (sys.WgT[0] / sum(mjMj)) * (mH2O + 2 * mO2 + 2 * mCO2 + mCO)
            + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
            + 2
            * sl.co2_melt(
                (CO2.Y * mCO2 * sys.P),
                CO2,
                (O2.Y * mO2 * sys.P),
                sys.T,
                sys.P,
                melt,
                name=run.C_MODEL,
            )
            + sl.co_melt((CO.Y * mCO * sys.P), sys.P, name=run.CO_MODEL)
        )

        # WtH(total)
        sys.atomicM["h"] = (
            2
            * cnst.m["h"]
            * (
                (sys.WgT[0] / sum(mjMj)) * (mH2O + mH2 + 2 * mCH4)
                + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
                + sl.h2_melt(mH2, H2, sys.P, melt, name=run.H2_MODEL)
                + 2 * sl.ch4_melt((CH4.Y * mCH4 * sys.P), sys.P, name=run.CH4_MODEL)
            )
        )

        # WtC
        sys.atomicM["c"] = cnst.m["c"] * (
            (sys.WgT[0] / sum(mjMj)) * (mCO + mCO2 + mCH4)
            + sl.co2_melt(
                (CO2.Y * mCO2 * sys.P),
                CO2,
                (O2.Y * mO2 * sys.P),
                sys.T,
                sys.P,
                melt,
                name=run.C_MODEL,
            )
            + sl.ch4_melt((CH4.Y * mCH4 * sys.P), sys.P, name=run.CH4_MODEL)
            + sl.co_melt((CO.Y * mCO * sys.P), sys.P, name=run.CO_MODEL)
            + melt.graph_current
        )

        sys.atomicM["s"] = 0.0

        sys.atomicM["n"] = 0.0

        lists = [gas.mH2O, gas.mO2, gas.mH2, gas.mCO, gas.mCO2, gas.mCH4]
        empty_list = [gas.mS2, gas.mSO2, gas.mH2S, gas.mN2]
        values = [mH2O, mO2, mH2, mCO, mCO2, mCH4]

    elif run.GAS_SYS == "SOH":
        mH2O, mO2, mH2 = get_molfrac(sys.P, fugacities, gamma)[:3]

        # Quadratic solve for SO2

        def quadratic():
            a = SO2.Y**2 / (((sys.K["K5"] * O2.Y * mO2) ** 2) * S2.Y * sys.P)
            b = 1 + (
                (sys.K["K4"] * H2O.Y * mH2O * SO2.Y)
                / (sys.K["K5"] * H2S.Y * ((O2.Y * mO2) ** 1.5) * np.sqrt(sys.P))
            )
            c = -(1 - mH2 - mO2 - mH2O)

            mSO2 = [x for x in np.roots([a, b, c]) if x > 0 and x < 1]
            return mSO2[0]

        mSO2 = quadratic()

        mS2 = (SO2.Y * mSO2) ** 2 / ((sys.K["K5"] * O2.Y * mO2) ** 2 * S2.Y * sys.P)

        mH2S = (1 - mO2 - mS2 - mSO2) / (
            1
            + (H2S.Y / (sys.K["K1"] * sys.K["K4"] * H2.Y * (S2.Y * mS2 * sys.P) ** 0.5))
            + (
                (H2S.Y * (O2.Y * mO2) ** 0.5)
                / (H2O.Y * sys.K["K4"] * (S2.Y * mS2) ** 0.5)
            )
        )

        lst = {"o2": mO2, "h2": mH2, "h2o": mH2O, "so2": mSO2, "s2": mS2, "h2s": mH2S}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

        # WtO
        sys.atomicM["o"] = cnst.m["o"] * (
            (sys.WgT[0] / sum(mjMj)) * (mH2O + 2 * mO2 + 2 * mSO2)
            + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
        )

        # WtH
        sys.atomicM["h"] = (
            2
            * cnst.m["h"]
            * (
                (sys.WgT[0] / sum(mjMj)) * (mH2O + mH2 + mH2S)
                + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
                + sl.h2_melt(mH2, H2, sys.P, melt, name=run.H2_MODEL)
            )
        )

        # WtS
        sys.atomicM["s"] = cnst.m["s"] * (
            (sys.WgT[0] / sum(mjMj)) * (2 * mS2 + mH2S + mSO2)
            + sl.sulfide_melt(
                (S2.Y * mS2 * sys.P),
                (O2.Y * mO2 * sys.P),
                sys.P,
                sys.T,
                melt,
                name=run.SULFIDE_CAPACITY,
            )
            + sl.sulfate_melt(
                (S2.Y * mS2 * sys.P),
                (O2.Y * mO2 * sys.P),
                sys.P,
                sys.T,
                melt,
                run,
                name=run.SULFATE_CAPACITY,
            )
        )

        sys.atomicM["c"] = 0.0

        sys.atomicM["n"] = 0.0

        lists = [gas.mH2O, gas.mO2, gas.mH2, gas.mS2, gas.mSO2, gas.mH2S]
        empty_list = [gas.mCO, gas.mCO2, gas.mCH4, gas.mN2]
        values = [mH2O, mO2, mH2, mS2, mSO2, mH2S]

    elif run.GAS_SYS == "COHS":
        mH2O, mO2, mH2, mS2, mSO2, mH2S = (
            get_molfrac(sys.P, fugacities, gamma)[:3]
            + get_molfrac(sys.P, fugacities, gamma)[6:]
        )

        mCH4 = (1 - mH2 - mO2 - mH2O - mS2 - mSO2 - mH2S) / (
            1
            + (
                (sys.K["K3"] * CH4.Y * (O2.Y * mO2) ** 2)
                / (CO2.Y * (H2O.Y * mH2O) ** 2)
            )
            + (
                (sys.K["K3"] * CH4.Y * (O2.Y * mO2) ** 2)
                / (
                    sys.K["K2"]
                    * CO.Y
                    * (O2.Y * mO2 * sys.P) ** 0.5
                    * (H2O.Y * mH2O) ** 2
                )
            )
        )

        mCO2 = (sys.K["K3"] * CH4.Y * mCH4 * (O2.Y * mO2) ** 2) / (
            (H2O.Y * mH2O) ** 2 * CO2.Y
        )

        mCO = CO2.Y * mCO2 / (sys.K["K2"] * CO.Y * (O2.Y * mO2 * sys.P) ** 0.5)

        if run.GRAPHITE_SATURATED is True:
            melt.graph_current = run.GRAPHITE_START / cnst.m["c"]

        elif run.GRAPHITE_SATURATED is False and melt.graphite_sat is True:
            msgs.graphite_warn_saturation(
                melt, (CO2.Y * mCO2 * sys.P), (O2.Y * mO2 * sys.P), CO2
            )

        lst = {
            "o2": mO2,
            "h2": mH2,
            "h2o": mH2O,
            "co": mCO,
            "co2": mCO2,
            "ch4": mCH4,
            "so2": mSO2,
            "s2": mS2,
            "h2s": mH2S,
        }
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

        # WtO
        sys.atomicM["o"] = cnst.m["o"] * (
            (sys.WgT[0] / sum(mjMj)) * (mH2O + 2 * mO2 + 2 * mCO2 + mCO + 2 * mSO2)
            + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
            + 2
            * sl.co2_melt(
                (CO2.Y * mCO2 * sys.P),
                CO2,
                (O2.Y * mO2 * sys.P),
                sys.T,
                sys.P,
                melt,
                name=run.C_MODEL,
            )
            + sl.co_melt((CO.Y * mCO * sys.P), sys.P, name=run.CO_MODEL)
        )

        # WtH
        sys.atomicM["h"] = (
            2
            * cnst.m["h"]
            * (
                (sys.WgT[0] / sum(mjMj)) * (mH2O + mH2 + 2 * mCH4 + mH2S)
                + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
                + sl.h2_melt(mH2, H2, sys.P, melt, name=run.H2_MODEL)
                + 2 * sl.ch4_melt((CH4.Y * mCH4 * sys.P), sys.P, name=run.CH4_MODEL)
            )
        )

        # WtC
        sys.atomicM["c"] = cnst.m["c"] * (
            (sys.WgT[0] / sum(mjMj)) * (mCO + mCO2 + mCH4)
            + sl.co2_melt(
                (CO2.Y * mCO2 * sys.P),
                CO2,
                (O2.Y * mO2 * sys.P),
                sys.T,
                sys.P,
                melt,
                name=run.C_MODEL,
            )
            + sl.co_melt((CO.Y * mCO * sys.P), sys.P, name=run.CO_MODEL)
            + sl.ch4_melt((CH4.Y * mCH4 * sys.P), sys.P, name=run.CH4_MODEL)
            + melt.graph_current
        )

        # WtS
        sys.atomicM["s"] = cnst.m["s"] * (
            (sys.WgT[0] / sum(mjMj)) * (2 * mS2 + mH2S + mSO2)
            + sl.sulfide_melt(
                (S2.Y * mS2 * sys.P),
                (O2.Y * mO2 * sys.P),
                sys.P,
                sys.T,
                melt,
                name=run.SULFIDE_CAPACITY,
            )
            + sl.sulfate_melt(
                (S2.Y * mS2 * sys.P),
                (O2.Y * mO2 * sys.P),
                sys.P,
                sys.T,
                melt,
                run,
                name=run.SULFATE_CAPACITY,
            )
        )

        sys.atomicM["n"] = 0.0

        lists = [
            gas.mH2O,
            gas.mO2,
            gas.mH2,
            gas.mCO,
            gas.mCO2,
            gas.mCH4,
            gas.mS2,
            gas.mSO2,
            gas.mH2S,
        ]
        empty_list = [gas.mN2]
        values = [mH2O, mO2, mH2, mCO, mCO2, mCH4, mS2, mSO2, mH2S]

    elif run.GAS_SYS == "COHSN":
        mH2O, mO2, mH2, mS2, mSO2, mH2S, mN2 = (
            get_molfrac(sys.P, fugacities, gamma)[:3]
            + get_molfrac(sys.P, fugacities, gamma)[6:]
        )

        mCH4 = (1 - mH2 - mO2 - mH2O - mS2 - mSO2 - mH2S - mN2) / (
            1
            + (
                (sys.K["K3"] * CH4.Y * (O2.Y * mO2) ** 2)
                / (CO2.Y * (H2O.Y * mH2O) ** 2)
            )
            + (
                (sys.K["K3"] * CH4.Y * (O2.Y * mO2) ** 2)
                / (
                    sys.K["K2"]
                    * CO.Y
                    * (O2.Y * mO2 * sys.P) ** 0.5
                    * (H2O.Y * mH2O) ** 2
                )
            )
        )

        mCO2 = (sys.K["K3"] * CH4.Y * mCH4 * (O2.Y * mO2) ** 2) / (
            (H2O.Y * mH2O) ** 2 * CO2.Y
        )

        mCO = CO2.Y * mCO2 / (sys.K["K2"] * CO.Y * (O2.Y * mO2 * sys.P) ** 0.5)

        if run.GRAPHITE_SATURATED is True:
            melt.graph_current = run.GRAPHITE_START / cnst.m["c"]

        elif run.GRAPHITE_SATURATED is False and melt.graphite_sat is True:
            msgs.graphite_warn_saturation(
                melt, (CO2.Y * mCO2 * sys.P), (O2.Y * mO2 * sys.P), CO2
            )

        lst = {
            "o2": mO2,
            "h2": mH2,
            "h2o": mH2O,
            "co": mCO,
            "co2": mCO2,
            "ch4": mCH4,
            "so2": mSO2,
            "s2": mS2,
            "h2s": mH2S,
            "n2": mN2,
        }
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

        # WtO
        sys.atomicM["o"] = cnst.m["o"] * (
            (sys.WgT[0] / sum(mjMj)) * (mH2O + 2 * mO2 + 2 * mCO2 + mCO + 2 * mSO2)
            + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
            + 2
            * sl.co2_melt(
                (CO2.Y * mCO2 * sys.P),
                CO2,
                (O2.Y * mO2 * sys.P),
                sys.T,
                sys.P,
                melt,
                name=run.C_MODEL,
            )
            + sl.co_melt((CO.Y * mCO * sys.P), sys.P, name=run.CO_MODEL)
        )

        # WtH
        sys.atomicM["h"] = (
            2
            * cnst.m["h"]
            * (
                (sys.WgT[0] / sum(mjMj)) * (mH2O + mH2 + 2 * mCH4 + mH2S)
                + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
                + sl.h2_melt(mH2, H2, sys.P, melt, name=run.H2_MODEL)
                + 2 * sl.ch4_melt((CH4.Y * mCH4 * sys.P), sys.P, name=run.CH4_MODEL)
            )
        )

        # WtC
        sys.atomicM["c"] = cnst.m["c"] * (
            (sys.WgT[0] / sum(mjMj)) * (mCO + mCO2 + mCH4)
            + sl.co2_melt(
                (CO2.Y * mCO2 * sys.P),
                CO2,
                (O2.Y * mO2 * sys.P),
                sys.T,
                sys.P,
                melt,
                name=run.C_MODEL,
            )
            + sl.ch4_melt((CH4.Y * mCH4 * sys.P), sys.P, name=run.CH4_MODEL)
            + sl.co_melt((CO.Y * mCO * sys.P), sys.P, name=run.CO_MODEL)
            + melt.graph_current
        )

        # WtS
        sys.atomicM["s"] = cnst.m["s"] * (
            (sys.WgT[0] / sum(mjMj)) * (2 * mS2 + mH2S + mSO2)
            + sl.sulfide_melt(
                (S2.Y * mS2 * sys.P),
                (O2.Y * mO2 * sys.P),
                sys.P,
                sys.T,
                melt,
                name=run.SULFIDE_CAPACITY,
            )
            + sl.sulfate_melt(
                (S2.Y * mS2 * sys.P),
                (O2.Y * mO2 * sys.P),
                sys.P,
                sys.T,
                melt,
                run,
                name=run.SULFATE_CAPACITY,
            )
        )

        # WtN
        sys.atomicM["n"] = cnst.m["n"] * (
            (sys.WgT[0] / sum(mjMj)) * (2 * mN2)
            + sl.n_melt(mN2, (O2.Y * mO2 * sys.P), sys.P, name=run.N_MODEL)
        )

        lists = [
            gas.mH2O,
            gas.mO2,
            gas.mH2,
            gas.mCO,
            gas.mCO2,
            gas.mCH4,
            gas.mS2,
            gas.mSO2,
            gas.mH2S,
            gas.mN2,
        ]
        empty_list = []
        values = [mH2O, mO2, mH2, mCO, mCO2, mCH4, mS2, mSO2, mH2S, mN2]

    # add gas speciation to lists to act as initial guess for first 'proper' P step
    for ls, val in zip(lists, values):
        ls.append(val)

    for ls in empty_list:
        ls.append(np.float64(0))

    # Round P down to nearest bar ready for solver to find WgT at 'starting' pressure.
    sys.P = cnvs.round_down(sys.P, sys.P_step)
    run.P_START = sys.P

    return P_sat, values, gamma, mols, melt.graphite_sat


def satp_writeout(sys, melt, gas, P, values, gamma, mols, graph_sat=False):
    """
    Returns a string in the same style as the output csv file for the
    saturation pressure.

    Parameters
    ----------
    sys : ThermoSystem class
        Active instance of the ThermoSystem class
    melt : Melt class
        Active instance of the Melt class
    gas : Gas class
        Active instance of the Gas class
    P : float
        The volatile saturation pressure
    values : list of floats
        The composition of the gas phase in mole fractions
    gamma : list of floats
        list of the gas phase fugacity coefficients
    mols : list of Molecule classes
        A list of the active Molecule classes
    graph_sat : bool, optional
        If True, indicates that the melt was graphite saturated at the
        point of volatile saturation, by default False

    Returns
    -------
    string
        A single string to be put at the top of the output csv file.
    """

    if sys.run.GAS_SYS == "OH":
        H2O, O2, H2 = mols
        mH2O, mO2, mH2 = tuple(values)
        h2oy, o2y, h2y = gamma[:3]

        wts = cnvs.mols2wts(H2O=mH2O, O2=mO2, H2=mH2)
        empty = ["CO2", "CO", "CH4", "SO2", "S2", "H2S", "N2"]
        for i in empty:
            wts[i] = 0.0
        M = cnvs.mean_mol_wt(H2O=mH2O, O2=mO2, H2=mH2)

    elif sys.run.GAS_SYS == "COH":
        H2O, O2, H2, CO, CO2, CH4 = mols
        mH2O, mO2, mH2, mCO, mCO2, mCH4 = tuple(values)
        h2oy, o2y, h2y, coy, co2y, ch4y = gamma[:6]

        wts = cnvs.mols2wts(H2O=mH2O, O2=mO2, H2=mH2, CO=mCO, CO2=mCO2, CH4=mCH4)
        empty = ["SO2", "S2", "H2S", "N2"]
        for i in empty:
            wts[i] = 0.0
        M = cnvs.mean_mol_wt(H2O=mH2O, O2=mO2, H2=mH2, CO=mCO, CO2=mCO2, CH4=mCH4)

    elif sys.run.GAS_SYS == "SOH":
        H2O, O2, H2, S2, SO2, H2S = mols
        mH2O, mO2, mH2, mS2, mSO2, mH2S = tuple(values)
        h2oy, o2y, h2y, s2y, so2y, h2sy = gamma[:3] + gamma[6:9]

        wts = cnvs.mols2wts(H2O=mH2O, O2=mO2, H2=mH2, S2=mS2, SO2=mSO2, H2S=mH2S)
        empty = ["CO2", "CO", "CH4", "N2"]
        for i in empty:
            wts[i] = 0.0
        M = cnvs.mean_mol_wt(H2O=mH2O, O2=mO2, H2=mH2, S2=mS2, SO2=mSO2, H2S=mH2S)

    elif sys.run.GAS_SYS == "COHS":
        H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S = mols
        mH2O, mO2, mH2, mCO, mCO2, mCH4, mS2, mSO2, mH2S = tuple(values)
        h2oy, o2y, h2y, coy, co2y, ch4y, s2y, so2y, h2sy = gamma[:9]

        wts = cnvs.mols2wts(
            H2O=mH2O,
            O2=mO2,
            H2=mH2,
            CO=mCO,
            CO2=mCO2,
            CH4=mCH4,
            S2=mS2,
            SO2=mSO2,
            H2S=mH2S,
        )
        empty = ["N2"]
        for i in empty:
            wts[i] = 0.0
        M = cnvs.mean_mol_wt(
            H2O=mH2O,
            O2=mO2,
            H2=mH2,
            CO=mCO,
            CO2=mCO2,
            CH4=mCH4,
            S2=mS2,
            SO2=mSO2,
            H2S=mH2S,
        )

    elif sys.run.GAS_SYS == "COHSN":
        H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2 = mols
        mH2O, mO2, mH2, mCO, mCO2, mCH4, mS2, mSO2, mH2S, mN2 = tuple(values)
        h2oy, o2y, h2y, coy, co2y, ch4y, s2y, so2y, h2sy, n2y = gamma

        wts = cnvs.mols2wts(
            H2O=mH2O,
            O2=mO2,
            H2=mH2,
            CO=mCO,
            CO2=mCO2,
            CH4=mCH4,
            S2=mS2,
            SO2=mSO2,
            H2S=mH2S,
            N2=mN2,
        )
        M = cnvs.mean_mol_wt(
            H2O=mH2O,
            O2=mO2,
            H2=mH2,
            CO=mCO,
            CO2=mCO2,
            CH4=mCH4,
            S2=mS2,
            SO2=mSO2,
            H2S=mH2S,
            N2=mN2,
        )

    rho_melt = melt.rho(P=(P * 1e5))

    GvF = 1 / (
        1
        + ((M * P * 1e5 * (1 - sys.WgT[0])) / (cnst.R * sys.T * rho_melt * sys.WgT[0]))
    )

    rho_gas = ((M / 1000) * P * 1e5) / (cnst.R * sys.T)
    rho_bulk = (rho_gas * GvF) + (rho_melt * (1 - GvF))

    if sys.run.ATOMIC_MASS_SET is False:
        melt_h2o = sys.run.WTH2O_START
    else:
        melt_h2o = (
            sl.h2o_melt(mH2O, H2O, P, name=sys.run.H2O_MODEL, Y=h2oy) * cnst.m["h2o"]
        )
    melt_h2 = sl.h2_melt(mH2, H2, P, melt, name=sys.run.H2_MODEL, Y=h2y) * cnst.m["h2"]

    if (
        sys.run.GAS_SYS == "COHS"
        or sys.run.GAS_SYS == "COH"
        or sys.run.GAS_SYS == "COHSN"
    ):
        if sys.run.ATOMIC_MASS_SET is False and graph_sat is False:
            melt_co2 = sys.run.WTCO2_START
        else:
            melt_co2 = (
                sl.co2_melt(
                    (co2y * mCO2 * P),
                    CO2,
                    (o2y * mO2 * P),
                    sys.T,
                    P,
                    melt,
                    name=sys.run.C_MODEL,
                )
                * cnst.m["co2"]
            )

        if sys.run.CO_MODEL != "None":
            melt_co = (
                sl.co_melt((coy * mCO * P), P, name=sys.run.CO_MODEL) * cnst.m["co"]
            )
        else:
            melt_co = 0.0

        if sys.run.CH4_MODEL != "None":
            melt_ch4 = (
                sl.ch4_melt((ch4y * mCH4 * P), P, name=sys.run.CH4_MODEL)
                * cnst.m["ch4"]
            )
        else:
            melt_ch4 = 0.0

        if graph_sat is True:
            graphite = (
                (sys.atomicM["c"] / cnst.m["c"])
                - (
                    ((wts["CO2"] * sys.WgT[0]) + melt_co2) / cnst.m["co2"]
                    + (wts["CO"] * sys.WgT[0] + melt_co) / cnst.m["co"]
                    + (wts["CH4"] * sys.WgT[0] + melt_ch4) / cnst.m["ch4"]
                )
            ) * cnst.m["c"]
        else:
            graphite = 0.0

    else:
        melt_co2 = 0.0
        melt_co = 0.0
        melt_ch4 = 0.0
        graphite = 0.0

        mCO, mCO2, mCH4 = 0.0, 0.0, 0.0
        coy, co2y, ch4y = 0.0, 0.0, 0.0

    if (
        sys.run.GAS_SYS == "COHS"
        or sys.run.GAS_SYS == "SOH"
        or sys.run.GAS_SYS == "COHSN"
    ):
        fo2 = o2y * mO2 * P
        melt.cm_dry, F = melt.iron_fraction(
            np.log(fo2), ppa=P * 1e5
        )  # Set the Fe2/Fe3 ratio and dry melt chemistry for the sulfide capacity

        # use mS2 to work back and find amount of S2- in melt
        s2_melt = (
            sl.sulfide_melt(
                (s2y * mS2 * P),
                (o2y * mO2 * P),
                P,
                sys.T,
                melt,
                name=sys.run.SULFIDE_CAPACITY,
            )
            * cnst.m["s"]
        )  # sulfide_melt returns *number of moles* not mass
        s6_melt = (
            sl.sulfate_melt(
                (s2y * mS2 * P),
                (o2y * mO2 * P),
                P,
                sys.T,
                melt,
                sys.run,
                name=sys.run.SULFATE_CAPACITY,
            )
            * cnst.m["s"]
        )
        melt_s = s2_melt + s6_melt

    else:
        mS2, mSO2, mH2S = 0, 0, 0
        s2y, so2y, h2sy = 0, 0, 0

        fo2 = o2y * mO2 * P
        melt.cm_dry, F = melt.iron_fraction(np.log(fo2), ppa=P * 1e5)

        melt_s = 0
        s2_melt = 0
        s6_melt = 0

    if sys.run.GAS_SYS == "COHSN":
        if sys.run.ATOMIC_MASS_SET is False:
            melt_n = sys.run.NITROGEN_START
        else:
            melt_n = (
                sl.n_melt(mN2, (o2y * mO2 * P), P, name=sys.run.N_MODEL) * cnst.m["n"]
            )

    else:
        melt_n, mN2, n2y = 0.0, 0.0, 0.0

    atomico = cnst.m["o"] * (
        ((sys.WgT[0] * wts["H2O"]) + melt_h2o) / cnst.m["h2o"]
        + (2 * sys.WgT[0] * wts["O2"]) / cnst.m["o2"]
        + (sys.WgT[0] * wts["CO"] + melt_co) / cnst.m["co"]
        + (2 * ((sys.WgT[0] * wts["CO2"]) + melt_co2)) / cnst.m["co2"]
        + (2 * (sys.WgT[0] * wts["SO2"])) / cnst.m["so2"]
    )

    fo2_buffer = sys.run.FO2_buffer

    data = {
        "P": P,
        fo2_buffer: cnvs.generate_fo2_buffer(sys, (o2y * mO2 * P), P),
        "fo2": fo2,
        "F": F,
        "rho_bulk": rho_bulk,
        "rho_melt": rho_melt,
        "Exsol_vol%": GvF * 100,
        "Gas_wt": sys.WgT[0] * 100,
        "mol_mass": M,
        "mH2O": gas.mH2O[0],
        "mH2": gas.mH2[0],
        "mO2": gas.mO2[0],
        "mCO2": gas.mCO2[0],
        "mCO": gas.mCO[0],
        "mCH4": gas.mCH4[0],
        "mSO2": gas.mSO2[0],
        "mH2S": gas.mH2S[0],
        "mS2": gas.mS2[0],
        "mN2": gas.mN2[0],
        "wH2O": wts["H2O"],
        "wH2": wts["H2"],
        "wO2": wts["O2"],
        "wCO2": wts["CO2"],
        "wCO": wts["CO"],
        "wCH4": wts["CH4"],
        "wSO2": wts["SO2"],
        "wH2S": wts["H2S"],
        "wS2": wts["S2"],
        "wN2": wts["N2"],
        "fH2O": (h2oy * mH2O * P),
        "fH2": (h2y * mH2 * P),
        "fCO2": (co2y * mCO2 * P),
        "fCO": (coy * mCO * P),
        "fCH4": (ch4y * mCH4 * P),
        "fSO2": (so2y * mSO2 * P),
        "fH2S": (h2sy * mH2S * P),
        "fS2": (s2y * mS2 * P),
        "fN2": (n2y * mN2 * P),
        "mCO2/CO": gas.mCO2[0] / gas.mCO[0] if "C" in sys.run.GAS_SYS else np.nan,
        "mCO2/H2O": gas.mCO2[0] / gas.mH2O[0] if "C" in sys.run.GAS_SYS else np.nan,
        "mCO2/SO2": gas.mCO2[0] / gas.mSO2[0] if "S" in sys.run.GAS_SYS else np.nan,
        "mH2S/SO2": gas.mH2S[0] / gas.mSO2[0] if "S" in sys.run.GAS_SYS else np.nan,
        "H2O_melt": melt_h2o * 100,
        "H2_melt": melt_h2 * 100,
        "CO2_melt": melt_co2 * 100,
        "CO_melt": melt_co * 100,
        "CH4_melt": melt_ch4 * 100,
        "graph_melt": graphite * 100,
        "S2-_melt": s2_melt * 100,
        "S6+_melt": s6_melt * 100,
        "Stot_melt": melt_s * 100,
        "N_melt": melt_n * 100,
        "tot_H": sys.atomicM["h"] * 1000000,
        "tot_C": sys.atomicM["c"] * 1000000,
        "tot_O_gas": atomico * 1000000,
        "tot_O": cnvs.atomicM_calc(sys, melt, gas, "o_tot", 0) * 1000000,
        "tot_S": sys.atomicM["s"] * 1000000,
        "tot_N": sys.atomicM["n"] * 1000000,
    }

    return data
