"""
Writes the results of the run to a file and
contains options to produce a graph of the results.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import evo.conversions as cnvt
import evo.sat_pressure as sat


def get_data(sys, gas, melt, P):
    """
    Creates the dictionary structure containing output data

    Parameters
    ----------
    sys : ThermoSystem class
        Active instance of the ThermoSystem class
    gas : Gas class
        Active instance of the Gas class
    melt : Melt class
        Active instance of the Melt class
    P : [float]
        Pressure path (bar)
    """

    fo2_buffer = sys.run.FO2_buffer

    data = {
        "P": P,
        fo2_buffer: [
            cnvt.generate_fo2_buffer(sys, np.exp(gas.fo2[i]), P[i])
            for i in range(len(P))
        ],
        "fo2": [np.exp(gas.fo2[i]) for i in range(len(P))],
        "F": melt.F[1:] if not sys.run.SINGLE_STEP else [melt.F[-1]],
        "Fe3FeT": melt.Fe3FeT[1:] if not sys.run.SINGLE_STEP else [melt.Fe3FeT[-1]],
        "rho_bulk": sys.rho,
        "rho_melt": melt.rho_store,
        "Exsol_vol%": [(sys.GvF[i] * 100) for i in range(len(P))],
        "Gas_wt": [(sys.WgT[i + 1] * 100) for i in range(len(P))],
        "mol_mass": gas.M,
        "mH2O": gas.mH2O[1:] if not sys.run.SINGLE_STEP else [gas.mH2O[-1]],
        "mH2": gas.mH2[1:] if not sys.run.SINGLE_STEP else [gas.mH2[-1]],
        "mO2": gas.mO2[1:] if not sys.run.SINGLE_STEP else [gas.mO2[-1]],
        "mCO2": gas.mCO2[1:] if not sys.run.SINGLE_STEP else [gas.mCO2[-1]],
        "mCO": gas.mCO[1:] if not sys.run.SINGLE_STEP else [gas.mCO[-1]],
        "mCH4": gas.mCH4[1:] if not sys.run.SINGLE_STEP else [gas.mCH4[-1]],
        "mSO2": gas.mSO2[1:] if not sys.run.SINGLE_STEP else [gas.mSO2[-1]],
        "mH2S": gas.mH2S[1:] if not sys.run.SINGLE_STEP else [gas.mH2S[-1]],
        "mS2": gas.mS2[1:] if not sys.run.SINGLE_STEP else [gas.mS2[-1]],
        "mN2": gas.mN2[1:] if not sys.run.SINGLE_STEP else [gas.mN2[-1]],
        "wH2O": gas.Wt["H2O"],
        "wH2": gas.Wt["H2"],
        "wO2": gas.Wt["O2"],
        "wCO2": gas.Wt["CO2"],
        "wCO": gas.Wt["CO"],
        "wCH4": gas.Wt["CH4"],
        "wSO2": gas.Wt["SO2"],
        "wH2S": gas.Wt["H2S"],
        "wS2": gas.Wt["S2"],
        "wN2": gas.Wt["N2"],
        "fH2O": gas.f["H2O"],
        "fH2": gas.f["H2"],
        "fCO2": gas.f["CO2"],
        "fCO": gas.f["CO"],
        "fCH4": gas.f["CH4"],
        "fSO2": gas.f["SO2"],
        "fH2S": gas.f["H2S"],
        "fS2": gas.f["S2"],
        "fN2": gas.f["N2"],
        "mCO2/CO": (
            np.nan
            if "C" not in sys.run.GAS_SYS
            else (
                [gas.mCO2[i + 1] / gas.mCO[i + 1] for i in range(len(P))]
                if not sys.run.SINGLE_STEP
                else [gas.mCO2[-1] / gas.mCO[-1]]
            )
        ),
        "mCO2/H2O": (
            np.nan
            if "C" not in sys.run.GAS_SYS
            else (
                [gas.mCO2[i + 1] / gas.mH2O[i + 1] for i in range(len(P))]
                if not sys.run.SINGLE_STEP
                else [gas.mCO2[-1] / gas.mH2O[-1]]
            )
        ),
        "mCO2/SO2": (
            np.nan
            if "S" not in sys.run.GAS_SYS
            else (
                [gas.mCO2[i + 1] / gas.mSO2[i + 1] for i in range(len(P))]
                if not sys.run.SINGLE_STEP
                else [gas.mCO2[-1] / gas.mSO2[-1]]
            )
        ),
        "mH2S/SO2": (
            np.nan
            if "S" not in sys.run.GAS_SYS
            else (
                [gas.mH2S[i + 1] / gas.mSO2[i + 1] for i in range(len(P))]
                if not sys.run.SINGLE_STEP
                else [gas.mH2S[-1] / gas.mSO2[-1]]
            )
        ),
        "H2O_melt": [melt.h2o[i] * 100 for i in range(len(P))],
        "H2_melt": [melt.h2[i] * 100 for i in range(len(P))],
        "CO2_melt": [melt.co2[i] * 100 for i in range(len(P))],
        "CO_melt": [melt.co[i] * 100 for i in range(len(P))],
        "CH4_melt": [melt.ch4[i] * 100 for i in range(len(P))],
        "graph_melt": [melt.graphite[i] * 100 for i in range(len(P))],
        "S2-_melt": [melt.sulfide[i] * 100 for i in range(len(P))],
        "S6+_melt": [melt.sulfate[i] * 100 for i in range(len(P))],
        "Stot_melt": [melt.s[i] * 100 for i in range(len(P))],
        "N_melt": [melt.n[i] * 100 for i in range(len(P))],
        "tot_H": [
            cnvt.atomicM_calc(sys, melt, gas, "h", i) * 1000000 for i in range(len(P))
        ],
        "tot_C": [
            cnvt.atomicM_calc(sys, melt, gas, "c", i) * 1000000 for i in range(len(P))
        ],
        "tot_O_gas": [
            cnvt.atomicM_calc(sys, melt, gas, "o", i) * 1000000 for i in range(len(P))
        ],
        "tot_O": [
            cnvt.atomicM_calc(sys, melt, gas, "o_tot", i) * 1000000
            for i in range(len(P))
        ],
        "tot_S": [
            cnvt.atomicM_calc(sys, melt, gas, "s", i) * 1000000 for i in range(len(P))
        ],
        "tot_N": [
            cnvt.atomicM_calc(sys, melt, gas, "n", i) * 1000000 for i in range(len(P))
        ],
    }

    return data


def writeout_file(sys, gas, melt, P, crashed=False):
    """
    At the end of a run, writes data to a single output CSV file

    Parameters
    ----------
    sys : ThermoSystem class
        Active instance of the ThermoSystem class
    gas : Gas class
        Active instance of the Gas class
    melt : Melt class
        Active instance of the Melt class
    P : list
        Pressure path (bar)
    crashed: bool, False
        If EVo has crashed before reaching its final pressure,
        crashed=true will append 'CRASHED' to the datafile name.
    """

    data = get_data(sys, gas, melt, P)

    df = pd.DataFrame(data)

    filepath = sys.run.results_folder
    if not filepath.exists():
        filepath.mkdir()

    if not crashed:
        file_name = (
            f"dgs_output_{sys.run.COMPOSITION}_{sys.run.GAS_SYS}_{sys.run.RUN_TYPE}"
            f"_{sys.T:.0f}K.csv"
        )
    else:
        file_name = (
            f"dgs_output_CRASHED_{sys.run.COMPOSITION}_{sys.run.GAS_SYS}"
            f"_{sys.run.RUN_TYPE}_{sys.T:.0f}K.csv"
        )

    output_path = filepath / file_name

    if sys.run.FIND_SATURATION is True or sys.run.ATOMIC_MASS_SET is True:
        df_sat = pd.DataFrame(
            sat.satp_writeout(
                sys,
                melt,
                gas,
                sys.sat_conditions[0],
                sys.sat_conditions[1],
                sys.sat_conditions[2],
                sys.sat_conditions[3],
                graph_sat=sys.sat_conditions[4],
            ),
            index=[0],
        )

        df = pd.concat([df_sat, df])

    df.to_csv(output_path, index=False)
    return df


def writeout_figs(sys, melt, gas, out, P):
    """
    Selects the figures to produce at the end of a successful run.

    If an output file is not used to select the required figures, all
    options are plotted.

    Parameters
    ----------
    sys : ThermoSystem class
        Active instance of the ThermoSystem class
    melt : Melt class
        Active instance of the Melt class
    gas : Gas class
        Active instance of the Gas class
    out : Output class or None
        Active instance of the Output class, or None
    P : list of float
        List of the pressures each step was calculated at (bar)
    """

    filepath = sys.run.results_folder
    # filelist = glob.glob(str(filepath / "*.png"))
    filelist = list(filepath.glob("*.png"))
    # Removes previous files so if output specification is changed
    # there is no confusion as to up to date files.

    for file in filelist:
        file.unlink()
    if (
        out is not None
    ):  # If an output file listing requested figures has been included:
        if out.Plt_melt_species:
            plot_meltspecies(melt, P, filepath)

        if out.Plt_gas_species_wt:
            plot_gasspecies_wt(gas, P, filepath)

        if out.Plt_gas_species_mol:
            plot_gasspecies_mol(gas, P, filepath)

        if out.Plt_gas_fraction:
            plot_gasfraction(sys, P, filepath)

        if out.Plt_fo2_dFMQ:
            plot_fo2FMQ(melt, gas, P, filepath)

        if out.Plt_gas_fraction is not None:
            pass
            # plot

    else:  # else plot every option
        plot_meltspecies(melt, P, filepath)
        plot_gasspecies_wt(gas, P, filepath)
        plot_gasspecies_mol(gas, P, filepath)
        plot_gasfraction(sys, P, filepath)
        plot_fo2FMQ(melt, gas, P, filepath)


# fO2 and dFMQ
def plot_fo2FMQ(melt, gas, P, path):
    """Plots the fO2 relative to FMQ against pressure"""
    plt.plot(P, melt.fmq)
    plt.xscale("log")
    plt.xlabel("Pressure (bars)")
    plt.ylabel(r"$\Delta$ FMQ")
    plt.savefig(path / "FMQ.png")
    plt.close()

    fo2 = []
    for x in gas.fo2:
        fo2.append(np.log10(np.exp(x)))

    plt.plot(P, fo2)
    # plt.xscale('log')
    plt.xlabel("Pressure (bars)")
    plt.ylabel("log(10) fO2")
    plt.savefig(path / "fO2.png")
    plt.close()


# Gas speciation mol
def plot_gasspecies_mol(gas, P, path):
    """Plots the gas speciation (mole fraction) vs pressure"""

    plt.plot(gas.mH2O[1:], P, c="darkblue")
    plt.plot(gas.mH2[1:], P, c="steelblue")
    plt.plot(gas.mO2[1:], P, c="firebrick")
    plt.annotate("H2O", xy=[gas.mH2O[-1], P[-1]], color="darkblue")
    plt.annotate("H2", xy=[gas.mH2[-1], P[-1]], color="steelblue")
    plt.annotate("O2", xy=[gas.mO2[-1], P[-1]], color="firebrick")
    if (
        gas.sys.run.GAS_SYS == "COH"
        or gas.sys.run.GAS_SYS == "COHS"
        or gas.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(gas.mCO2[1:], P, c="saddlebrown")
        plt.plot(gas.mCO[1:], P, c="darkgreen")
        plt.plot(gas.mCH4[1:], P, c="forestgreen")
        plt.annotate("CO2", xy=[gas.mCO2[-1], P[-1]], color="saddlebrown")
        plt.annotate("CO", xy=[gas.mCO[-1], P[-1]], color="darkgreen")
        plt.annotate("CH4", xy=[gas.mCH4[-1], P[-1]], color="forestgreen")
    if (
        gas.sys.run.GAS_SYS == "SOH"
        or gas.sys.run.GAS_SYS == "COHS"
        or gas.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(gas.mSO2[1:], P, c="maroon")
        plt.plot(gas.mS2[1:], P, c="goldenrod")
        plt.plot(gas.mH2S[1:], P, c="pink")
        plt.annotate("SO2", xy=[gas.mSO2[-1], P[-1]], color="maroon")
        plt.annotate("S2", xy=[gas.mS2[-1], P[-1]], color="goldenrod")
        plt.annotate("H2S", xy=[gas.mH2S[-1], P[-1]], color="pink")
    if gas.sys.run.GAS_SYS == "COHSN":
        plt.plot(gas.mN2[1:], P, c="grey")
        plt.annotate("N2", xy=[gas.mN2[-1], P[-1]], color="grey")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel(f"Speciation in a {gas.sys.run.GAS_SYS} gas (mol frac)")
    plt.ylabel("Pressure (bar)")
    plt.savefig(path / "speciation(mol).png")
    plt.close()


# Gas speciation wt%


def plot_gasspecies_wt(gas, P, path):
    """Plots the gas speciation (weight fraction) vs pressure"""

    plt.plot(gas.Wt["H2O"], P, c="darkblue")
    plt.plot(gas.Wt["H2"], P, c="steelblue")
    plt.plot(gas.Wt["O2"], P, c="firebrick")
    plt.annotate("H2O", xy=[gas.Wt["H2O"][-1], P[-1]], color="darkblue")
    plt.annotate("H2", xy=[gas.Wt["H2"][-1], P[-1]], color="steelblue")
    plt.annotate("O2", xy=[gas.Wt["O2"][-1], P[-1]], color="firebrick")
    if (
        gas.sys.run.GAS_SYS == "COH"
        or gas.sys.run.GAS_SYS == "COHS"
        or gas.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(gas.Wt["CO2"], P, c="saddlebrown")
        plt.plot(gas.Wt["CO"], P, c="darkgreen")
        # plt.plot(gas.Wt['CH4'], P, c='mediumseagreen')
        plt.annotate("CO2", xy=[gas.Wt["CO2"][-1], P[-1]], color="saddlebrown")
        plt.annotate("CO", xy=[gas.Wt["CO"][-1], P[-1]], color="darkgreen")
        # plt.annotate('CH4', xy=[gas.mCH4[-1], 0.9], color='mediumseagreen')
    if (
        gas.sys.run.GAS_SYS == "SOH"
        or gas.sys.run.GAS_SYS == "COHS"
        or gas.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(gas.Wt["SO2"], P, c="maroon")
        plt.plot(gas.Wt["S2"], P, c="goldenrod")
        plt.plot(gas.Wt["H2S"], P, c="pink")
        plt.annotate("SO2", xy=[gas.Wt["SO2"][-1], P[-1]], color="maroon")
        plt.annotate("S2", xy=[gas.Wt["S2"][-1], P[-1]], color="goldenrod")
        plt.annotate("H2S", xy=[gas.Wt["H2S"][-1], P[-1]], color="pink")
    if gas.sys.run.GAS_SYS == "COHSN":
        plt.plot(gas.Wt["N2"], P, c="grey")
        plt.annotate("N2", xy=[gas.Wt["N2"][-1], P[-1]], color="grey")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel(f"Gas phase speciation of a {gas.sys.run.GAS_SYS} system (wt %)")
    plt.ylabel("Pressure (bar)")
    plt.savefig(path / "speciation(wt).png")
    plt.close()


# Melt volatile wt%


def plot_meltspecies(melt, P, path):
    """Plots the volatile content of the melt vs pressure"""

    plt.plot(melt.h2o, P, c="darkblue")
    plt.plot(melt.h2, P, c="steelblue")
    plt.annotate("H2O", xy=[melt.h2o[-1], P[-1]], color="darkblue")
    plt.annotate("H2", xy=[melt.h2[-1], P[-1]], color="steelblue")
    if (
        melt.sys.run.GAS_SYS == "COH"
        or melt.sys.run.GAS_SYS == "COHS"
        or melt.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(melt.co2, P, c="saddlebrown")
        plt.annotate("CO2", xy=[melt.co2[-1], P[-1]], color="saddlebrown")
    if (
        melt.sys.run.GAS_SYS == "SOH"
        or melt.sys.run.GAS_SYS == "COHS"
        or melt.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(melt.s, P, c="maroon")
        plt.annotate("S", xy=[melt.s[-1], P[-1]], color="maroon")
    if melt.sys.run.GAS_SYS == "COHSN":
        plt.plot(melt.n, P, c="grey")
        plt.annotate("N", xy=[melt.n[-1], P[-1]], color="grey")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel("Melt volatile content (wt%)")
    plt.ylabel("Pressure (bar)")
    plt.savefig(path / "meltspecies.png")
    plt.close()


# Exsolved gas mass and volume


def plot_gasfraction(sys, P, path):
    """Plots the gas volume fraction vs pressure"""

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.plot(cnvt.frac2perc(sys.WgT[1:]), P)
    ax1.invert_yaxis()
    ax1.set_xlabel("Exsolved gas wt%")
    ax1.set_ylabel("Pressure (bar)")
    ax2.plot(cnvt.frac2perc(sys.GvF), P)
    ax2.set_xlabel("Exsolved gas volume %")
    plt.savefig(path / "exsolved_gas.png")
