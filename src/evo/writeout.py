"""
Writes the results of the run to a file and
contains options to produce a graph of the results.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import evo.conversions as cnvt
import evo.plots as plots
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


def writeout_figs(df, filepath, out=None):
    """
    Selects the figures to produce at the end of a successful run.

    If an output file is not used to select the required figures, all
    options are plotted.

    Parameters
    ----------
    df:
        DataFrame containing the output data
    filepath : Path
        Path to the desired output folder
    out : Output class or None
        Active instance of the Output class, or None
    """

    filelist = list(filepath.glob("*.png"))
    # Removes previous files so if output specification is changed
    # there is no confusion as to up to date files.

    for file in filelist:
        file.unlink()

    plot_map = {
        "plot_melt_species": ("melt_species.png", plots.plot_meltspecies),
        "plot_gas_species_wt": ("speciation(wt).png", plots.plot_gasspecies_wt),
        "plot_gas_species_mol": ("speciation(mol).png", plots.plot_gasspecies_mol),
        "plot_gas_fraction": ("exsolved_gas.png", plots.plot_gasfraction),
        "plot_fo2_dFMQ": ("fo2_fmq.png", plots.plot_fo2FMQ),
    }

    chose_plots = {
        k: getattr(out, k) if out is not None else True for k in plot_map.keys()
    }

    for key, (filename, plot_func) in plot_map.items():
        if not chose_plots[key]:
            continue
        fig = plot_func(df)  # Call the plotting function
        fig.savefig(filepath / filename)
        plt.close(fig)  # Ensure the figure is closed to free memory
