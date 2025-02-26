"""
dgs.py  PL: UPDATE THESE DESCRIPTIONS

--- SPECIES LIST ---

C    CO   CO_2 CH_4    | carbon, carbon monoxide, carbon dioxide, methane
H_2  H_2O H_2S         | hydrogen, water, hydrogen sulfide
SO_2 S_2               | sulfur dioxide, sulfur
O_2  OCS               | oxygen, carbonyl sulfide
FeO  Fe_2O_3           | ferrous iron, ferric iron
N2                     | nitrogen

-- MODULES ---

density - calculate density of the melt
ferric - calculate fo2 from ferric iron proportions
fixed_weights - calculate the saturation point of a melt from it's volatile content
        given as atomic weights of species in the melt
init_cons - calculates the initial conditions in a standard run, given gas weight
             fraction and some subset of gas or melt properties
messages - stores warnings and error messages
multirun - a template for running parameters sweeps with EVo
sat_pressure - calculates the pressure of volatile saturation
solubility_laws - stores the solubility laws for all species
solver - calculates the final system composition at every pressure step
solvgas - calculate fugacity coefficients and equilibrium coefficients for gas phase
            equilibria

--- LITERATURE SOURCES ---

model description ---
Gaillard & Scaillet 2014, EPSL
Burgisser & Scaillet 2007
Moretti & Papale 2004
Holloway 1987
Burgisser 2015
Liggins et al 2020, 2021

solubility models ---
Iacono-marziano et al. 2012, EPSL
Gaillard & Scaillet 2009
Behrens and Gaillard 2006 (CO CH4 solubility)
Morizet et al. 2010 (CO CH4 solubility)
Hirschman 2012
Righter 2009


--- MODEL DESCRIPTION ---

Equilibrium constants and mass balance method
solve homogenoud and heterogeneous equilibria simultaneously at each step

"""

# ------------------------------------------------------------------------
# IMPORTS
# ------------------------------------------------------------------------

# python main [ ensure these are on your system ]
import time
from pathlib import Path

import numpy as np

import evo.conversions as cnvs
import evo.solver as solver

# bundled scripts
from evo.readin import readin, set_init_chem
from evo.writeout import writeout_figs, writeout_file

# ------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------


def run_evo(f_chem, f_env, f_out, folder="outputs"):
    """Main function for EVo.

    Call to run the model.

    Parameters
    ----------
    f_chem : str
        Path to the chemistry input file
    f_env : str
        Path to the environment input file
    f_out : str or NoneType
        Path to the file describing the required outputs, None if not used
    folder : str
        Path to the folder to write the results to
    """
    start = time.time()

    print("Reading in from:")
    print("Chemistry file:", f_chem)
    print("Environment file:", f_env)
    print("Output file:", f_out, "\n")

    # Instantiate the run, thermosystem, melt and output objects
    run, sys, melt, out = readin(f_chem, f_env, f_out)
    run.results_folder = Path(folder)

    print("Set parameters:")
    run.print_conditions()

    # Setup the initial gas + melt chemistry, calculate the atomic masses of elements
    run, sys, melt, gas, molecules = set_init_chem(run, sys, melt)

    print("\nSystem chemistry (wt%)")
    sm = 0.0
    for k, v in cnvs.mol2wt(melt.cm_dry).items():
        print(k, ":", f"{v * 100:0.3g}")
        sm += v * 100
    print("Total =", round(sm, 2), "\n")

    # print output specifications if made
    if f_out:
        out.print_conditions()

    # -------------------------------------------------------------------
    # CALCULATION LOOP --------------------------------------------------
    # -------------------------------------------------------------------

    if run.SINGLE_STEP is True:
        solver.decompress(run, sys, melt, gas, molecules)
        sys.P_track.append(sys.P)

        print(
            f"The total gas weight percent is {sys.WgT[-1] * 100:.2f}% and "
            f"the gas volume fraction is {sys.GvF[-1] * 100:.2f}%"
        )
        print(
            f"The gas phase at {sys.P:.0f} bar is partitioned according to molar "
            "fraction as: \n "
            f"O2: {gas.mO2} \n "
            f"H2: {gas.mH2} \n "
            f"H2O: {gas.mH2O} \n "
            f"CO2: {gas.mCO2} \n "
            f"CO: {gas.mCO} \n "
            f"CH4: {gas.mCH4} \n "
            f"S2: {gas.mS2} \n "
            f"SO2: {gas.mSO2} \n "
            f"H2S: {gas.mH2S} \n "
            f"N2: {gas.mN2} \n"
        )
        writeout_file(sys, gas, melt, sys.P_track)

    elif run.SINGLE_STEP is False:
        dp = (
            int(abs(np.floor(np.log10(np.abs(run.P_STOP))))) + 1
            if run.P_STOP < 1e-4
            else 5
        )
        while np.round(sys.P, decimals=dp) >= run.P_STOP:
            sys.P = np.round(sys.P, decimals=dp)
            sys.P_track.append(sys.P)
            solver.decompress(run, sys, melt, gas, molecules)  # does one pressure step
            sys.pressure_step()

        if (
            run.FIND_SATURATION is True or run.ATOMIC_MASS_SET is True
        ) and sys.sat_conditions[0] < run.P_STOP:
            exit(
                f"Error: The saturation pressure ({sys.sat_conditions[0]:.3} bar) is "
                f"lower than P_STOP ({run.P_STOP:.3} bar).\n"
                "Please either lower P_STOP or change the melt initial fo2/melt "
                "volatile content.\n"
                "Exiting..."
            )

        print(
            f"The pressure is {run.P_STOP:.2f} bar, "
            f"the total gas weight percent is {sys.WgT[-1] * 100:.2f}% "
            f"and the gas volume fraction is {sys.GvF[-1] * 100:.2f}% "
            "\nThe run has finished."
        )

        print(
            "The gas is partitioned as (mol%): \n "
            f"O2: {gas.mO2[-1] * 100:.4e} \n "
            f"H2: {gas.mH2[-1] * 100:.4f} \n "
            f"H2O: {gas.mH2O[-1] * 100:.4f} \n "
            f"CO2: {gas.mCO2[-1] * 100:.4f} \n "
            f"CO: {gas.mCO[-1] * 100:.4f} \n "
            f"CH4: {gas.mCH4[-1] * 100:.4e} \n "
            f"S2: {gas.mS2[-1] * 100:.4f} \n "
            f"SO2: {gas.mSO2[-1] * 100:.4f} \n "
            f"H2S: {gas.mH2S[-1] * 100:.4f} \n "
            f"N2: {gas.mN2[-1] * 100:.4f} \n"
        )

        print(
            "The first, middle, last gas molar fractions are: \n "
            f"O2: {gas.mO2[0]}   {gas.mO2[(len(gas.mO2) // 2)]}   {gas.mO2[-1]}\n"
            f"H2: {gas.mH2[0]}   {gas.mH2[(len(gas.mH2) // 2)]}   {gas.mH2[-1]}\n"
            f"H2O: {gas.mH2O[0]}   {gas.mH2O[(len(gas.mH2O) // 2)]}   {gas.mH2O[-1]}\n"
            f"CO2: {gas.mCO2[0]}   {gas.mCO2[(len(gas.mCO2) // 2)]}   {gas.mCO2[-1]}\n"
            f"CO: {gas.mCO[0]}   {gas.mCO[(len(gas.mCO) // 2)]}   {gas.mCO[-1]}\n"
            f"CH4: {gas.mCH4[0]}   {gas.mCH4[(len(gas.mCH4) // 2)]}   {gas.mCH4[-1]}\n"
            f"S2: {gas.mS2[0]}   {gas.mS2[(len(gas.mS2) // 2)]}   {gas.mS2[-1]}\n"
            f"SO2: {gas.mSO2[0]}   {gas.mSO2[(len(gas.mSO2) // 2)]}   {gas.mSO2[-1]}\n"
            f"H2S: {gas.mH2S[0]}   {gas.mH2S[(len(gas.mH2S) // 2)]}   {gas.mH2S[-1]}\n"
            f"N2: {gas.mN2[0]}   {gas.mN2[(len(gas.mN2) // 2)]}   {gas.mN2[-1]}\n"
        )

        end = time.time()
        print("Run time is ", end - start)

        df = writeout_file(sys, gas, melt, sys.P_track)
        writeout_figs(df, run.results_folder, out=out)

        return df
