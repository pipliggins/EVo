# readin.py

# ------------------------------------------------------------------------
# IMPORTS
# ------------------------------------------------------------------------

import numpy as np
import yaml
from numpy import exp, log10
from numpy import log as ln

import evo.conversions as cnvs
import evo.fixed_weights as fw
import evo.messages as msgs
import evo.sat_pressure as sat
import evo.solvgas as sg

# bundled scripts
from evo.dgs_classes import Gas, Melt, Molecule, Output, RunDef, ThermoSystem

# ------------------------------------------------------------------------
# FUNCTION DEFINTIONS
# ------------------------------------------------------------------------


def readin_env(f):
    """
    Instantiates the RunDef and ThermoSystem classes.

    Sets the parameters of the RunDef and ThermoSystem classes according
    to the setup parameters listed in the environment input file.

    Parameters
    ----------
    f : file object
        The open environment file

    Returns
    -------
    run : RunDef class
        The instantiated RunDef class
    sys : ThermoSystem class
        The instantiated ThermoSystem class
    """

    # creates a dictionary of the environment parameters from env.yaml
    x = yaml.full_load(f)

    # setup run definitions
    run = RunDef()
    run.param_set(x)  # calls the param_set function from the RunDef() class

    # use run definitions to initialise system state
    sys = ThermoSystem(run)

    if run.FO2_buffer_SET is True:
        if run.FO2_buffer == "FMQ":
            sys.FO2_buffer = run.FO2_buffer_START
            sys.FO2 = ln(
                10 ** cnvs.fmq2fo2(sys.FO2_buffer, sys.T, sys.P, run.FMQ_MODEL)
            )  # ln fO2.

        elif run.FO2_buffer == "IW":
            sys.FO2_buffer = run.FO2_buffer_START
            sys.FO2 = ln(10 ** cnvs.iw2fo2(sys.FO2_buffer, sys.T, sys.P, run.FMQ_MODEL))

        elif run.FO2_buffer == "NNO":
            sys.FO2_buffer = run.FO2_buffer_START
            sys.FO2 = ln(
                10 ** cnvs.nno2fo2(sys.FO2_buffer, sys.T, sys.P, run.FMQ_MODEL)
            )

        run.FO2_SET = True

    elif run.FO2_SET is True:
        sys.FO2 = np.log(run.FO2_START)
        sys.FO2_buffer = cnvs.fo2_2fmq(log10(exp(sys.FO2)), sys.T, sys.P, run.FMQ_MODEL)

    if run.GRAPHITE_SATURATED is True:
        if run.C_MODEL == "eguchi2018":
            pass
        else:
            exit(
                "Error: If melt is graphite saturated, "
                "please use the eguchi2018 model as the C_MODEL setting."
            )

    sys.gas_system(run)  # creates the molecule list (self.SC)

    return run, sys


#  -------------------------------------------------------------------------
def readin_chem(f, run, sys):
    """
    Instantiates the Melt class, reading in the chemistry file.

    Sets the parameters of the Melt class according to the chemistry
    input file. Checks that only recognised species are listed in the
    melt composition, and that all the information required to run the
    requested setup is provided.

    Parameters
    ----------
    f : file object
        The open chemistry file
    run : RunDef class
        Active instance of the RunDef class
    sys : ThermoSystem class
        Active instance of the ThermoSystem class

    Returns
    -------
    melt : Melt class
        The instantiated Melt class
    """

    data = yaml.full_load(f)

    ele_names = []
    chems = []

    # define allowed components
    allowed_names = [
        "sio2",
        "tio2",
        "al2o3",
        "feo",
        "fe2o3",
        "mno",
        "mgo",
        "cao",
        "na2o",
        "k2o",
        "p2o5",
        "nio",
        "cr2o3",
    ]

    for chem, frac in data.items():
        if chem.lower() in allowed_names:
            ele_names.append(chem.lower())
            chems.append(frac)
        else:
            exit(
                f"{chem} does not match an element in the allowed list "
                "(case insensitive):\n{allowed_names}"
            )

    # test for iron/ fO2 definition
    if "feo" in ele_names and "fe2o3" not in ele_names and run.FO2_SET is not True:
        if run.GAS_SYS == "OH" and (run.FH2_SET):
            pass
        elif run.FH2_SET and (run.FH2O_SET or run.WTH2O_SET):
            pass
        else:
            exit("Error: Only FeO given without setting dgs_FO2 = True.")

    elif "fe2o3" in ele_names and "feo" not in ele_names and run.FO2_SET is not True:
        if run.GAS_SYS == "OH" and (run.FH2_SET):
            pass
        elif run.FH2_SET and (run.FH2O_SET or run.WTH2O_SET):
            pass
        else:
            exit("Error: Only Fe2O3 given without setting dgs_FO2 = True.")

    elif "feo" not in ele_names and "fe2o3" not in ele_names:
        exit("Error: No iron in system")

    elif run.FO2_SET is True and "feo" in ele_names and "fe2o3" in ele_names:
        exit("Error: fO2 and iron proportions specified, only give one.")

    # populate melt object
    melt = Melt(run, sys)

    melt.chem_init(ele_names, chems)

    return melt


#  -------------------------------------------------------------------------
def readin_output(f):
    """
    Instantiates the Output class using the output file.

    Parameters
    ----------
    f : file object
        The open environment file

    Returns
    -------
    out : Output class
        The instantiated Output class
    """

    x = yaml.full_load(f)

    out = Output()
    out.set_outputs(x)

    return out


#  ------------------------------------------------------------------------
def run_melt_match(run, melt):
    """
    Asserts that the composition name given matches the melt SiO2 content.

    Parameters
    ----------
    run : RunDef class
        Active instance of the RunDef class
    melt : Melt class
        Active instance of the Melt class
    """

    if run.COMPOSITION == "basalt":
        sio2_lim = [45, 55]
    elif run.COMPOSITION == "phonolite":
        sio2_lim = [52, 63]
    if run.COMPOSITION == "rhyolite":
        sio2_lim = [65, 80]

    magma = melt.Cw()

    if magma["sio2"] > sio2_lim[0] and magma["sio2"] < sio2_lim[1]:
        pass
    else:
        # Creates a terminal message checking if they want to continue with their setup.
        msgs.run_chem_mismatch(run.COMPOSITION, sio2_lim, magma["sio2"])


#  -------------------------------------------------------------------------
def readin(f_chem, f_env, *args):
    """
    Opens the input files and sets up major classes with the input data.

    Checks that correct input data has been provided for the requested
    run type.

    Parameters
    ----------
    f_chem : string
        Path to the chemistry input file
    f_env : string
        Path to the environment input file
    *args : string
        optional path to the output options input file

    Returns
    -------
    run : RunDef class
        The active instance of the RunDef class
    sys : ThermoSystem class
        The active instance of the ThermoSystem class
    melt : Melt class
        The active instance of the Melt class
    out : Output class
        The active instance of the Output class
    """

    # environment
    with open(f_env) as f:
        run, sys = readin_env(f)

    # check stepsize requirements
    if run.DP_MIN > run.DP_MAX:
        raise ValueError(
            "Please make the minimum stepsize DP_MIN <= maximum stepsize DP_MAX."
        )

    if run.RUN_TYPE == "open" and run.DP_MAX != run.DP_MIN:
        raise ValueError(
            "Open system degassing is path dependent.\n"
            "For internal consistency, please set the DP_MAX and DP_MIN "
            "pressure steps to be equal (we suggest <1 bar for the run to complete; "
            "0.5 bar is usually sufficient)"
        )

    # chemistry
    with open(f_chem) as f:
        melt = readin_chem(f, run, sys)

    # check the melt SiO2 content matches the run composition definition
    run_melt_match(run, melt)

    if run.FIND_SATURATION is False:
        msgs.vol_setup_standard(run, sys)
    elif run.FIND_SATURATION is True:
        msgs.vol_setup_saturation(run, sys)

    # outputs
    for arg in args:
        if arg:
            with open(arg) as f:
                out = readin_output(f)
        else:
            out = None

    return run, sys, melt, out


# -------------------------------------------------------------------------
def set_init_chem(run, sys, melt):
    """
    Wrapper function to setup initial conditions of the volcanic system.

    Instantiates the Molecule class for each of the gas phase species.
    Calculates the equilibrium constants for the system temperature, the
    initial gas gas composition, the starting pressure (if the saturation
    pressure is requested) and sets the atomic volatile masses for mass
    balance conservation. Normalises the melt oxide composition with the
    volatile content and records the mass of total Fe in the system.

    Parameters
    ----------
    run : RunDef class
        The active instance of the RunDef class
    sys : ThermoSystem class
        The active instance of the ThermoSystem class
    melt : Melt class
        The active instance of the Melt class

    Returns
    -------
    run : RunDef class
        The active instance of the RunDef class
    sys : ThermoSystem class
        The active instance of the ThermoSystem class
    melt : Melt class
        The active instance of the Melt class
    gas : Gas class
        The active instance of the Gas class
    molecules : tuple of Molecule classes
        A tuple of all the Molecule class instances. Setup to be passed
        through all functions, preserving ordering.
    """

    # Instantiates the Molecule class with each of the species that should be present
    # depending on the gas system which the user has specified.

    H2O = Molecule(run, sys, melt, "H2O")
    O2 = Molecule(run, sys, melt, "O2")
    H2 = Molecule(run, sys, melt, "H2")

    if run.GAS_SYS == "OH":
        molecules = H2O, O2, H2

    elif run.GAS_SYS == "SOH":
        S2 = Molecule(run, sys, melt, "S2")
        SO2 = Molecule(run, sys, melt, "SO2")
        H2S = Molecule(run, sys, melt, "H2S")

        molecules = H2O, O2, H2, S2, SO2, H2S

    else:
        if "C" in run.GAS_SYS:
            CO = Molecule(run, sys, melt, "CO")
            CO2 = Molecule(run, sys, melt, "CO2")
            CH4 = Molecule(run, sys, melt, "CH4")

            molecules = H2O, O2, H2, CO, CO2, CH4

            if "S" in run.GAS_SYS:
                S2 = Molecule(run, sys, melt, "S2")
                SO2 = Molecule(run, sys, melt, "SO2")
                H2S = Molecule(run, sys, melt, "H2S")

                molecules = H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S

                if "N" in run.GAS_SYS:
                    N2 = Molecule(run, sys, melt, "N2")

                    molecules = H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2

    gas = Gas(sys)

    sys.K = sg.get_K(sys, molecules)  # Calculates and stores all relevant K values

    sg.set_Y(sys, molecules)  # Gets all relevant activity coefficients

    # Finds the saturation pressure for the melt composition given,
    # sets sys.P to this and WgT to 0.001%
    if run.FIND_SATURATION is True:
        sys.sat_conditions = sat.sat_pressure(run, sys, gas, melt, molecules)

    # Finds the saturation pressure for the melt composition given,
    # sets sys.P to this and WgT to 0.001%
    elif run.ATOMIC_MASS_SET is True:
        sys.sat_conditions = fw.sat_pressure(run, sys, gas, melt, molecules)

    else:
        sys.get_atomic_mass(
            run, gas, melt, molecules
        )  # Gets atomic mass values and initial partitioning

    # Normalises the wt % fractions of oxides in the melt with the volatile species
    # and FeO(t).
    sys.norm_with_gas(melt)

    # Gets the total WtO after the melt (and therefore Fe) has been normalised to
    # the gas content
    sys.get_wto_tot(melt)

    return run, sys, melt, gas, molecules
