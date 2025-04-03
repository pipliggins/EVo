# dgs_classes

import copy
import importlib.resources
import inspect
import itertools
import re
import sys
from pathlib import PosixPath
from typing import Literal, Optional

import numpy as np
from pydantic import BaseModel, field_validator, model_validator

import evo.constants as cnst
import evo.conversions as cnvs
import evo.density as density
import evo.init_cons as ic
import evo.messages as msgs
import evo.solubility_laws as sl
import evo.solvgas as sg


class RunDef(BaseModel):
    """
    Stores the parameters of the degassing run

    Attributes
    ----------
    COMPOSITION : {'basalt', 'phonolite', 'rhyolite'}
        Name of magma composition
    RUN_TYPE : {'closed', 'open'}
        Run type - closed keeps gas and melt in contact, open removes an
        aliquot of the gas after every pressure step
    SINGLE_STEP : bool, default = False
        Choose whether to calculate the composition at a single pressure
    FIND_SATURATION : bool, default = False
        Choose to calculate the volatile saturation point rather than set
        a starting pressure
    GAS_SYS : {'OH', 'COH', 'SOH', 'COHS', 'COHSN'}
        Sets the volatile elements of interest
    FE_SYSTEM : bool, default = True
        If True, oxygen is exchanged between iron in the melt and the gas
    S_SAT_WARN : bool, default = True
        If true, stops a run is the melt becomes sulfide saturated.
    T_START : float
        System temperature, K
    P_START : float
        Starting pressure, bar
    DP_MIN : float
        Minimum pressure step size, bar. Must be equal to `DP_MAX` if
        `RUN_TYPE` = 'open'
    DP_MAX : float
        Maximum pressure step size, bar. Must be equal to `DP_MIN` if
        `RUN_TYPE` = 'open'
    MASS : float
        System mass, grams
    WgT : float
        Initial gas mass fraction (not percent). Not required if either
        `FIND_SATURATION` or `ATOMIC_MASS_SET` = True.
    LOSS_FRAC : float, <1
        The fraction of the total gas in the system to be removed if
        `RUN_TYPE` = 'open'.
    DENSITY_MODEL : {'spera2000'}
        Select a model to calculate the magma density with.
    FO2_MODEL : {'kc1991', 'r2013'}
        Select the model to calculate the relationship between the ferric/
        ferrous ratio and the melt fO2.
    FMQ_MODEL : {'frost1991'}
        Select the model to calculate fO2 mineral buffers with
    H2O_MODEL : {'burguisser2015'}
        Select the H2O solubility model
    H2_MODEL : {'gaillard2003', 'burguisser2015'}
        Select the H2 solubility model
    C_MODEL : {'burguisser2015', 'eguchi2018'}
        Select the CO2 solubility model
    CO_MODEL : {None, 'armstrong2015'}
        Select the CO solubility model
    CH4_MODEL : {None, 'ardia2013'}
        Select the CH4 solubility model
    SULFIDE_CAPACITY : {'oneill2020', 'oneill2002'}
        Select the sulfide capacity model
    SULFATE_CAPACITY : {'nash2019'}
        Select the sulfate capacity model
    SCSS : {'liu2007'}
        Select the model to calculate the sulfide capacity at sulfide
        saturation
    N_MODEL : {'libourel2003'}
        Select the N2 solubility model
    FO2_buffer : {'FMQ', 'IW', 'NNO'}
        Mineral buffer to cite fO2 relative to
    FO2_buffer_START : float
        System fO2 relative to `FO2_buffer`.
    FO2_START : float
        Set the absolute fO2. Do not use in conjunction with `FO2_buffer_START`.
    FH2_START : float
        Set the H2 fugacity as a starting condition (bar)
    FH2O_START : float
        Set the H2O fugacity as a starting condition (bar)
    FCO2_START : float
        Set the CO2 fugacity as a starting condition (bar)
    ATOMIC_MASS_SET : bool, default False
        Choose to set the mass of each volatile element in the system and calculate the
        volatile saturation point, rather than set a starting pressure
    ATOMIC_H, ATOMIC_C, ATOMIC_S, ATOMIC_N : float
        Mass fraction (parts per million) of the corresponding element in
        the system. Only use alongside `ATOMIC_MASS_START`.
    WTH2O_START, WTCO2_START, SULFUR_START, NITROGEN_START : float
        Initial weight FRACTION (not %) of the species dissolved in the
        melt at the starting pressure.
    GRAPHITE_SATURATED  : bool, default = False
        If True, the melt starts graphite saturated, with a mass of
        graphite specified in `GRAPHITE_START`.
    GRAPHITE_START : float
        Initial melt graphite content as a weight fraction, use if
        `GRAPHITE_SATURATED` is True
    """

    COMPOSITION: Literal["basalt", "phonolite", "rhyolite"] = "basalt"
    RUN_TYPE: Literal["closed", "open"] = "closed"
    SINGLE_STEP: bool = False
    FIND_SATURATION: bool = False
    GAS_SYS: str = "OH"
    FE_SYSTEM: bool = True
    OCS: bool = False
    S_SAT_WARN: bool = True

    # Physical parameters
    T_START: float = 1473.15
    P_START: float = 3000.0
    P_STOP: float = 1.0
    DP_MIN: float = 1.0
    DP_MAX: float = 100.0
    MASS: float = 100.0
    WgT: float = 0.001
    LOSS_FRAC: float = 0.99

    # Model options
    DENSITY_MODEL: Literal["spera2000"] = "spera2000"
    FO2_MODEL: Literal["kc1991", "r2013"] = "kc1991"
    FMQ_MODEL: Literal["frost1991"] = "frost1991"
    H2O_MODEL: Literal["burguisser2015"] = "burguisser2015"
    H2_MODEL: Literal["gaillard2003", "burguisser2015"] = "gaillard2003"
    C_MODEL: Literal["burguisser2015", "eguchi2018"] = "burguisser2015"
    CO_MODEL: Optional[Literal["armstrong2015"]] = None
    CH4_MODEL: Optional[Literal["ardia2013"]] = None
    SULFIDE_CAPACITY: Literal["oneill2020", "oneill2002"] = "oneill2020"
    SULFATE_CAPACITY: Literal["nash2019"] = "nash2019"
    SCSS: Literal["liu2007"] = "liu2007"
    N_MODEL: Literal["libourel2003"] = "libourel2003"

    # fO2 and fugacity controls
    FO2_buffer_SET: bool = False
    FO2_buffer: Literal["FMQ", "IW", "NNO"] = "FMQ"
    FO2_buffer_START: float = 0.0
    FO2_SET: bool = False
    FO2_START: float = 0.0
    FH2_SET: bool = False
    FH2_START: float = 0.24
    FH2O_SET: bool = False
    FH2O_START: float = 1000.0
    FCO2_SET: bool = False
    FCO2_START: float = 0.01

    # Atomic masses - ppm
    ATOMIC_MASS_SET: bool = False
    ATOMIC_H: float = 550.0
    ATOMIC_C: float = 200.0
    ATOMIC_S: float = 4000.0
    ATOMIC_N: float = 10.0

    # Initial melt volatile concentrations
    WTH2O_SET: bool = False
    WTH2O_START: float = 0.03
    WTCO2_SET: bool = False
    WTCO2_START: float = 0.01
    SULFUR_SET: bool = False
    SULFUR_START: float = 0.001
    NITROGEN_SET: bool = False
    NITROGEN_START: float = 0.01
    GRAPHITE_SATURATED: bool = False
    GRAPHITE_START: float = 0.0

    # Optional output parameter
    results_folder: Optional[PosixPath] = None

    @model_validator(mode="after")
    def validate_config(self) -> "RunDef":
        if self.OCS:
            raise ValueError(
                "Error: OCS functionality has not been released yet. Sorry!"
            )

        if self.MASS <= 0:
            raise ValueError("The system has no mass value")

        if self.WgT <= 0 and not self.FIND_SATURATION and not self.ATOMIC_MASS_SET:
            raise ValueError(
                "A gas weight fraction of greater than 0 must be entered "
                "if the saturation pressure is not being calculated."
            )

        return self

    @field_validator("GAS_SYS")
    @classmethod
    def validate_gas_sys(cls, v: str) -> str:
        ALLOWED_VARIANTS = ["OH", "COH", "SOH", "COHS", "COHSN"]
        # Normalize: uppercase, sort letters
        normalized = "".join(sorted(v.upper()))
        for allowed in ALLOWED_VARIANTS:
            if normalized == "".join(sorted(allowed)):
                return allowed  # Return canonical form (e.g. "COHS")
        raise ValueError(
            f"Invalid GAS_SYS: {v}. Must be a permutation of one of {ALLOWED_VARIANTS}"
        )

    def print_conditions(self):
        print("".join(f"{k} = {v}\n" for k, v in self.model_dump().items()))


#  -------------------------------------------------------------------------
class ThermoSystem:
    """
    Thermodynamic properties of the whole system (gas and melt).

    Attributes
    ----------

    run : RunDef class
        Active instance of the RunDef class
    T : float
        System temperature (K), inherited from RunDef
    P : float
        Current pressure (bar)
    P_step : float
        The current size of the pressure step to take (bar)
    P_track : list of floats
        Stores the pressures at which solutions have been found
    M : float
        The system mass, inherited from RunDef
    rho : list of floats
        The current bulk density of the melt+gas at every step
    sat_conditions : tuple
        Stores the saturation conditions (P_sat, gas composition, fugacity
        coefficients, melt content)
    graph_unsat_rerun : bool, default = False
        Records whether a re-run of a pressure is being performed because
        the melt became graphite under-saturated.
    FO2 : float
        current absolute fO2 STORED AS LN(FO2)
    FO2_buffer : float
        Stores the fO2 of the system relative to the mineral buffer stored
        in RunDef.
    SC : list of str
        Stores the system chemistry - the gas species being modelled.
    K : dict
        Stores the the equilibrium constants for gas phase equilibria
    atomicM : dict
        Masses of each volatile element - weight fractions.
    WgT : list of floats
        stores the total gas weight fraction at every pressure step.
    GvF : list of floats
        Stores the gas volume fraction at every pressure step
    Cw : dict
        Stores the SYSTEM mass fractions of oxides and volatile elements,
        normalised to 1. Iron stored as Fe so oxygen isn't counted twice.
    """

    def __init__(self, run):
        self.run = run
        self.T = self.run.T_START  # K
        self.P = self.run.P_START  # bar
        self.P_step = self.run.DP_MAX  # bar
        self.P_track = []  # stores the pressure steps run at
        self.M = self.run.MASS  # g
        self.rho = []  # bulk density
        self.sat_conditions = 0  # store saturation pressure + gas speciation
        self.graph_unsat_rerun = False

        # intialise for subsequent setting
        self.FO2 = None  # fo2 # STORED AS A LN(FO2) SO DO e**SELF.FO2 WHERE IT'S NEEDED
        self.FO2_buffer = None
        self.SC = []  # System chemistry - i.e. gas species present.
        self.K = {}  # Stores K values (equilibrium constants)
        self.atomicM = {}  # Atomic masses for mass balance - weight fraction
        self.WgT = [self.run.WgT]  # Total gas wt%, initially set by user.
        self.GvF = []  # Stores the gas volume fraction at each step
        self.Cw = {}
        self.OCS = False  # Can be flipped after the last pressure step for a final run.

    @property
    def Ppa(self):
        """Converts system pressure in bar, to pressure in pascal (Pa)"""
        return self.P * 1e5

    # Selects the system chemistry from the system the user wants to examine.
    def gas_system(self, run):
        """
        Selects the system chemistry given the elements the user wishes to track.

        Parameters
        ----------
        run : RunDef class
            Active instance of the RunDef class

        Returns
        -------
        self.SC : list of str
            List of all the gas phase species.
        """

        run.GAS_SYS = run.GAS_SYS.upper()

        def perm(string, length):
            return list(map("".join, itertools.permutations(string, length)))

        OH = perm("OH", 2)
        COH = perm("COH", 3)
        SOH = perm("SOH", 3)
        COHS = perm("COHS", 4)
        COHSN = perm("COHSN", 5)

        if run.GAS_SYS in OH:
            run.GAS_SYS = "OH"
            self.SC = ["H2O", "O2", "H2"]
        elif run.GAS_SYS in COH:
            run.GAS_SYS = "COH"
            self.SC = ["H2O", "O2", "H2", "CO2", "CO", "CH4"]
        elif run.GAS_SYS in SOH:
            run.GAS_SYS = "SOH"
            self.SC = ["H2O", "O2", "H2", "S2", "H2S", "SO2"]
        elif run.GAS_SYS in COHS:
            run.GAS_SYS = "COHS"
            self.SC = ["H2O", "O2", "H2", "CO2", "CO", "CH4", "S2", "H2S", "SO2"]
        elif run.GAS_SYS in COHSN:
            run.GAS_SYS = "COHSN"
            self.SC = ["H2O", "O2", "H2", "CO2", "CO", "CH4", "S2", "H2S", "SO2", "N2"]
        return self.SC

    def get_atomic_mass(self, run, gas, melt, mols):
        """Finds the masses of volatile elements in the system."""

        if run.GAS_SYS == "OH":
            ic.oh(self, run, melt, gas, mols)
        elif run.GAS_SYS == "COH":
            ic.coh(self, run, melt, gas, mols)
        elif run.GAS_SYS == "SOH":
            ic.soh(self, run, melt, gas, mols)
        elif run.GAS_SYS == "COHS":
            ic.cohs(self, run, melt, gas, mols)
        elif run.GAS_SYS == "COHSN":
            ic.cohsn(self, run, melt, gas, mols)

    def norm_with_gas(self, melt):  # for iron equilibration.
        """
        Normalises the melt with the gas phase ready to equilibrate the oxygen.

        Updates the Cw attribute of ThermoSystem to be oxides + gases
        normalised to 1. Finds the ferric/ferrous ratio and adds this
        to melt.F

        Parameters
        ----------
        melt : Melt class
            Active instance of the Melt class
        """

        cons = 0  # holds constant volatile species and iron content
        C_s = {}  # Holds the values to be normalised
        C_w = melt.Cw()  # wt %

        for key in self.atomicM:
            self.Cw[key] = self.atomicM[key]
            cons += self.atomicM[key]

        # The oxide composition of the melt partitioned as wt %
        # PL: This doesn't convert fe2o3 to feo at the moment, it probably should?
        # Need to catch scenario where only fe2o3 is provided?
        for ele in C_w:
            if ele == "feo" or ele == "fe2o3":
                cons += C_w[ele] / 100
                self.Cw[ele] = C_w[ele] / 100
            else:
                # oxide weights only, without volatiles or iron, ready to be scaled.
                C_s[ele] = C_w[ele] / 100

        # returns the 9 major oxides as weight %, with feoT and volatiles held constant.
        oxides = cnvs.norm(C_s, (1 - cons))

        self.Cw.update(
            oxides
        )  # adds the normalised oxide values to the total system mass frac dictionary

        melt.cm_dry, F = melt.iron_fraction(
            self.FO2
        )  # this is now a duplicate of the one inside init cons, hopefully.
        melt.F.append(F)

    def get_wto_tot(self, melt):
        """
        Calculates the total amount of O: volatile system + iron.

        Finds the mass fraction of elemental iron, O stored in Fe within
        the melt and total O in the system. Updates self.Cw so it has
        feo/fe2o3 ratio.

        Parameters
        ----------
        melt : Melt class
            Active instance of the Melt class

        Returns
        -------
        self.atomicM : dict
            Dict of the volatile element masses, updated with oxygen values
        """

        Cm = cnvs.wt2mol(self.Cw)

        for k, v in melt.cm_dry.items():
            Cm[k] = v

        Cm = cnvs.norm(Cm, cons=["o", "h", "c", "s", "n"])

        Cm["fe"] = Cm["feo"] + 2 * Cm["fe2o3"]
        m_ofe = (
            Cm["feo"] + 3 * Cm["fe2o3"]
        )  # The mole frac of oxygen tied up in the melt.

        Cm.pop("feo", None)
        Cm.pop("fe2o3", None)

        Cm["ogas"] = Cm["o"]
        Cm["ofe"] = m_ofe
        Cm.pop("o", None)

        Cm = cnvs.norm(Cm, cons=["fe", "ogas", "ofe", "h", "c", "s", "n"])
        self.Cw = dict(cnvs.mol2wt(Cm))

        self.atomicM["fe"] = self.Cw["fe"]
        self.atomicM["ofe"] = self.Cw["ofe"]
        self.atomicM["o_tot"] = self.atomicM["o"] + self.Cw["ofe"]

        return self.atomicM

    def fe_equil(self, melt, mo2, O2):
        """
        Calculate the amount of oxygen in the melt for use in solver.

        Parameters
        ----------
        melt : Melt class
            Active instance of the Melt class
        mo2 : float
            Current mol fraction of O2 in the gas phase
        O2 : Molecule class
            O2 instance of the Molecule class

        Returns
        -------
        ofe : float
            Mass of volatile elemental O in the melt.
        """

        fo2 = np.log(O2.Y * mo2 * self.P)

        # returns anhydrous silicate melt comp mol fraction, Fe2/Fe3.
        F = melt.iron_fraction(fo2)[1]

        ofe = (
            cnst.m["o"]
            * (self.atomicM["fe"] / cnst.m["fe"])
            * ((1 + 3 * F) / (1 + 2 * F))
        )  # THIS STEP total weight fraction of atomic O in iron

        return ofe

    def fe_save(self, melt, mo2, O2):
        """
        Saveout the amount of O in the gas and melt.

        After the melt is in redox equilibrium with the gas phase, the
        amount of oxygen and the ferric/ferrous ratio is recorded.

        Parameters
        ----------
        melt : Melt class
            Active instance of the Melt class
        mo2 : float
            mole fraction of O2 in the gas phase
        O2 : Molecule class
            O2 instance of the Molecule class
        """

        fo2 = np.log(O2.Y * mo2 * self.P)

        # returns anhydrous silicate melt comp mol fraction, Fe2/Fe3.
        F = melt.iron_fraction(fo2)[1]

        ofe = (
            cnst.m["o"]
            * (self.atomicM["fe"] / cnst.m["fe"])
            * ((1 + 3 * F) / (1 + 2 * F))
        )  # THIS STEP total weight fraction of atomic O in iron

        self.atomicM["ofe"] = ofe
        self.atomicM["o"] = self.atomicM["o_tot"] - ofe
        melt.ofe.append(ofe)
        melt.F.append(F)

    def mass_conservation_reset(self, melt, gas):
        """
        Deletes results of current run to redo at a higher pressure.

        Called if the element masses are no longer being conserved. Deletes
        the results of the last run and reduces the step size to help
        convergence towards the correct solution, ready to run again with
        the same initial conditions at a higher pressure.

        Parameters
        ----------
        melt : Melt class
            Active instance of the Melt class
        gas : Gas class
            Active instance of the Gas class
        """

        sys_lsts = [self.WgT, self.GvF, self.rho]
        gas_lsts = [
            gas.fo2,
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
            gas.M,
        ]
        melt_lsts = [
            melt.fmq,
            melt.h2o,
            melt.h2,
            melt.co2,
            melt.co,
            melt.ch4,
            melt.graphite,
            melt.sulfide,
            melt.sulfate,
            melt.s,
            melt.n,
            melt.F,
            melt.ofe,
        ]
        gas_wts_f = ["H2O", "O2", "H2", "CO", "CO2", "CH4", "S2", "SO2", "H2S", "N2"]

        for x in sys_lsts + gas_lsts + melt_lsts:
            del x[-1]

        for x in gas_wts_f:
            del gas.Wt[x][-1]
            if x != "O2":
                del gas.f[x][-1]

        if melt.graphite_sat is True:
            melt.graph_current = melt.graphite[-1]

        if self.P_step / 10 >= self.run.DP_MIN:
            del self.P_track[-1]
            self.P = self.P + self.P_step
            self.P_step = self.P_step / 10
            print("mass conservation reset")
        else:
            del self.P_track[-1]
            msgs.closed_earlyexit(self, gas, melt)
            exit(
                "Error: Mass is no longer being conserved. "
                "Please reduce the minimum pressure stepsize and try again."
                "\nExiting..."
            )

    def pressure_step(self):
        """Reduces the pressure step size by x10 as it approaches `P_STOP`"""

        # If taking another pressure step won't take you below P_STOP, continue.
        if self.P - self.P_step > self.run.P_STOP:
            self.P -= self.P_step

        # If taking another step puts you below the min pressure and step/10 is less
        # than the minimum stepsize, set stepsize to min.
        elif (
            self.P - self.P_step / 10 > self.run.P_STOP
            and self.P_step / 10 <= self.run.DP_MIN
        ):
            self.P_step = self.run.DP_MIN
            self.P -= self.P_step

        else:
            while (
                self.P_step > self.run.DP_MIN
            ):  # Loop to get to the next possible pressure step
                if (
                    self.P - self.P_step / 10 > self.run.P_STOP
                ):  # If step/10 doesn't put you below min pressure
                    self.P_step = self.P_step / 10
                    self.P -= self.P_step
                    break
                else:
                    self.P_step = self.P_step / 10  # otherwise /10 and loop back around
            else:
                self.P -= (
                    self.P_step
                )  # This will trigger the outer loop when it sees self.P < P_STOP.

    def variable_step(self):
        """Reduces the step size by x10 when called during decompression"""

        if self.P_step / 10 >= self.run.DP_MIN:
            self.P = self.P + self.P_step - (self.P_step / 10)
            self.P_step = self.P_step / 10
            print(f"Convergence issue; reducing stepsize to {self.P_step}")
            return None
        else:
            return (
                "Minimum pressure step reached without achieving convergence "
                f"at {self.P} bar."
            )

    def rho_bulk(self, melt, gas):
        """Returns the density of the bulk magma (gas+melt combined)"""

        return (gas.rho() * self.GvF[-1]) + (
            melt.rho(P=self.P * 1e5) * (1 - self.GvF[-1])
        )


#  -------------------------------------------------------------------------
class Molecule:
    """
    Defines the properties of each molecular species in the gas phase

    Parameters
    ----------
    run : RunDef class
        The active instance of the RunDef class
    sys : ThermoSystem class
        The active instance of the ThermoSystem class
    melt : Melt class
        The active instance of the Melt class
    mol : str
        The name of the molecule being instantated

    Attributes
    ----------
    Mol : str
        Molecule name
    M : float
        Molecular mass (g/mol)
    PTcrit : dict
        Critical T's (K) and P's (bar) for calculating fugacity coefficients
        using the method from Shi and Saxena (1992)
    solCon : dict
        Contains the 'a' and 'b' solubility constants for use in solubility
        laws from Burguisser & Scaillet (2015).
    sys : ThermoSystem class
        Active instance of the ThermoSystem class
    melt : Melt class
        The active instance of the Melt class
    delG : float
        Gibbs free energy of formation
    Y : float
        Current fugacity coefficient
    y_constants : dict
        Stores any constants required for calculating the fugacity
        coefficient.
    """

    def __init__(self, run, sys, melt, mol):
        self.Mol = str(mol)  # Molecule name

        # initialise for subsequent settings
        self.M = cnst.m[self.Mol.lower()]  # Molecular mass

        # Critical T's (K) and P's (bar) respectively, for fugacity coefficients.
        self.PTcrit = {"T": cnst.PTcrit[self.Mol][0], "P": cnst.PTcrit[self.Mol][1]}

        if run.COMPOSITION == "basalt":
            if self.Mol == "S2":
                pass
            else:
                self.solCon = {
                    "a": cnst.basalt["a"][str(self.Mol)],
                    "b": cnst.basalt["b"][str(self.Mol)],
                }  # Solubility constant a and b

        elif run.COMPOSITION == "phonolite":
            lim = [1073.15, 1523.15]  # Tested temperature dependent solubility limits

            if self.Mol == "H2O" and lim[0] < sys.T < lim[1]:
                self.solCon = {
                    "a": cnst.phonolite["a"][str(self.Mol)](cnvs.K2C(sys.T)),
                    "b": cnst.phonolite["b"][str(self.Mol)](cnvs.K2C(sys.T)),
                }  # T-dependent solubility constant a and b

            elif self.Mol == "H2O" and (sys.T < lim[0] or sys.T > lim[1]):
                msgs.solubility_temp(sys.T, lim)
                self.solCon = {
                    "a": 4.69e-4,
                    "b": 0.6441,
                }  # T-independent solubility constant a and b
            else:
                self.solCon = {
                    "a": cnst.phonolite["a"][str(self.Mol)],
                    "b": cnst.phonolite["b"][str(self.Mol)],
                }  # Solubility constant a and b

        elif run.COMPOSITION == "rhyolite":
            lim = [1063.15, 1283.15]

            if (self.Mol == "H2O" or self.Mol == "CO2") and lim[0] < sys.T < lim[1]:
                self.solCon = {
                    "a": cnst.rhyolite["a"][str(self.Mol)](cnvs.K2C(sys.T)),
                    "b": cnst.rhyolite["b"][str(self.Mol)](cnvs.K2C(sys.T)),
                }  # Solubility constant a and b

            elif self.Mol == "H2O" and (sys.T < lim[0] or sys.T > lim[1]):
                msgs.solubility_temp(sys.T, lim)
                self.solCon = {
                    "a": 5.4344e-4,
                    "b": 0.6344,
                }  # T-independent solubility constant a and b

            elif self.Mol == "CO2" and (sys.T < lim[0] or sys.T > lim[1]):
                self.solCon = {
                    "a": 4.9358e-7,
                    "b": 1.04896,
                }  # T-independent solubility constant a and b

            else:
                self.solCon = {
                    "a": cnst.rhyolite["a"][str(self.Mol)],
                    "b": cnst.rhyolite["b"][str(self.Mol)],
                }  # Solubility constant a and b

        self.sys = sys
        self.melt = melt
        self.delG = 0  # Gibbs Free energy of formation
        self.Y = 0  # Fugacity coefficient
        self.y_constants = (
            {}
        )  # Stores any constants needed for the activity coefficient calculation

    def get_G(self, T):
        """
        Gets the Gibbs free energy of formation.

        Calculates the GfF for a given temperature, to use in calculations
        of equilibrium constants. Reads in critical data for the species
        from the 'data' folder, linearly interpolating between temperatures.

        Parameters
        ----------
        T : float
            Current temperature (K)

        Returns
        -------
        self.delG : float
            Gibbs free energy of formation
        """

        filename = f"{self.Mol}.txt"

        data_file = importlib.resources.files("evo") / "data" / filename

        with data_file.open("r", encoding="utf-8") as path:
            Temp_ref = []  # To store temperatures for interpolation
            del_G_ref = []  # To store values of G for interpolation
            inter_count = 0

            # Iterate through file to find correct values according to T provided.
            for aRow in path:
                values = aRow.split("\t")
                if not aRow.startswith("#"):  # Skips descriptive headers
                    if T - float(values[0]) == 0:
                        # Exact match found, no interpolation needed
                        self.delG = float(values[6])
                        return self.delG
                    elif inter_count == 1:
                        # Adds in the upper table values for interpolation
                        Temp_ref.append(float(values[0]))
                        del_G_ref.append(float(values[6]))
                        break
                    elif T - float(values[0]) <= 99:
                        Temp_ref.append(float(values[0]))
                        del_G_ref.append(float(values[6]))
                        inter_count += 1

        # Perform interpolation if needed
        if len(Temp_ref) == 2:
            self.delG = np.interp(T, Temp_ref, del_G_ref)

        return self.delG


#  -------------------------------------------------------------------------
class Melt:
    """
    Creates an object for the silicate melt, stores mass and chemical comp

    Parameters
    ----------
    run : RunDef class
        The active instance of the RunDef class
    sys : ThermoSystem class
        The active instance of the ThermoSystem class

    Attributes
    ----------
    graphite_sat : bool, default = False
        Records whether the melt is currently graphite saturated
    Cs : dict
        Stores the major oxide chemistry as real mass values (g).
        All iron as FeO.
    cm_dry : dict
        Dry major oxide composition as mole fractions, with feo+fe2o3
    cm_density : dict
        Stores the normalised major oxide comp plus H2O and CO2 as mole
        fractions, for use calculating the melt density.
    cm_volatiles : dict
        Stores the full melt composition, major oxides plus volatiles, as
        mole fractions.
    fmq : list of floats
        Stores the fO2 of the melt at every pressure step relative to the
        FMQ buffer
    h2o, h2, co2, co, ch4, sulfide, sulfate, n : list of floats
        Stores the melt weight fraction of each species after every
        pressure step
    s : list of floats
        Stores the total melt S weight fraction after every P step
    F : list of floats
        Stores the ferric/ferrous ratio of the melt after every P step
    ofe : list of floats
        Store the mass of O in iron oxides after every P step
    graph_current : float
        The number of moles in the system in the current timestep
    graphite : list of floats
        Stores the weight fraction of graphite in the melt after every P step
    rho_store : list of floats
        Stores the density of the volatile free melt after every P step
    """

    def __init__(self, run, sys):
        self.run = run
        self.sys = sys
        self.graphite_sat = False
        self.Cs = {}
        self.cm_dry = {}
        self.cm_density = {}
        self.cm_volatiles = {}
        self.fmq = []
        self.h2o = []
        self.h2 = []
        self.co2 = []
        self.co = []
        self.ch4 = []
        self.sulfide = []  # stores s2-
        self.sulfate = []  # stores s6+
        self.s = []
        self.n = []  # stores the melt n content, independent of speciation.
        self.F = []  # stores the mfe2o3/mFeO ratio at each step.
        self.ofe = []
        self.graph_current = (
            0  # number of moles of graphite in the system at current timestep
        )
        self.graphite = []  # stores graphite mass after each timesetep
        self.rho_store = []

    def chem_init(self, eles, chem):
        """
        Loads in major oxide composition from chemistry file.

        Gets the weight percent and actual masses of the oxide species
        in the silicate melt, adds in fe2o3 as a required species if
        missing, checks the ferric/ferrous model is appropriate for
        the amount of iron in the melt and calculates the fO2 if the
        ferric/ferrous ratio has been given.

        Parameters
        ----------
        eles : list of str
            The oxide names present in the silicate melt
        chem : list of floats
            The proportion of each species in the melt as weight PERCENT
        """

        # intialise constant masses of oxides vectors
        # as record of conserved system property
        self.C_s = {}  # actual oxide masses (g) NOT MASS FRACTION

        tmp_Cw = {}  # weight percent (the value which is read in by the file)

        length = 0  # counts through the values stored for the oxide composition
        sm = 0  # sums for average calculations

        # readin raw composition data
        for e in eles:
            tmp_Cw.update({e: chem[length]})
            sm += tmp_Cw[e]
            length += 1

        # normalise to 100 and define initial system element mass
        for e in eles:
            tmp_Cw[e] = tmp_Cw[e] * 100.0 / sm
            self.C_s[e] = tmp_Cw[e] / 100.0 * self.run.MASS  # actual oxide mass

        # add Fe+++ in if component has not been specified
        if "fe2o3" not in tmp_Cw:
            tmp_Cw["fe2o3"] = 0.0
            self.C_s["fe2o3"] = 0.0

        self.Cs = copy.copy(self.C_s)  # These are actual masses in g

        # This won't be completely accurate if fe2o3 is high
        # because they're not normalised; but close enough.
        if (
            self.run.FO2_MODEL == "r2013"
            and (tmp_Cw["feo"] + (tmp_Cw["fe2o3"] / cnst.m["fe2o3"]) * cnst.m["feo"])
            <= 15.0
        ):
            msgs.fo2_model_mismatch(
                (tmp_Cw["feo"] + (tmp_Cw["fe2o3"] / cnst.m["fe2o3"]) * cnst.m["feo"]),
                self.run,
            )
        elif (
            self.run.FO2_MODEL == "kc1991"
            and (tmp_Cw["feo"] + (tmp_Cw["fe2o3"] / cnst.m["fe2o3"]) * cnst.m["feo"])
            > 15.0
        ):
            msgs.fo2_model_mismatch(
                (tmp_Cw["feo"] + (tmp_Cw["fe2o3"] / cnst.m["fe2o3"]) * cnst.m["feo"]),
                self.run,
            )

        # if Fe3 and Fe2 exist, then need to determine system fO2;
        # also updates the Fe mass in system.
        if "feo" in eles and "fe2o3" in eles and self.sys.FO2 is None:
            self.sys.FO2 = np.log(
                cnvs.c2fo2(self.Cm(), self.sys.T, self.sys.Ppa, self.run.FO2_MODEL)
            )

        # store element names comprising liquid
        self.eles = self.Cs.keys()

    def iron_fraction(self, fo2, ppa=None):
        """
        Returns the melt composition and ferric/ferrous ratio for a given fO2

        Parameters
        ----------
        fo2 : float
            Ln(fO2)
        ppa : float, optional
            Current pressure in pascal, by default None. If None, pressure
            is taken from the ThermoSystem object

        Returns
        -------
        composition : dict
            major oxide composition as mole fractions,
              updated with new ferric/ferrous ratio
        F : float
            ferric/ferrous ratio
        """

        if ppa is None:
            ppa = self.sys.Ppa

        composition = self.Cm()  # molar composition just silicate melt

        # Outputs MOL fe2o3 / feo. Runs without dissolved volatiles to get the fraction,
        # then normalised with below.
        F = cnvs.fo2_2F(composition, self.sys.T, ppa, fo2, self.sys.run.FO2_MODEL)

        mFeoT = composition["feo"]

        mfe2o3 = F * mFeoT / (1 + 0.8998 * F)

        mfeo = mfe2o3 / F

        composition["feo"] = mfeo

        composition["fe2o3"] = mfe2o3

        composition = cnvs.norm(composition)
        return composition, F

    def Cw(self):
        "Converts ACTUAL MASS to WEIGHT PERCENT."
        return cnvs.norm(self.Cs, 100.0)

    def Cm(self):
        "Converts to ACTUAL MASS to MOLAR FRACTION"

        Cw = cnvs.norm(
            self.Cs, 1.0
        )  # converts real mass to mass fraction for the wt2mol conversion.
        return cnvs.wt2mol(Cw)

    def cm_singleO(self, mol_oxides):
        """
        Returns the melt composition as mol fracs on a single oxygen basis

        Parameters
        ----------
        mol_oxides : dict
            Melt composition as mole fractions

        Returns
        -------
        cation_formula : dict
            Melt composition as single oxygen mole fractions, with names
            as cations only.
        """

        def num_oxygen(x):
            """Returns the number of oxygen atoms in a species"""
            if x in ["na2o", "k2o", "mgo", "cao", "mno", "feo", "nio"]:
                return 1

            elif x in ["sio2", "tio2"]:
                return 2

            elif x in ["al2o3", "cr2o3", "fe2o3"]:
                return 3

            elif x == "p2o5":
                return 5

        def num_anion(x):
            """Returns the number of anion atoms in a species"""
            if x in ["sio2", "tio2", "mgo", "cao", "mno", "feo", "nio"]:
                return 1

            elif x in ["al2o3", "cr2o3", "fe2o3", "p2o5", "na2o", "k2o"]:
                return 2

        total_o = 0
        mol_cations = {}

        for name, x in mol_oxides.items():
            mol_cations[name] = x * num_anion(name)
            total_o += x * num_oxygen(name)

        formula = {}

        for name, x in mol_cations.items():
            formula[name] = x / total_o

        cation_formula = {}

        for x in formula.keys():
            if x == "feo":
                cation_formula["fe2"] = formula[x]
            elif x == "fe2o3":
                cation_formula["fe3"] = formula[x]
            else:
                no_num = re.sub(r"\d", "", x)
                new_name = re.sub("o", "", no_num)
                cation_formula[new_name] = formula[x]

        return cation_formula

    def formula_weight(self, fo2, P):
        """
        Calculates the single-O formula weight of the volatile-free melt.

        Finds the fe2/fe3 ratio first.

        Parameters
        ----------
        fO2 : float
            Absolute oxygen fugacity (bar)
        P : float
            Pressure (bar)

        Returns
        -------
        fwm : float
            The mean formula weight (g/mol)
        """

        mol_oxides = self.iron_fraction(np.log(fo2), ppa=P * 1e5)[0]

        formula = self.cm_singleO(mol_oxides)

        fwm = 0

        # strip down to cation name
        for x in formula.keys():
            no_num = re.sub(r"\d", "", x)
            new_name = re.sub("o", "", no_num)

            fwm += cnst.m[new_name] * formula[x]

        fwm += cnst.m["o"]

        return fwm

    def rho(self, T=None, P=None):
        """calculates the density of the volatile-free silicate melt"""
        if T is None:
            T = self.sys.T
        if P is None:
            P = self.sys.Ppa

        if self.cm_density:
            composition = self.cm_density
        else:
            composition = self.Cm()

        # K, Pa
        if self.run.DENSITY_MODEL == "spera2000":
            return density.den_calc_spera2000(composition, T, P)

    def sulfide_sat(self, comp):
        """
        Calculates the sulfur content at sulfide saturation of the melt

        Parameters
        ----------
        comp : dict
            The composition of the melt as mole fractions, including the
            dissolved volatile content.

        Returns
        -------
        float
            The sulfur content at sulfide saturation, in parts per million

        References
        ----------
        Liu, Y., Samaha, N-T., Baker, D.R. (2007) Sulfur concentration at
        sulfide saturation (SCSS) in magmatic silicate melts. GCA.
        """

        if self.run.SCSS == "liu2007":
            "Liu et al 2007"

            def mfm():
                Na = comp["na2o"] * 2
                K = comp["k2o"] * 2
                Ca = comp["cao"]
                Mg = comp["mgo"]
                Fe2 = comp["feo"]
                Si = comp["sio2"]
                Al = comp["al2o3"] * 2
                Fe3 = comp["fe2o3"] * 2

                return (Na + K + 2 * (Ca + Mg + Fe2)) / (Si * (Al + Fe3))

            MFM = mfm()

            mH2O = comp["h2o"]  # need to get mol frac of H2O IN MELT
            mFeO = comp["feo"]

            return np.exp(
                11.35251
                - 4454.6 / self.sys.T
                - 0.0319 * (self.sys.P / self.sys.T)
                + 0.71006 * np.log(MFM)
                - 1.98063 * (MFM * mH2O)
                + 0.21867 * np.log(mH2O)
                + 0.36192 * np.log(mFeO)
            )

    def melt_dicts(self):
        """
        Generates dicts containing the mole fraction composition of
        melt+volatiles, and one to calculate the density with.

        Returns
        -------
        melt_wet_mol : dict
            Melt + volatile composition as mole fractions,
             for melt.cm_volatiles
        density_mol : dict
            Melt + CO2 + H2O composition as mole fractions,
             for melt.cm_density
        """
        dry_wt = cnvs.mol2wt(self.cm_dry)
        dry_wt["wgt"] = self.sys.WgT[-1]

        dry_wt["h2o"] = self.h2o[-1]
        dry_wt["h2"] = self.h2[-1]
        dry_wt["co2"] = self.co2[-1]
        dry_wt["co"] = self.co[-1]
        dry_wt["ch4"] = self.ch4[-1]
        dry_wt["s"] = self.s[-1]
        dry_wt["n"] = self.n[-1]

        sys_wet_wt = cnvs.norm(
            dry_wt
        )  # PL: check this at a later date, seems superfluous?

        melt_wet_wt = dict(sys_wet_wt)
        melt_wet_wt.pop("wgt", None)
        melt_wet_wt = cnvs.norm(melt_wet_wt)
        melt_wet_mol = cnvs.wt2mol(
            melt_wet_wt
        )  # mol fraction comp of the melt, including all dissolved volatile species

        # Double check density isn't like sulphur where the melt is anhydrous,
        # plus wt water on top. PL: check co/ch4 don't need to go in?
        density = dict(sys_wet_wt)
        density.pop("wgt", None)
        density.pop("h2", None)
        density.pop("co", None)
        density.pop("ch4", None)
        density.pop("s", None)
        density.pop("n", None)
        density = cnvs.norm(density)
        density_mol = cnvs.wt2mol(
            density
        )  # includes all silicate elements plus dissolved H2O and CO2.

        return melt_wet_mol, density_mol

    def melt_composition(self, gas, mols):
        """
        Major function to run at end of P step, updates all results stores.

        Runs after a successful pressure step. Calculates and stores the
        new volatile content of the melt, saves the fO2, ferric/ferrous
        ratio and new major oxide fractions, and checks for both sulfide
        saturation and whether the sulfate fraction exceeds 10% (at which
        point the run ends).

        Parameters
        ----------
        gas : Gas class
            Active instance of the Gas class
        mols : list of Molecule classes
            A list of the active Molecule classes
        """

        # set the molecule names, and calculate the melt volatile weight percents.
        if self.run.GAS_SYS == "OH":
            H2O, O2, H2 = mols
        elif self.run.GAS_SYS == "COH":
            H2O, O2, H2, CO, CO2, CH4 = mols
        elif self.run.GAS_SYS == "SOH":
            H2O, O2, H2, S2, SO2, H2S = mols
        elif self.run.GAS_SYS == "COHS":
            H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S = mols
        elif self.run.GAS_SYS == "COHSN":
            H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2 = mols

        self.h2o.append(
            sl.h2o_melt(gas.mH2O[-1], H2O, self.sys.P, name=self.sys.run.H2O_MODEL)
            * cnst.m["h2o"]
        )
        self.h2.append(
            sl.h2_melt(gas.mH2[-1], H2, self.sys.P, self, name=self.sys.run.H2_MODEL)
            * cnst.m["h2"]
        )

        if (
            self.sys.run.GAS_SYS == "COHS"
            or self.sys.run.GAS_SYS == "COH"
            or self.sys.run.GAS_SYS == "COHSN"
        ):
            self.co2.append(
                sl.co2_melt(
                    (CO2.Y * gas.mCO2[-1] * self.sys.P),
                    CO2,
                    (O2.Y * gas.mO2[-1] * self.sys.P),
                    self.sys.T,
                    self.sys.P,
                    self,
                    name=self.sys.run.C_MODEL,
                )
                * cnst.m["co2"]
            )

            if self.sys.run.CO_MODEL != "None":
                self.co.append(
                    sl.co_melt(
                        (CO.Y * gas.mCO[-1] * self.sys.P),
                        self.sys.P,
                        name=self.sys.run.CO_MODEL,
                    )
                    * cnst.m["co"]
                )
            else:
                self.co.append(0.0)

            if self.sys.run.CH4_MODEL != "None":
                self.ch4.append(
                    sl.ch4_melt(
                        (CH4.Y * gas.mCH4[-1] * self.sys.P),
                        self.sys.P,
                        name=self.sys.run.CH4_MODEL,
                    )
                    * cnst.m["ch4"]
                )
            else:
                self.ch4.append(0.0)

            if self.graphite_sat is True:
                self.graphite.append(self.graph_current * cnst.m["c"])
            else:
                self.graphite.append(0.0)

        else:
            self.co2.append(0.0)
            self.co.append(0.0)
            self.ch4.append(0.0)
            self.graphite.append(0.0)

        if (
            self.sys.run.GAS_SYS == "COHS"
            or self.sys.run.GAS_SYS == "SOH"
            or self.sys.run.GAS_SYS == "COHSN"
        ):
            fo2 = O2.Y * gas.mO2[-1] * self.sys.P
            self.cm_dry, F = self.iron_fraction(np.log(fo2))

            if self.run.FE_SYSTEM is False:
                # otherwise this is done in fe_save
                self.F.append(F)

            self.sulfide.append(
                sl.sulfide_melt(
                    (S2.Y * gas.mS2[-1] * self.sys.P),
                    fo2,
                    self.sys.P,
                    self.sys.T,
                    self,
                    name=self.sys.run.SULFIDE_CAPACITY,
                )
                * cnst.m["s"]
            )
            self.sulfate.append(
                sl.sulfate_melt(
                    (S2.Y * gas.mS2[-1] * self.sys.P),
                    fo2,
                    self.sys.P,
                    self.sys.T,
                    self,
                    self.sys.run,
                    name=self.sys.run.SULFATE_CAPACITY,
                )
                * cnst.m["s"]
            )
            self.s.append(self.sulfide[-1] + self.sulfate[-1])

        else:
            if self.run.FE_SYSTEM is False:
                # otherwise this is done in self.fe_save
                fo2 = O2.Y * gas.mO2[-1] * self.sys.P
                self.cm_dry, F = self.iron_fraction(np.log(fo2))
                self.F.append(F)

            self.sulfide.append(0.0)
            self.sulfate.append(0.0)
            self.s.append(0.0)

        if self.sys.run.GAS_SYS == "COHSN":
            self.n.append(
                sl.n_melt(
                    gas.mN2[-1],
                    (O2.Y * gas.mO2[-1] * self.sys.P),
                    self.sys.P,
                    name=self.sys.run.N_MODEL,
                )
                * cnst.m["n"]
            )
        else:
            self.n.append(0.0)

        self.cm_volatiles, self.cm_density = self.melt_dicts()

        self.rho_store.append(self.rho())  # Saves the melt density at each step

        if self.sys.run.S_SAT_WARN is True and self.run.GAS_SYS in [
            "SOH",
            "COHS",
            "COHSN",
        ]:
            melt_mols = (
                self.cm_volatiles
            )  # The melt phase normalised to 1, including dissolved volatiles

            # Sulfide + sulfate content of the melt in ppm
            S = (
                sl.sulfide_melt(
                    (S2.Y * gas.mS2[-1] * self.sys.P),
                    (O2.Y * gas.mO2[-1] * self.sys.P),
                    self.sys.P,
                    self.sys.T,
                    self,
                    name=self.sys.run.SULFIDE_CAPACITY,
                )
                * cnst.m["s"]
                * 1e6
            )
            S6 = (
                sl.sulfate_melt(
                    (S2.Y * gas.mS2[-1] * self.sys.P),
                    (O2.Y * gas.mO2[-1] * self.sys.P),
                    self.sys.P,
                    self.sys.T,
                    self,
                    self.sys.run,
                    name=self.sys.run.SULFATE_CAPACITY,
                )
                * cnst.m["s"]
                * 1e6
            )

            SCSS = self.sulfide_sat(
                melt_mols
            )  # gives the sulphur content in ppm for sulphide saturation

            if S > SCSS:
                msgs.scss_warn(
                    SCSS, S, self.sys, gas, self
                )  # If sulfide content is above SCSS, writeout and exit.
            elif S6 / (S + S6) >= 0.1:
                msgs.sulfate_warn(
                    S6 / (S + S6), self.sys, gas, self
                )  # If the S6 content is >10%, writes out current data and closes
                # as so4 currently doesn't degass.


#  -------------------------------------------------------------------------


class Gas:
    """Stores the properties and results for the gas as it decompresses.
    Holds the results of each decompression step.

    Attributes
    ----------
    sys : ThermoSystem object
        Active instance of the ThermoSystem class
    fo2 : string of floats
        ln(fo2) after every pressure step
    mO2, mH2O, mH2, mCO, mCO2, mCH4, mSO2, mH2S, mS2, mN2 : list of floats
        Stores the mole fractions of each species in the gas phase after
        every pressure step
    atomicM : dict
        Stores the mass fraction of each element in the gas phase
    wt : dict
        Temporary store for the results of converting mole frac to wt frac
    Wt : dict
        Stores the gas weight fraction after every pressure step
    f : dict
        Stores the fugacity of each species in the gas phase after every
        pressure step
    M : list of floats
        Stores the mean molecular mass of the gas phase after every P step
    """

    def __init__(self, sys):
        self.sys = sys
        self.fo2 = []
        self.mO2 = []  # Stores mO2 value after each pressure step
        self.mH2O = []
        self.mH2 = []
        self.mCO = []
        self.mCO2 = []
        self.mCH4 = []
        self.mSO2 = []
        self.mH2S = []
        self.mS2 = []
        self.mN2 = []
        self.atomicM = {}
        self.wt = (
            {}
        )  # Temporary store for the results of converting mole frac to weight frac
        self.Wt = {}  # Store of weight fractions at each pressure step
        self.f = {}  # Store the fugacity of a species at each step.
        self.M = []  # Stores the mean molecular weight of the gas phase

    def get_WgT(self, melt, mols):
        """Calculates the gas weight fraction. NOT PERCENT."""
        if self.sys.run.GAS_SYS == "OH":
            H2O = mols
            wgt = (
                (self.sys.atomicM["o"] / cnst.m["o"])
                - sl.h2o_melt(
                    self.mH2O[-1], H2O, self.sys.P, name=self.sys.run.H2O_MODEL
                )
            ) / (
                (self.mH2O[-1] + 2 * self.mO2[-1])
                / (
                    self.mH2[-1] * cnst.m["h2"]
                    + self.mH2O[-1] * cnst.m["h2o"]
                    + self.mO2[-1] * cnst.m["o2"]
                )
            )

        elif self.sys.run.GAS_SYS == "COH":
            H2O, O2, H2, CO, CO2, CH4 = mols
            lst = {
                "o2": self.mO2[-1],
                "h2": self.mH2[-1],
                "h2o": self.mH2O[-1],
                "co": self.mCO[-1],
                "co2": self.mCO2[-1],
                "ch4": self.mCH4[-1],
            }

            mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

            wgt = (
                (
                    self.sys.atomicM["c"] / cnst.m["c"]
                    - sl.co2_melt(
                        (CO2.Y * self.mCO2[-1] * self.sys.P),
                        CO2,
                        (O2.Y * self.mO2[-1] * self.sys.P),
                        self.sys.T,
                        self.sys.P,
                        melt,
                        name=self.sys.run.C_MODEL,
                    )
                    - sl.co_melt(
                        (CO.Y * self.mCO[-1] * self.sys.P),
                        self.sys.P,
                        name=self.sys.run.CO_MODEL,
                    )
                    - sl.ch4_melt(
                        (CH4.Y * self.mCH4[-1] * self.sys.P),
                        self.sys.P,
                        name=self.sys.run.CH4_MODEL,
                    )
                    - melt.graph_current
                )
                / (self.mCO[-1] + self.mCO2[-1] + self.mCH4[-1])
            ) * sum(mjMj)

        elif self.sys.run.GAS_SYS == "SOH":
            H2O, O2, H2, S2, SO2, H2S = mols
            lst = {
                "o2": self.mO2[-1],
                "h2": self.mH2[-1],
                "h2o": self.mH2O[-1],
                "s2": self.mS2[-1],
                "so2": self.mSO2[-1],
                "h2s": self.mH2S[-1],
            }

            mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

            wgt = (
                (
                    self.sys.atomicM["h"] / (2 * cnst.m["h"])
                    - (
                        sl.h2o_melt(
                            self.mH2O[-1], H2O, self.sys.P, name=self.sys.run.H2O_MODEL
                        )
                        + sl.h2_melt(
                            self.mH2[-1],
                            H2,
                            self.sys.P,
                            melt,
                            name=self.sys.run.H2_MODEL,
                        )
                    )
                )
                / (self.mH2O[-1] + self.mH2[-1] + self.mH2S[-1])
            ) * sum(mjMj)

        elif self.sys.run.GAS_SYS == "COHS":
            H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S = mols
            lst = {
                "o2": self.mO2[-1],
                "h2": self.mH2[-1],
                "h2o": self.mH2O[-1],
                "co": self.mCO[-1],
                "co2": self.mCO2[-1],
                "ch4": self.mCH4[-1],
                "s2": self.mS2[-1],
                "so2": self.mSO2[-1],
                "h2s": self.mH2S[-1],
            }

            mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

            wgt = (
                (
                    self.sys.atomicM["c"] / cnst.m["c"]
                    - sl.co2_melt(
                        (CO2.Y * self.mCO2[-1] * self.sys.P),
                        CO2,
                        (O2.Y * self.mO2[-1] * self.sys.P),
                        self.sys.T,
                        self.sys.P,
                        melt,
                        name=self.sys.run.C_MODEL,
                    )
                    - sl.co_melt(
                        (CO.Y * self.mCO[-1] * self.sys.P),
                        self.sys.P,
                        name=self.sys.run.CO_MODEL,
                    )
                    - sl.ch4_melt(
                        (CH4.Y * self.mCH4[-1] * self.sys.P),
                        self.sys.P,
                        name=self.sys.run.CH4_MODEL,
                    )
                    - melt.graph_current
                )
                / (self.mCO[-1] + self.mCO2[-1] + self.mCH4[-1])
            ) * sum(mjMj)

        elif self.sys.run.GAS_SYS == "COHSN":
            H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2 = mols
            lst = {
                "o2": self.mO2[-1],
                "h2": self.mH2[-1],
                "h2o": self.mH2O[-1],
                "co": self.mCO[-1],
                "co2": self.mCO2[-1],
                "ch4": self.mCH4[-1],
                "s2": self.mS2[-1],
                "so2": self.mSO2[-1],
                "h2s": self.mH2S[-1],
                "n2": self.mN2[-1],
            }

            mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

            wgt = (
                (
                    self.sys.atomicM["c"] / cnst.m["c"]
                    - sl.co2_melt(
                        (CO2.Y * self.mCO2[-1] * self.sys.P),
                        CO2,
                        (O2.Y * self.mO2[-1] * self.sys.P),
                        self.sys.T,
                        self.sys.P,
                        melt,
                        name=self.sys.run.C_MODEL,
                    )
                    - sl.co_melt(
                        (CO.Y * self.mCO[-1] * self.sys.P),
                        self.sys.P,
                        name=self.sys.run.CO_MODEL,
                    )
                    - sl.ch4_melt(
                        (CH4.Y * self.mCH4[-1] * self.sys.P),
                        self.sys.P,
                        name=self.sys.run.CH4_MODEL,
                    )
                    - melt.graph_current
                )
                / (self.mCO[-1] + self.mCO2[-1] + self.mCH4[-1])
            ) * sum(mjMj)

        if wgt < 0.0:
            raise ValueError("Gas weight fraction is negative.")
        else:
            return wgt

    def get_vol_frac(self, melt):
        """
        Returns the gas volume fraction of the system.

        Eq 15 of Burguisser & Scaillet (2015).

        Parameters
        ----------
        melt : Melt class
            Active instance of the Melt class

        Returns
        -------
        GvF : float
            The volume fraction of gas in the system.
        """

        if self.sys.run.GAS_SYS == "OH":
            self.wt = cnvs.mols2wts(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1])
            empty_lists = {
                "CO2": 0.0,
                "CO": 0.0,
                "CH4": 0.0,
                "S2": 0.0,
                "SO2": 0.0,
                "H2S": 0.0,
                "N2": 0.0,
            }
            M = cnvs.mean_mol_wt(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1])
        elif self.sys.run.GAS_SYS == "COH":
            empty_lists = {"S2": 0.0, "SO2": 0.0, "H2S": 0.0, "N2": 0.0}
            self.wt = cnvs.mols2wts(
                H2O=self.mH2O[-1],
                O2=self.mO2[-1],
                H2=self.mH2[-1],
                CO=self.mCO[-1],
                CO2=self.mCO2[-1],
                CH4=self.mCH4[-1],
            )
            M = cnvs.mean_mol_wt(
                H2O=self.mH2O[-1],
                O2=self.mO2[-1],
                H2=self.mH2[-1],
                CO=self.mCO[-1],
                CO2=self.mCO2[-1],
                CH4=self.mCH4[-1],
            )
        elif self.sys.run.GAS_SYS == "SOH":
            empty_lists = {"CO2": 0.0, "CO": 0.0, "CH4": 0.0, "N2": 0.0}
            self.wt = cnvs.mols2wts(
                H2O=self.mH2O[-1],
                O2=self.mO2[-1],
                H2=self.mH2[-1],
                S2=self.mS2[-1],
                SO2=self.mSO2[-1],
                H2S=self.mH2S[-1],
            )
            M = cnvs.mean_mol_wt(
                H2O=self.mH2O[-1],
                O2=self.mO2[-1],
                H2=self.mH2[-1],
                S2=self.mS2[-1],
                SO2=self.mSO2[-1],
                H2S=self.mH2S[-1],
            )
        elif self.sys.run.GAS_SYS == "COHS":
            empty_lists = {"N2": 0.0}
            self.wt = cnvs.mols2wts(
                H2O=self.mH2O[-1],
                O2=self.mO2[-1],
                H2=self.mH2[-1],
                CO=self.mCO[-1],
                CO2=self.mCO2[-1],
                CH4=self.mCH4[-1],
                S2=self.mS2[-1],
                SO2=self.mSO2[-1],
                H2S=self.mH2S[-1],
            )
            M = cnvs.mean_mol_wt(
                H2O=self.mH2O[-1],
                O2=self.mO2[-1],
                H2=self.mH2[-1],
                CO=self.mCO[-1],
                CO2=self.mCO2[-1],
                CH4=self.mCH4[-1],
                S2=self.mS2[-1],
                SO2=self.mSO2[-1],
                H2S=self.mH2S[-1],
            )
        elif self.sys.run.GAS_SYS == "COHSN":
            empty_lists = None
            self.wt = cnvs.mols2wts(
                H2O=self.mH2O[-1],
                O2=self.mO2[-1],
                H2=self.mH2[-1],
                CO=self.mCO[-1],
                CO2=self.mCO2[-1],
                CH4=self.mCH4[-1],
                S2=self.mS2[-1],
                SO2=self.mSO2[-1],
                H2S=self.mH2S[-1],
                N2=self.mN2[-1],
            )
            M = cnvs.mean_mol_wt(
                H2O=self.mH2O[-1],
                O2=self.mO2[-1],
                H2=self.mH2[-1],
                CO=self.mCO[-1],
                CO2=self.mCO2[-1],
                CH4=self.mCH4[-1],
                S2=self.mS2[-1],
                SO2=self.mSO2[-1],
                H2S=self.mH2S[-1],
                N2=self.mN2[-1],
            )

        for key in self.wt:
            if key in self.Wt:
                self.Wt[key].append(self.wt[key])
            else:
                self.Wt[key] = [self.wt[key]]

        if empty_lists is not None:
            for key in empty_lists:
                if key in self.Wt:
                    self.Wt[key].append(empty_lists[key])
                else:
                    self.Wt[key] = [empty_lists[key]]

        GvF = 1 / (
            1
            + (
                (M * self.sys.P * 1e5 * (1 - self.sys.WgT[-1]))
                / (
                    cnst.R
                    * self.sys.T
                    * melt.rho(P=(self.sys.P * 1e5))
                    * self.sys.WgT[-1]
                )
            )
        )

        return GvF

    def get_fugacity(self, molecules, gas_phase):
        """Returns the fugacity of each species in the gas phase.

        Saves the fugacities directly to self.f. Inputs as [H2O, O2,...],
        [gas.mH2O[-1], gas.mO2[-1],...].

        Parameters
        ----------
        molecules : list of Molecule instances
            List of all the Molecule instances, corresponding to all
            species in the gas phase.
        gas_phase : list of floats
            List of gas phase mole fractions, in same order as `molecules`
        """

        for m, g in zip(molecules, gas_phase):
            if isinstance(m, str):
                if m in self.f:
                    self.f[m].append(0.0)
                else:
                    self.f[m] = [0.0]

            else:
                if m.Mol in self.f:
                    self.f[m.Mol].append(
                        m.Y * g * self.sys.P
                    )  # appends to the f dictionary for storage.
                else:
                    self.f[m.Mol] = [m.Y * g * self.sys.P]

    def rho(self):
        """Returns the density of the gas phase"""
        M = 0
        for mol in self.sys.SC:
            M += self.wt[mol] * cnst.m[mol.lower()]

        return ((M / 1000) * self.sys.P * 1e5) / (cnst.R * self.sys.T)

    def get_ys(self, system):
        """
        Returns the fugacity coefficients of all the gas phase species

        Calls the solvgas function `set_Y()`.

        Parameters
        ----------
        system : tuple of Molecule instances
            Tuple containing all the Molecule class instances,
            corresponding to all the species in the gas phase.
        """

        return sg.set_Y(self.sys, system)

    def get_atomic_mass(self, run, sys):
        """
        Calculates the mass fraction of each element in the gas phase only

        Parameters
        ----------
        run : RunDef class
            Active instance of the RunDef class
        sys : ThermoSystem class
            Active instance of the ThermoSystem class

        Warnings
        --------
        Currently unused, only valid for 'gas only' runs currently not
        implemented.
        """

        if run.GAS_SYS == "COHS" and sys.OCS is True:
            mO2 = self.mO2[-1]
            mH2O = self.mH2O[-1]
            mH2 = self.mH2[-1]
            mCO2 = self.mCO2[-1]
            mCO = self.mCO[-1]
            mCH4 = self.mCH4[-1]
            mSO2 = self.mSO2[-1]
            mS2 = self.mS2[-1]
            mH2S = self.mH2S[-1]

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
            self.atomicM["o"] = cnst.m["o"] * (
                (sys.WgT[-1] * mH2O) / sum(mjMj)
                + (2 * (sys.WgT[-1] * mO2) / sum(mjMj))
                + 2 * (sys.WgT[-1] * mSO2) / sum(mjMj)
                + (sys.WgT[-1] * mCO) / sum(mjMj)
                + 2 * (sys.WgT[-1] * mCO2) / sum(mjMj)
            )

            # WtH
            self.atomicM["h"] = (
                2
                * cnst.m["h"]
                * (
                    (sys.WgT[-1] * mH2O) / sum(mjMj)
                    + (sys.WgT[-1] * mH2) / sum(mjMj)
                    + (sys.WgT[-1] * mH2S) / sum(mjMj)
                    + 2 * ((sys.WgT[-1] * mCH4) / sum(mjMj))
                )
            )

            # WtC
            self.atomicM["c"] = cnst.m["c"] * (
                (sys.WgT[-1] * mCO) / sum(mjMj)
                + (sys.WgT[-1] * mCO2) / sum(mjMj)
                + (sys.WgT[-1] * mCH4) / sum(mjMj)
            )

            # WtS
            self.atomicM["s"] = cnst.m["s"] * (
                (2 * sys.WgT[-1] * mS2) / sum(mjMj)
                + (sys.WgT[-1] * mH2S) / sum(mjMj)
                + (sys.WgT[-1] * mSO2) / sum(mjMj)
            )

    def open_system(self, melt, fraction):
        """
        Reduces the gas fraction by removing an aliquot.

        Removes a gas fraction from the system in accordance with open
        system degassing. Fraction to be removed is controlled by the user
        in the environment file. If fraction = 0.9, 10% 0f the gas phase
        will remain after a pressure step.

        Then recalculates the masses of each volatile element in the
        system now some has been removed.

        Parameters
        ----------
        melt : Melt class
            Active instance of the Melt class
        fraction : float
            Fraction of the gas phase to be removed after every pressure step.
        """

        keep = 1 - fraction
        self.sys.WgT[-1] = self.sys.WgT[-1] * keep

        self.sys.atomicM["h"] = cnvs.atomicM_calc(
            self.sys, melt, self, "h", -1, WgT=self.sys.WgT[-1]
        )
        self.sys.atomicM["o"] = cnvs.atomicM_calc(
            self.sys, melt, self, "o", -1, WgT=self.sys.WgT[-1]
        )
        self.sys.atomicM["o_tot"] = cnvs.atomicM_calc(
            self.sys, melt, self, "o_tot", -1, WgT=self.sys.WgT[-1]
        )
        self.sys.atomicM["c"] = cnvs.atomicM_calc(
            self.sys, melt, self, "c", -1, WgT=self.sys.WgT[-1]
        )
        self.sys.atomicM["s"] = cnvs.atomicM_calc(
            self.sys, melt, self, "s", -1, WgT=self.sys.WgT[-1]
        )
        self.sys.atomicM["n"] = cnvs.atomicM_calc(
            self.sys, melt, self, "n", -1, WgT=self.sys.WgT[-1]
        )


#  -------------------------------------------------------------------------
class Output:
    """
    To request the desired outputs through the output control file

    Attributes
    ----------
    plot_melt_species : bool
        If True, the melt volatile content will be plotted against
        pressure at the end of the run.
    plot_gas_species_wt : bool
        If True, the composition of the gas phase as weight fractions
        will be plotted against pressure at the end of the run.
    plot_gas_species_mol : bool
        If True, the composition of the gas phase as mole fractions
        will be plotted against pressure at the end of the run.
    plot_gas_fraction : bool
        If True, the gas volume fraction will be plotted against pressure
    plot_fo2_dFMQ : bool
        If True, fO2 relative to FMQ will be plotted against pressure
    """

    n_par = 5

    def __init__(self):
        # Graphical outputs
        self.plot_melt_species = True
        self.plot_gas_species_wt = True
        self.plot_gas_species_mol = True
        self.plot_gas_fraction = True
        self.plot_fo2_dFMQ = False

    def set_outputs(self, params):
        """
        Sets class up with parameters from the file.

        Checks for invalid strings and values from the file.

        Parameters
        ----------
        params : dict
            dictionary of attribute:value pairs read in from yaml file.
        """

        for par, val in params.items():
            if self.zTest_Params(par):
                setattr(self, par, val)

            else:
                raise ValueError("Invalid output choice: %s" % par)

    # helper methods
    def zTest_Params(self, par):
        """Checks `par` is a valid attribute of the Output class"""

        for item in inspect.getmembers(self):
            if item[0] == par:
                return True
        return 0

    def print_conditions(self):
        """Prints the chosen outputs to the terminal after setup with file."""

        length = 0
        for item in inspect.getmembers(self):
            if not inspect.ismethod(self):
                if length < self.n_par:
                    print(item[0], ":", item[1])
                    length += 1
                else:
                    break
        print("\n")
