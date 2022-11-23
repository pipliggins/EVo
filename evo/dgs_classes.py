# dgs_classes

# ------------------------------------------------------------------------
# IMPORT
# ------------------------------------------------------------------------

import copy
import itertools
import inspect
import numpy as np
import re
import sys
from pathlib import Path

from evo import constants as cnst
from evo import conversions as cnvs
from evo import density
from evo import init_cons as ic
from evo import messages as msgs
from evo import solubility_laws as sl
from evo import solvgas as sg

# ------------------------------------------------------------------------
# CLASSES
# ------------------------------------------------------------------------

#  -------------------------------------------------------------------------
class RunDef:
    """ Creates an object defining the parameters of the degassing run.
    Contains system properties, start and stop T/P, FO2 and buffering.
    fo2 - Delta FMQ (Frost 1991, Rev. Min.)
    """

    # keep track of number of possible parameters (mainly for output purposes)
    n_par = int(40)
    
    # setup DEFAULT run parameters - will be overwritten by anything set in the environment file
    def __init__(self):

        #run type
        self.COMPOSITION = 'basalt'             # Magma composition
        self.RUN_TYPE = 'closed'                # Select the run type from closed system, open system or gas only.
        self.SINGLE_STEP = False                # Specifies single pressure or full decompression
        self.FIND_SATURATION = False            # Search for the saturation point based on melt composition; ignores P_Start and WgT
        self.GAS_SYS = "OH"                     # Sets the gas phase system user wants to look at
        self.FE_SYSTEM = False                  # Set to true to run iron equilibration steps
        self.OCS = False
        self.S_SAT_WARN = True                  # Turns the check on sulphur saturation on/off for a system including S and Fe_system on.

        # model parameters
        self.T_START = 1473.15                  # K
        self.P_START = 3000                     # bar
        self.P_STOP = 1.0                       # bar
        self.DP_MIN = 1.0                       # Minimum pressure step
        self.DP_MAX = 100.0                     # Maximum pressure step
        self.MASS = 100.0                       # System mass, g
        self.WgT = 0.001                        # Initial total gas FRACTION (not percent)
        self.LOSS_FRAC = 0.99                   # The fraction of the gas phase lost at each step during open system degassing

        # Select physical property and chemical solubility models
        self.DENSITY_MODEL = str('spera2000')
        self.FO2_MODEL = str('kc1991')
        self.FMQ_MODEL = str('frost1991')
        self.H2O_MODEL = str('burguisser2015')
        self.H2_MODEL = str('burguisser2015')
        self.C_MODEL = str('burguisser2015')
        self.CO_MODEL = str('None')
        self.CH4_MODEL = str('None')     
        self.SULFIDE_CAPACITY = str('oneill2020')
        self.SULFATE_CAPACITY = str('nash2019')
        self.SCSS = str('liu2007')
        self.N_MODEL = str('libourel2003')

        #initial gas fugacities
        self.FO2_buffer_SET = False
        self.FO2_buffer = 'FMQ'
        self.FO2_buffer_START = 0.0             #delta buffer value
        self.FO2_SET = False
        self.FO2_START = 0.0         # bar
        self.FH2_SET = False
        self.FH2_START = 0.24        # bar
        self.FH2O_SET = False
        self.FH2O_START = 1000.0       # bar
        self.FCO2_SET = False
        self.FCO2_START = 0.01

        # set to run using total atomic mass fractions - ppm
        self.ATOMIC_MASS_SET = False
        self.ATOMIC_H = 550.0
        self.ATOMIC_C = 200.0
        self.ATOMIC_S = 4000.0
        self.ATOMIC_N = 10.0

        # initial wt frac volatile in magma (at given T and P, not total)
        self.WTH2O_SET = False
        self.WTH2O_START = 0.03                 # initial melt FRACTION
        self.WTCO2_SET = False
        self.WTCO2_START = 0.01                 # initial melt FRACTION
        self.SULFUR_SET = False
        self.SULFUR_START = 0.001               # initial melt FRACTION
        self.NITROGEN_SET = False
        self.NITROGEN_START = 0.01              # initial melt FRACTION
        self.GRAPHITE_SATURATED = False
        self.GRAPHITE_START = 0.0               # initial melt graphite content


    def param_set(self, params):

        # expect dictionary of key:value pairs
        # force type based on what's already been set
        for par, val in params.items():
            if self.zTest_Params(par):
                setattr(self, par, val)

            else:
                sys.exit("Warning: %s not a valid calculation parameter" % par)

        # Check for valid parameters
        
        if self.OCS == True:
            exit("Error: OCS functionality has not been released yet. Sorry!")
        
        if self.COMPOSITION.lower() != 'basalt' and self.COMPOSITION.lower() != 'phonolite' and self.COMPOSITION.lower() != 'rhyolite':
            exit(f"Error: The composition '{self.COMPOSITION}' is invalid. Please enter 'basalt', 'phonolite', or 'rhyolite'.")
        else:
            self.COMPOSITION = self.COMPOSITION.lower()

        if self.RUN_TYPE.lower() == 'closed' or self.RUN_TYPE.lower() == 'open':
            self.RUN_TYPE == self.RUN_TYPE.lower()
        else:
            exit(f"RUN_TYPE '{self.RUN_TYPE}'' is not recognised. Please select either 'closed' or 'open'; gas-only is not available in this release.")
        
        if self.MASS <= 0:
            sys.exit("Error: The system has no mass value")
        
        if self.WgT <= 0:
            sys.exit("Error: A gas weight fraction of greater than 0 must be entered")


    # helper methods
    def zTest_Params(self, par):  # z is a fudge to put at end of print list, so l stops the 'method class' stuff printing

        for item in inspect.getmembers(self):
            if item[0] == par:
                return True
        return 0

    # Prints the run conditions to the console
    def print_conditions(self):
        
        l = 0
        for item in inspect.getmembers(self):
            if l < self.n_par:
                print(item[0], ':', item[1])
                l += 1
            else:
                break


#  -------------------------------------------------------------------------
class ThermoSystem:
    """ Thermodynamic properties of the whole system (gas and melt) """

    def __init__(self, run):
        self.run = run
        self.T = self.run.T_START       # K
        self.P = self.run.P_START       # bar
        self.P_step = self.run.DP_MAX   # bar
        self.P_track = []               # stores the pressure steps run at
        self.M = self.run.MASS          # g
        self.rho = []                   # bulk density
        self.sat_conditions = 0         # store saturation pressure + gas speciation
        self.graph_unsat_rerun = False  # If melt goes graphite undersaturated, track that this is a rerun to reset to graph_saturated if the stepsize needs to reduce

        # intialise for subsequent setting
        self.FO2 = None                 # fo2 # STORED AS A LN(FO2) SO DO e**SELF.FO2 WHERE IT'S NEEDED
        self.FO2_buffer = None          # stores the fo2 as a value relative to the chosen buffer, default is FMQ.
        self.FH2 = None                 # fH2
        self.SC = []                    # System chemistry - i.e. gas species present.
        self.K = {}                     # Stores K values (equilibrium constants)
        self.atomicM = {}               # Atomic masses for mass balance - weight fraction
        self.WgT = [self.run.WgT]       # Total gas wt%, initially set by user.
        self.GvF = []                   # Stores the gas volume fraction at each step
        self.Cw = {}                    # SYSTEM mass fractions of oxides and volatile elements (C, Ogas, OFe, H + S), normalised to 1. iron as fe rather than feo + fe2o3, so don't count the oxygen twice.
        self.OCS = False                # Can be flipped to True after the last pressure step for a final run.

    @property
    def Ppa(self):
        """
        Converts system pressure in bar, to pressure in pascal (Pa)
        """
        return self.P*1e5
    
    # Selects the system chemistry from the system the user wants to examine.
    def gas_system(self, run):

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

        if run.GAS_SYS == 'OH':
            ic.oh(self, run, melt, gas, mols)
        elif run.GAS_SYS == 'COH':
            ic.coh(self, run, melt, gas, mols)
        elif run.GAS_SYS =='SOH':
            ic.soh(self, run, melt, gas, mols)
        elif run.GAS_SYS == 'COHS':
            ic.cohs(self, run, melt, gas, mols)
        elif run.GAS_SYS == 'COHSN':
            ic.cohsn(self, run, melt, gas, mols)

    def norm_with_gas(self, melt, mols, gas):  # for iron equilibration.
        cons = 0  # holds the values that DON'T want to be normalised, i.e. volatile species and iron content
        C_s = {}  # Holds the values to be normalised
        C_w = melt.Cw()  # wt %

        # Adds the mass of volatile species in the system to cons as masses of each atomic element
        for key in self.atomicM:
            self.Cw[key] = self.atomicM[key]
            cons += self.atomicM[key]  # the total wt fraction of volatiles both in gas phase and dissolved.

        # The oxide composition of the melt partitioned as wt %
        for ele in C_w:
            if ele == 'feo' or ele == 'fe2o3': # PL: This doesn't convert fe2o3 to feo at the moment, it probably should? Need to catch scenario where only fe2o3 is provided?
                cons += C_w[ele]/100
                self.Cw[ele] = C_w[ele]/100
            else:
                # oxide weights only, without volatiles or iron in any form, ready to be scaled.
                C_s[ele] = C_w[ele]/100

        # returns the 9 major oxides as weight %, with feoT and volatiles held constant.
        oxides = cnvs.norm(C_s, (1 - cons))

        self.Cw.update(oxides)  # adds the new normalised oxide values to the total system mass frac dictionary

        melt.cm_dry, F = melt.iron_fraction(self.FO2) # this is now a duplicate of the one inside init cons, hopefully.
        melt.F.append(F)

    def get_wto_tot(self, melt):
        # This is after the norm with gas step, so factors in the gas equilibration.
        # update self.Cw so it has feo/fe2o3 ratio in, and then change the volatile storage so stores as o, ogas and ofe

        Cm = cnvs.wt2mol(self.Cw)
        
        for k, v in melt.cm_dry.items():
            Cm[k] = v
        
        Cm = cnvs.norm(Cm, cons=['o', 'h', 'c', 's', 'n'])
        
        F = melt.F[-1]

        Cm['fe'] = Cm['feo'] + 2*Cm['fe2o3']
        m_ofe = Cm['feo'] + 3*Cm['fe2o3']  # The mole frac of oxygen tied up in the melt.

        Cm.pop('feo', None)
        Cm.pop('fe2o3', None)    
        
        Cm['ogas'] = Cm['o']
        Cm['ofe'] = m_ofe
        Cm.pop('o', None)

        Cm = cnvs.norm(Cm, cons=['fe', 'ogas', 'ofe', 'h', 'c', 's', 'n'])
        self.Cw = dict(cnvs.mol2wt(Cm))
        
        self.atomicM['fe'] = self.Cw['fe']
        self.atomicM['ofe'] = self.Cw['ofe']
        self.atomicM['o_tot'] = self.atomicM['o'] + self.Cw['ofe']

        return self.atomicM
    
    def fe_equil(self, melt, mo2, O2):
        """
        Calculate the amount of oxygen in the melt for use in solver.

        Args:
            melt (class): active instance of the Melt class
            mo2 (mpfr): current mol fraction of O2 in the gas phase
            O2 (class): O2 instance of the Molecule class

        Returns:
            ofe (mpfr): atomic mass of O in the melt.
        """

        fo2 = np.log(O2.Y * mo2 * self.P)

        # returns anhydrous silicate melt comp mol fraction & F, the ratio of Fe2:Fe3. Only need F here.
        F = melt.iron_fraction(fo2)[1]  # Once iron equilibration is done, F needs to be stored.

        ofe = cnst.m['o'] * (self.atomicM['fe']/cnst.m['fe']) * ((1 +3*F)/(1 + 2*F))       # THIS STEP total weight fraction of atomic O in iron

        return ofe

    def fe_save(self, melt, mo2, O2):
        """
        Saveout the amount of O in the gas and melt respectively after a successful fe_equil run.
        """

        fo2 = np.log(O2.Y * mo2 * self.P)

        # returns anhydrous silicate melt comp mol fraction & F, the ratio of Fe2:Fe3. Only need F here.
        F = melt.iron_fraction(fo2)[1]  # Once iron equilibration is done, F needs to be stored.

        ofe = cnst.m['o'] * (self.atomicM['fe']/cnst.m['fe']) * ((1 +3*F)/(1 + 2*F))       # THIS STEP total weight fraction of atomic O in iron

        self.atomicM['ofe'] = ofe
        self.atomicM['o'] = self.atomicM['o_tot'] - ofe
        melt.ofe.append(ofe)
        melt.F.append(F)
    
    def mass_conservation_reset(self, melt, gas):
        """
        If mass is no longer conserved go back and rerun last step with lower stepsize.

        If the atomic masses are no longer being conserved, delete the results of the 
        last run and reduces the stepsize to help convergence towards the correct
        solution.

        Args:
            melt (class): active instance of the Melt class
            gas (class): active instance of the Gas class
        
        Returns:
            None
        """

        sys_lsts = [self.WgT, self.GvF, self.rho]
        gas_lsts = [gas.fo2, gas.mH2O, gas.mO2, gas.mH2, gas.mCO, gas.mCO2, gas.mCH4, gas.mS2, gas.mSO2, gas.mH2S, gas.mN2, gas.M]
        melt_lsts = [melt.fmq, melt.h2o, melt.h2, melt.co2, melt.co, melt.ch4, melt.graphite, melt.sulfide, melt.sulfate, melt.s, melt.n, melt.F, melt.ofe]
        gas_wts_f = ['H2O', 'O2', 'H2', 'CO', 'CO2', 'CH4', 'S2', 'SO2', 'H2S', 'N2']
        
        for x in sys_lsts + gas_lsts + melt_lsts:
            del x[-1]
        
        for x in gas_wts_f:
            del gas.Wt[x][-1]
            if x != 'O2':
                del gas.f[x][-1]

        if melt.graphite_sat == True:
            melt.graph_current = melt.graphite[-1]
        
        if self.P_step / 10 >= self.run.DP_MIN:
            del self.P_track[-1]
            self.P = self.P + self.P_step
            self.P_step = self.P_step / 10
            print('mass conservation reset')
        else:
            del self.P_track[-1]
            msgs.closed_earlyexit(self, gas, melt)
            exit('Error: Mass is no longer being conserved. Please reduce the minimum pressure stepsize and try again.\nExiting...')

    # Reduces the size of the pressure step by factors of 10 as it approaches the P_STOP value.
    def pressure_step(self):

        if self.P - self.P_step > self.run.P_STOP: # If taking another pressure step won't take you below the final pressure specified, continue.
            self.P -= self.P_step
        
        elif self.P - self.P_step/10 > self.run.P_STOP and self.P_step/10 <= self.run.DP_MIN: # If taking another step puts you below the min pressure and step/10 is less than the minimum stepsize, set stepsize to min.
            self.P_step = self.run.DP_MIN
            self.P -= self.P_step
        
        else:
            while self.P_step > self.run.DP_MIN:  # Loop to get to the next possible pressure step
                if self.P - self.P_step/10 > self.run.P_STOP:   # If step/10 doesn't put you below min pressure
                    self.P_step = self.P_step / 10
                    self.P -= self.P_step
                    break
                else:
                    self.P_step = self.P_step / 10  # otherwise /10 and loop back around
            else:
                self.P -= self.P_step  # This will trigger the outer loop when it sees self.P < P_STOP.
        
    def variable_step(self):
        # reduces the step size by x10 when called during decompression

        if self.P_step / 10 >= self.run.DP_MIN:
            self.P = self.P + self.P_step - (self.P_step/10)
            self.P_step = self.P_step / 10
            print(f'Convergence issue; reducing stepsize to {self.P_step}')
            return None
        else:
            return f'Minimum pressure step reached without achieving convergence at {self.P} bar.'

    def rho_bulk(self, melt, gas):
        return (gas.rho()*self.GvF[-1]) + (melt.rho(P=self.P*1e5)*(1-self.GvF[-1]))


#  -------------------------------------------------------------------------
class Molecule:
    """Defines the properties of each molecular species in the gas phase"""

    def __init__(self, run, sys, melt, mol):

        self.Mol = str(mol)  # Molecule name

        # initialise for subsequent settings
        self.M = cnst.m[self.Mol.lower()]  # Molecular mass

        # Critical T's (K) and P's (bar) respectively, for fugacity coefficient calculation. Shi and Saxena 1992
        self.PTcrit = {"T": cnst.PTcrit[self.Mol][0], "P": cnst.PTcrit[self.Mol][1]}

        if run.COMPOSITION == 'basalt':
            if self.Mol == 'S2':
                pass
            else:
                self.solCon = {"a": cnst.basalt['a'][str(self.Mol)], "b": cnst.basalt['b'][str(self.Mol)]} # Solubility constant a and b
        
        elif run.COMPOSITION == 'phonolite':
            lim = [1073.15, 1523.15] # Tested temperature dependent solubility limits 

            if self.Mol == 'H2O' and lim[0] < sys.T < lim[1]:
                self.solCon = {"a": cnst.phonolite['a'][str(self.Mol)](cnvs.K2C(sys.T)), "b": cnst.phonolite['b'][str(self.Mol)](cnvs.K2C(sys.T))} # T-dependent solubility constant a and b
            
            elif self.Mol == 'H2O' and (sys.T < lim[0] or sys.T > lim[1]):
                msgs.solubility_temp(sys.T, lim)
                self.solCon = {"a": 4.69e-4, "b": 0.6441} # T-independent solubility constant a and b
            else:
                self.solCon = {"a": cnst.phonolite['a'][str(self.Mol)], "b": cnst.phonolite['b'][str(self.Mol)]} # Solubility constant a and b
        
        elif run.COMPOSITION == 'rhyolite':
            lim = [1063.15, 1283.15]
            
            if (self.Mol == 'H2O' or self.Mol == 'CO2') and lim[0] < sys.T < lim[1]:               
                self.solCon = {"a": cnst.rhyolite['a'][str(self.Mol)](cnvs.K2C(sys.T)), "b": cnst.rhyolite['b'][str(self.Mol)](cnvs.K2C(sys.T))} # Solubility constant a and b
            
            elif self.Mol == 'H2O' and(sys.T < lim[0] or sys.T > lim[1]):
                msgs.solubility_temp(sys.T, lim)
                self.solCon = {"a": 5.4344e-4, "b": 0.6344} # T-independent solubility constant a and b
            
            elif self.Mol == 'CO2' and(sys.T < lim[0] or sys.T > lim[1]):
                # No message as only needs to happen once (pick H2O as in all runs unlike CO2); if out of T bounds and continuing, assume with fixed constants.
                self.solCon = {"a": 4.9358e-7, "b": 1.04896} # T-independent solubility constant a and b
            
            else:
                self.solCon = {"a": cnst.rhyolite['a'][str(self.Mol)], "b": cnst.rhyolite['b'][str(self.Mol)]} # Solubility constant a and b
        
        self.sys = sys
        self.melt = melt
        self.delG = 0           # Gibbs Free energy of formation
        self.Y = 0              # Fugacity coefficient
        self.y_constants = {}   # Stores any constants needed for the activity coefficient calculation

    def get_G(self, sys):
        """Gets the Gibbs free energy of formation for a given temperature for K calculation.
        Assumes that T is going to be constant - if not then read in ALL values of G for later pick + choose.
        """
        T = sys.T

        file_path = Path(__file__).parent / f"Data/{self.Mol}.txt"
        
        path = open(file_path, "r")                      # Opens the file for data for the molecule
        Temp_ref = []                                    # To store temperatures which need to be interpolated
        del_G_ref = []                                   # To store values of G which need to be interpolated
        inter_count = 0
        for aRow in path:                                # iterates through the file to find the correct values according to the temperature provided.
            values = aRow.split('\t')                    # indicates tab delimited
            if not aRow.startswith('#'):                 # takes out descriptive headers
                if T - float(values[0]) == 0:            # can find K with values in table, no interpolation needed
                    self.delG = float(values[6])         # adds name of mol. and delG. of form. to the del_G dictionary
                    break
                elif inter_count == 1:                   # adds in the upper table values for interpolation
                    Temp_ref.append(float(values[0]))
                    del_G_ref.append(float(values[6]))
                    break
                elif T - float(values[0]) <= 99:
                    Temp_ref.append(float(values[0]))
                    del_G_ref.append(float(values[6]))
                    inter_count += 1
        if len(Temp_ref) == 2:
            self.delG = np.interp(T, Temp_ref, del_G_ref)  # interpolates between 2 and returns value for Gf.
        path.close()
        return self.delG

#  -------------------------------------------------------------------------
class Melt:
    """Creates an object which is specifically the melt - stores mass and chemical composition"""

    def __init__(self, run, sys):
        self.run = run
        self.sys = sys
        self.graphite_sat = False
        self.Cs = {}  # stores magma chemistry as a real mass value (g). Can be converted to mol fraction using cm(), or weight percent using cw(). All iron as FeO.
        self.cm_dry = {}  # Stores silicate melt composition as a mol fraction including both feo and fe2o3
        self.cm_density = {}  # Stores the mol fraction composition of the silicate melt plus dissolved H2O and CO2 content for melt density calculation.
        self.cm_volatiles = {}  # Stores the mol fractions of the entire melt here, all volatiles plus silicate ions, normalised to 1
        self.fe_mass = 0
        self.fmq = []
        self.h2o = []
        self.h2 = []
        self.co2 = []
        self.co = []
        self.ch4 = []
        self.so2 = []
        self.h2s = []
        self.sulfide = []  # stores s2-
        self.sulfate = []  # stores s6+
        self.s = []  # stores the melt S content (mass of S in so2 and h2s dissolved) - now based on s2(g) -> s2-
        self.n = [] # stores the melt n content, independent of speciation.
        self.F = []  # stores the mfe2o3/mFeO ratio at each step.
        self.ofe = []
        self.graph_current = 0      # number of moles of graphite in the system at current timestep
        self.graphite = []          # stores graphite mass after each timesetep
        self.rho_store = []

    def chem_init(self, eles, chem):

        # intialise constant masses of oxides vectors
        # as record of conserved system property
        self.C_s = {}           # actual oxide masses (g) NOT MASS FRACTION

        tmp_Cw = {}            # weight percent (the value which is read in by the file)

        l = 0                   # counts through the values stored for the oxide composition
        sm = 0                  # sums for average calculations

        # readin raw composition data
        for e in eles:
            tmp_Cw.update({e: chem[l]})
            sm += tmp_Cw[e]
            l += 1

        # normalise to 100 and define initial system element mass
        for e in eles:
            tmp_Cw[e] = tmp_Cw[e] * 100.0 / sm
            self.C_s[e] = tmp_Cw[e] / 100.0 * self.run.MASS  # converts % to decimal and multiplies by the system mass to get actual oxide mass

        # add Fe+++ in if component has not been specified
        if 'fe2o3' not in tmp_Cw:
            tmp_Cw['fe2o3'] = 0.0
            self.C_s['fe2o3'] = 0.0

        self.Cs = copy.copy(self.C_s)  # These are actual masses in g

        # This won't be completely accurate is fe2o3 is high because they're not normalised; but close enough.
        if self.run.FO2_MODEL == "r2013" and (tmp_Cw['feo'] + (tmp_Cw['fe2o3']/cnst.m['fe2o3'])*cnst.m['feo']) <= 15.0:
            msgs.fo2_model_mismatch((tmp_Cw['feo'] + (tmp_Cw['fe2o3']/cnst.m['fe2o3'])*cnst.m['feo']), self.run)
        elif self.run.FO2_MODEL == "kc1991" and (tmp_Cw['feo'] + (tmp_Cw['fe2o3']/cnst.m['fe2o3'])*cnst.m['feo']) > 15.0:
            msgs.fo2_model_mismatch((tmp_Cw['feo'] + (tmp_Cw['fe2o3']/cnst.m['fe2o3'])*cnst.m['feo']), self.run)

        # if Fe+++ and Fe++ exist, then need to determine system fO2; also updates the Fe mass in system.
        if 'feo' in eles and 'fe2o3' in eles and self.sys.FO2 is None:
            self.sys.FO2 = np.log(cnvs.c2fo2(self.Cm(), self.sys.T, self.sys.Ppa, self.run.FO2_MODEL))

        # store element names comprising liquid
        self.eles = self.Cs.keys()

    def iron_fraction(self, fo2, ppa = None):

        if ppa == None:
            ppa = self.sys.Ppa
        
        composition = self.Cm() # molar composition just silicate melt; this matches other calculations of the fe2/fe3 ratio i've seen but not D-compress.

        F = cnvs.fo2_2F(composition, self.sys.T, ppa, fo2, self.sys.run.FO2_MODEL)  # Outputs MOL fe2o3 / feo. Runs without dissolved volatiles to get the fraction, then normalised with below.
        
        mFeoT = composition['feo']

        mfe2o3 = F*mFeoT/(1+0.8998*F)

        mfeo = mfe2o3/F

        composition['feo'] = mfeo
        
        composition['fe2o3'] = mfe2o3

        composition = cnvs.norm(composition)
        return composition, F

    def Cw(self):
        "Converts ACTUAL MASS to WEIGHT PERCENT."
        return cnvs.norm(self.Cs, 100.)

    def Cm(self):
        "Converts to ACTUAL MASS to MOLAR FRACTION"

        Cw = cnvs.norm(self.Cs, 1.)  # converts real mass to mass fraction for the wt2mol conversion.
        return cnvs.wt2mol(Cw)

    def cm_singleO(self, mol_oxides):
        """
        Converts major element oxides as wt frac, to the composition as mol fractions
        on a single oxygen basis.

        Args:
            mol_oxides (dict): Melt composition as mol fraction

        Returns:
            cation_formula (dict): Melt composition as single oxygen mole fractions, with names as cations only.
        """

        def num_oxygen(x):
            
            if x in ['na2o', 'k2o', 'mgo', 'cao', 'mno', 'feo', 'nio']:
                return 1
            
            elif x in ['sio2', 'tio2']:
                return 2

            elif x in ['al2o3', 'cr2o3', 'fe2o3']:
                return 3
                    
            elif x == 'p2o5':
                return 5
        
        def num_anion(x):
            if x in ['sio2', 'tio2', 'mgo', 'cao', 'mno', 'feo', 'nio']:
                return 1

            elif x in ['al2o3', 'cr2o3', 'fe2o3', 'p2o5', 'na2o', 'k2o']:
                return 2
        
        total_o = 0
        mol_cations = {}
        
        for name, x in mol_oxides.items():
            mol_cations[name] = x * num_anion(name)
            total_o += x * num_oxygen(name)
            
        formula = {}
        
        for name, x in mol_cations.items():
            formula[name] = x/total_o
        
        cation_formula = {}

        for x in formula.keys():
            if x == 'feo':
                cation_formula['fe2'] = formula[x]
            elif x == 'fe2o3':
                cation_formula['fe3'] = formula[x]
            else:
                no_num = re.sub('\d', '', x)
                new_name = re.sub('o', '', no_num)
                cation_formula[new_name] = formula[x]
        
        return cation_formula

    def formula_weight(self, fo2, P):
        """
        Calculates the single-O formula weight of the volatile-free melt.
        Finds the fe2/fe3 content first.

        Args:
            fO2 (float): Oxygen fugacity (bar)
            P (float): Pressure (bar)
        
        Returns:
            fwm (float): The mean formula weight (g/mol)
        """

        mol_oxides = self.iron_fraction(np.log(fo2), ppa = P*1e5)[0]

        formula = self.cm_singleO(mol_oxides)

        fwm = 0

        # strip down to cation name
        for x in formula.keys():
            no_num = re.sub('\d', '', x)
            new_name = re.sub('o', '', no_num)
            
            fwm += cnst.m[new_name] * formula[x]
        
        fwm += cnst.m['o']

        return fwm
    
    def rho(self,T=None,P=None):
        """calculates the density of the silicate melt phase (without any dissolved volatiles)"""
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
        "Comp is the dictionary containing the normalised mol fraction of the melt phase, including dissolved volatiles."
        
        if self.run.SCSS == 'liu2007':
            "Liu et al 2007"
            
            def mfm():
                
                Na = comp['na2o']*2
                K = comp['k2o']*2
                Ca = comp['cao']
                Mg = comp['mgo']
                Fe2 = comp['feo']
                Si = comp['sio2']
                Al = comp['al2o3']*2
                Fe3 = comp['fe2o3']*2

                return (Na + K + 2*(Ca+Mg+Fe2))/(Si*(Al+Fe3))
            
            MFM = mfm()       

            mH2O = comp['h2o']  # need to get mol frac of H2O IN MELT
            mFeO = comp['feo']
            
            return np.exp(11.35251 - 4454.6/self.sys.T - 0.0319*(self.sys.P/self.sys.T) + 0.71006*np.log(MFM) - 1.98063*(MFM*mH2O) + 0.21867*np.log(mH2O) + 0.36192*np.log(mFeO))

    def melt_dicts(self):
        dry_wt = cnvs.mol2wt(self.cm_dry)
        dry_wt['wgt'] = self.sys.WgT[-1]

        dry_wt['h2o'] = self.h2o[-1]
        dry_wt['h2'] = self.h2[-1]
        dry_wt['co2'] = self.co2[-1]
        dry_wt['co'] = self.co[-1]
        dry_wt['ch4'] = self.ch4[-1]
        dry_wt['s'] = self.s[-1]
        dry_wt['n'] = self.n[-1]
        
        sys_wet_wt = cnvs.norm(dry_wt)

        melt_wet_wt = dict(sys_wet_wt)
        melt_wet_wt.pop('wgt', None)
        melt_wet_wt = cnvs.norm(melt_wet_wt)
        melt_wet_mol = cnvs.wt2mol(melt_wet_wt) # Dictionary containing the mol fraction comp of the melt, including all dissolved volatile species

        # Double check density isn't like sulphur where the melt is anhydrous, plus wt water on top. PL: check co/ch4 don't need to go in?
        density = dict(sys_wet_wt)
        density.pop('wgt', None)
        density.pop('h2', None)
        density.pop('co', None)
        density.pop('ch4', None)
        density.pop('s', None)
        density.pop('n', None)
        density = cnvs.norm(density)  
        density_mol = cnvs.wt2mol(density)  # Now a mol frac dict including all silicate elements plus dissolved H2O and CO2.

        return melt_wet_mol, density_mol
    
    def melt_composition(self, gas, mols):
        """Calculates the amount of volatile in the melt, creates a normalised melt composition of minerals and dissolved volatiles,
        then checks for sulphide (and eventually carbon) saturation."""        
        
        # set the molecule names, and calculate the melt volatile weight percents.
        if self.run.GAS_SYS == 'OH':
            H2O, O2, H2 = mols
        elif self.run.GAS_SYS == 'COH':
            H2O, O2, H2, CO, CO2, CH4 = mols
        elif self.run.GAS_SYS == 'SOH':
            H2O, O2, H2, S2, SO2, H2S = mols
        elif self.run.GAS_SYS == "COHS":
            H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S = mols
        elif self.run.GAS_SYS == "COHSN":
            H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2 = mols

        self.h2o.append(sl.h2o_melt(gas.mH2O[-1], H2O, self.sys.P, name=self.sys.run.H2O_MODEL)*cnst.m['h2o'])
        self.h2.append(sl.h2_melt(gas.mH2[-1], H2, self.sys.P, self, name = self.sys.run.H2_MODEL)*cnst.m['h2'])
        
        if self.sys.run.GAS_SYS == "COHS" or self.sys.run.GAS_SYS == "COH" or self.sys.run.GAS_SYS == "COHSN":
            #self.co2.append(sl.co2_melt(gas.mCO2[-1], CO2, self.sys.P, name=self.sys.run.C_MODEL)*cnst.m['co2'])
            self.co2.append(sl.co2_melt((CO2.Y*gas.mCO2[-1]*self.sys.P), CO2, (O2.Y*gas.mO2[-1]*self.sys.P), self.sys.T, self.sys.P, self, name=self.sys.run.C_MODEL)*cnst.m['co2'])

            if self.sys.run.CO_MODEL != 'None':
                self.co.append(sl.co_melt((CO.Y*gas.mCO[-1]*self.sys.P), self.sys.P, name = self.sys.run.CO_MODEL)*cnst.m['co'])
            else:
                self.co.append(0.0)
            
            if self.sys.run.CH4_MODEL != 'None':
                self.ch4.append(sl.ch4_melt((CH4.Y*gas.mCH4[-1]*self.sys.P), self.sys.P, name = self.sys.run.CH4_MODEL)*cnst.m['ch4'])
            else:
                self.ch4.append(0.0)

            if self.graphite_sat == True:
                self.graphite.append(self.graph_current*cnst.m['c'])
            else:
                self.graphite.append(0.0)

        else:
            self.co2.append(0.0)
            self.co.append(0.0)
            self.ch4.append(0.0)
            self.graphite.append(0.0)

        if self.sys.run.GAS_SYS == "COHS" or self.sys.run.GAS_SYS == "SOH" or self.sys.run.GAS_SYS == "COHSN":
            fo2 = (O2.Y * gas.mO2[-1] * self.sys.P)
            self.cm_dry, F = self.iron_fraction(np.log(fo2))
            self.F.append(F)

            self.sulfide.append(sl.sulfide_melt((S2.Y * gas.mS2[-1] * self.sys.P), fo2, self.sys.P, self.sys.T, self, name=self.sys.run.SULFIDE_CAPACITY)*cnst.m['s'])
            self.sulfate.append(sl.sulfate_melt((S2.Y * gas.mS2[-1] * self.sys.P), fo2, self.sys.P, self.sys.T, self, self.sys.run, name=self.sys.run.SULFATE_CAPACITY)*cnst.m['s'])
            self.s.append(self.sulfide[-1] + self.sulfate[-1])
            
        else:
            fo2 = (O2.Y * gas.mO2[-1] * self.sys.P)
            self.cm_dry, F = self.iron_fraction(np.log(fo2))
            self.F.append(F)

            self.sulfide.append(0.0)
            self.sulfate.append(0.0)            
            self.s.append(0.0)
        
        if self.sys.run.GAS_SYS == "COHSN":
            self.n.append(sl.n_melt(gas.mN2[-1], (O2.Y*gas.mO2[-1]*self.sys.P), self.sys.P, name=self.sys.run.N_MODEL)*cnst.m['n'])
        else:
            self.n.append(0.0)
            
        self.cm_volatiles, self.cm_density = self.melt_dicts()

        self.rho_store.append(self.rho())  # Saves the melt density at each step          
            
        if self.sys.run.S_SAT_WARN == True and self.run.GAS_SYS in ['SOH', 'COHS', 'COHSN']:

            melt_mols = self.cm_volatiles   # The melt phase normalised to 1, including dissolved volatiles, as mol fractions

            # Sulfide + sulfate content of the melt in ppm
            S = sl.sulfide_melt((S2.Y * gas.mS2[-1] * self.sys.P), (O2.Y * gas.mO2[-1] * self.sys.P), self.sys.P, self.sys.T, self, name=self.sys.run.SULFIDE_CAPACITY)*cnst.m['s']*1e6
            S6 = sl.sulfate_melt((S2.Y * gas.mS2[-1] * self.sys.P), (O2.Y * gas.mO2[-1] * self.sys.P), self.sys.P, self.sys.T, self, self.sys.run, name=self.sys.run.SULFATE_CAPACITY)*cnst.m['s']*1e6

            SCSS = self.sulfide_sat(melt_mols)  # gives the sulphur content in ppm for sulphide saturation 

            if S > SCSS:
                msgs.scss_warn(SCSS, S, self.sys, gas, self)  # If sulfide content is above SCSS, writeout and exit.
            elif S6/(S+S6) >= 0.1:
                msgs.sulfate_warn(S6/(S+S6), self.sys, gas, self)   # If the S6 content is >10%, writes out current data and closes as so4 currently doesn't degass.

#  -------------------------------------------------------------------------
class Gas:
    """Stores the properties of the exsolved gas as it decompresses. Holds the results of each decompression step. """

    def __init__(self, sys):
        self.sys = sys
        self.fo2 = []
        self.mO2 = []   # Stores mO2 value after each pressure step
        self.mH2O = []
        self.mH2 = []
        self.mCO = []
        self.mCO2 = []
        self.mCH4 = []
        self.mSO2 = []
        self.mH2S = []
        self.mS2 = []
        self.mOCS = []
        self.mN2 = []
        self.atomicM = {}
        self.wt = {}    # Temporary store for the results of converting mole frac to weight frac
        self.Wt = {}    # Store of weight fractions at each pressure step
        self.f = {}     # Store the fugacity of a species at each step.
        self.M = []     # Stores the mean molecular weight of the gas phase


    def get_WgT(self, melt, mols):
        """Calculates the gas weight fraction. NOT PERCENT. """
        if self.sys.run.GAS_SYS == 'OH':
            H2O = mols
            wgt = ((self.sys.atomicM['o'] / cnst.m['o']) - sl.h2o_melt(self.mH2O[-1], H2O, self.sys.P, name=self.sys.run.H2O_MODEL))/((self.mH2O[-1] + 2 * self.mO2[-1])/ (self.mH2[-1]*cnst.m['h2']+self.mH2O[-1]*cnst.m['h2o']+self.mO2[-1]*cnst.m['o2']))

        elif self.sys.run.GAS_SYS == 'COH':
            H2O, O2, H2, CO, CO2, CH4 = mols
            lst = {'o2': self.mO2[-1], 'h2': self.mH2[-1], 'h2o': self.mH2O[-1], 'co': self.mCO[-1], 'co2': self.mCO2[-1], 'ch4': self.mCH4[-1]}
            
            mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

            wgt = ((self.sys.atomicM['c']/cnst.m['c'] - sl.co2_melt((CO2.Y*self.mCO2[-1]*self.sys.P), CO2, (O2.Y*self.mO2[-1]*self.sys.P), self.sys.T, self.sys.P, melt, name=self.sys.run.C_MODEL) - sl.co_melt((CO.Y*self.mCO[-1]*self.sys.P), self.sys.P, name=self.sys.run.CO_MODEL) - sl.ch4_melt((CH4.Y*self.mCH4[-1]*self.sys.P), self.sys.P, name=self.sys.run.CH4_MODEL) - melt.graph_current)/(self.mCO[-1] + self.mCO2[-1] + self.mCH4[-1])) * sum(mjMj)

        elif self.sys.run.GAS_SYS == 'SOH':
            H2O, O2, H2, S2, SO2, H2S = mols
            lst = {'o2': self.mO2[-1], 'h2': self.mH2[-1], 'h2o': self.mH2O[-1], 's2': self.mS2[-1],
                   'so2': self.mSO2[-1], 'h2s': self.mH2S[-1]}
            
            mjMj = [lst[ele] * cnst.m[ele] for ele in lst]
            
            wgt = ((self.sys.atomicM['h']/(2*cnst.m['h']) - (sl.h2o_melt(self.mH2O[-1], H2O, self.sys.P, name=self.sys.run.H2O_MODEL) + sl.h2_melt(self.mH2[-1], H2, self.sys.P, melt, name=self.sys.run.H2_MODEL)))/(self.mH2O[-1] + self.mH2[-1] + self.mH2S[-1])) * sum(mjMj)

        elif self.sys.run.GAS_SYS == 'COHS':
            H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S = mols
            lst = {'o2': self.mO2[-1], 'h2': self.mH2[-1], 'h2o': self.mH2O[-1], 'co': self.mCO[-1], 'co2': self.mCO2[-1], 'ch4': self.mCH4[-1], 's2': self.mS2[-1], 'so2': self.mSO2[-1], 'h2s': self.mH2S[-1]}
            
            mjMj = [lst[ele] * cnst.m[ele] for ele in lst]
            
            wgt = ((self.sys.atomicM['c']/cnst.m['c'] - sl.co2_melt((CO2.Y*self.mCO2[-1]*self.sys.P), CO2, (O2.Y*self.mO2[-1]*self.sys.P), self.sys.T, self.sys.P, melt, name=self.sys.run.C_MODEL) - sl.co_melt((CO.Y*self.mCO[-1]*self.sys.P), self.sys.P, name=self.sys.run.CO_MODEL) - sl.ch4_melt((CH4.Y*self.mCH4[-1]*self.sys.P), self.sys.P, name=self.sys.run.CH4_MODEL) - melt.graph_current)/(self.mCO[-1] + self.mCO2[-1] + self.mCH4[-1])) * sum(mjMj)

        elif self.sys.run.GAS_SYS == 'COHSN':
            H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2 = mols
            lst = {'o2': self.mO2[-1], 'h2': self.mH2[-1], 'h2o': self.mH2O[-1], 'co': self.mCO[-1], 'co2': self.mCO2[-1], 'ch4': self.mCH4[-1], 's2': self.mS2[-1], 'so2': self.mSO2[-1], 'h2s': self.mH2S[-1], 'n2':self.mN2[-1]}
            
            mjMj = [lst[ele] * cnst.m[ele] for ele in lst]
            
            wgt = ((self.sys.atomicM['c']/cnst.m['c'] - sl.co2_melt((CO2.Y*self.mCO2[-1]*self.sys.P), CO2, (O2.Y*self.mO2[-1]*self.sys.P), self.sys.T, self.sys.P, melt, name=self.sys.run.C_MODEL) - sl.co_melt((CO.Y*self.mCO[-1]*self.sys.P), self.sys.P, name=self.sys.run.CO_MODEL) - sl.ch4_melt((CH4.Y*self.mCH4[-1]*self.sys.P), self.sys.P, name=self.sys.run.CH4_MODEL) - melt.graph_current)/(self.mCO[-1] + self.mCO2[-1] + self.mCH4[-1])) * sum(mjMj)

        if wgt < 0.0:
            raise ValueError('Gas weight fraction is negative.')
        else:
            return wgt

    def get_vol_frac(self, melt):
        # EQ 15, D-Compress paper. Returns the gas volume fraction of the system (a).

        if self.sys.run.GAS_SYS == 'OH':
            self.wt = cnvs.mols2wts(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1])
            empty_lists = {'CO2':0.0, 'CO':0.0, 'CH4':0.0, 'S2':0.0, 'SO2':0.0, 'H2S':0.0, 'N2':0.0}
            M = cnvs.mean_mol_wt(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1])
        elif self.sys.run.GAS_SYS == 'COH':
            empty_lists = {'S2':0.0, 'SO2':0.0, 'H2S':0.0, 'N2':0.0}
            self.wt = cnvs.mols2wts(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1], CO=self.mCO[-1], CO2=self.mCO2[-1], CH4=self.mCH4[-1])
            M = cnvs.mean_mol_wt(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1], CO=self.mCO[-1], CO2=self.mCO2[-1], CH4=self.mCH4[-1])
        elif self.sys.run.GAS_SYS == 'SOH':
            empty_lists = {'CO2':0.0, 'CO':0.0, 'CH4':0.0, 'N2':0.0}
            self.wt = cnvs.mols2wts(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1], S2=self.mS2[-1], SO2=self.mSO2[-1], H2S=self.mH2S[-1])
            M = cnvs.mean_mol_wt(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1], S2=self.mS2[-1], SO2=self.mSO2[-1], H2S=self.mH2S[-1])
        elif self.sys.run.GAS_SYS == 'COHS':
            empty_lists = {'N2':0.0}
            self.wt = cnvs.mols2wts(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1], CO=self.mCO[-1], CO2=self.mCO2[-1], CH4=self.mCH4[-1], S2=self.mS2[-1], SO2=self.mSO2[-1], H2S=self.mH2S[-1])
            M = cnvs.mean_mol_wt(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1], CO=self.mCO[-1], CO2=self.mCO2[-1], CH4=self.mCH4[-1], S2=self.mS2[-1], SO2=self.mSO2[-1], H2S=self.mH2S[-1])
        elif self.sys.run.GAS_SYS == 'COHSN':
            empty_lists = None
            self.wt = cnvs.mols2wts(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1], CO=self.mCO[-1], CO2=self.mCO2[-1], CH4=self.mCH4[-1], S2=self.mS2[-1], SO2=self.mSO2[-1], H2S=self.mH2S[-1], N2=self.mN2[-1])
            M = cnvs.mean_mol_wt(H2O=self.mH2O[-1], O2=self.mO2[-1], H2=self.mH2[-1], CO=self.mCO[-1], CO2=self.mCO2[-1], CH4=self.mCH4[-1], S2=self.mS2[-1], SO2=self.mSO2[-1], H2S=self.mH2S[-1], N2=self.mN2[-1])


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
        
        GvF = 1 / (1 + ((M * self.sys.P * 1e5 * (1 - self.sys.WgT[-1])) / (cnst.R * self.sys.T * melt.rho(P=(self.sys.P*1e5)) * self.sys.WgT[-1])))

        return GvF

    def get_fugacity(self, molecules, gas_phase):
        """returns the fugacity of species in the gas phase.
        Data provided as a list of molecule instances (e.g. H2O) and their corresponding mol fraction values (e.g. gas.mH2O[-1])"""

        for m, g in zip(molecules, gas_phase):  # zips together an instance of the molecule class to it's corresponding list of stored mol frac values
            if isinstance(m, str):
                if m in self.f:
                    self.f[m].append(0.0)
                else:
                    self.f[m] = [0.0]
            
            else:
                if m.Mol in self.f:
                    self.f[m.Mol].append(m.Y * g * self.sys.P) # appends to the f dictionary for storage.
                else:
                    self.f[m.Mol] = [m.Y * g * self.sys.P]

    def rho(self):
        M = 0
        for mol in self.sys.SC:
            M += self.wt[mol] * cnst.m[mol.lower()]

        return ((M/1000)*self.sys.P*1e5) / (cnst.R * self.sys.T)

    def get_ys(self, system):
        # Calls the solvgas find_Y function to recalculate activity constants for each species in the system in question
        return sg.set_Y(self.sys, system)

    def get_atomic_mass(self, run, sys):

        if run.GAS_SYS == 'COHS' and sys.OCS == True:

            mO2 = self.mO2[-1]
            mH2O = self.mH2O[-1]
            mH2 = self.mH2[-1]
            mCO2 = self.mCO2[-1]
            mCO = self.mCO[-1]
            mCH4 = self.mCH4[-1]
            mSO2 = self.mSO2[-1]
            mS2 = self.mS2[-1]
            mH2S = self.mH2S[-1]
            
            lst = {'o2':mO2, 'h2':mH2, 'h2o':mH2O, 'co':mCO, 'co2':mCO2, 'ch4':mCH4, 'so2':mSO2, 's2':mS2, 'h2s':mH2S}
            mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

            # WtO
            self.atomicM['o'] = cnst.m['o'] * ((sys.WgT[-1] * mH2O) / sum(mjMj) +
                                            (2 * (sys.WgT[-1] * mO2) / sum(mjMj)) + 2 * (sys.WgT[-1] * mSO2) / sum(mjMj) +
                                            (sys.WgT[-1] * mCO) / sum(mjMj) + 2 * (sys.WgT[-1] * mCO2) / sum(mjMj))

            # WtH
            self.atomicM['h'] = 2 * cnst.m['h'] * ((sys.WgT[-1] * mH2O) / sum(mjMj) +
                                                (sys.WgT[-1] * mH2) / sum(mjMj) +
                                                (sys.WgT[-1] * mH2S) / sum(mjMj) +
                                                2*((sys.WgT[-1] * mCH4) / sum(mjMj)))

            # WtC
            self.atomicM['c'] = cnst.m['c'] * ((sys.WgT[-1] * mCO) / sum(mjMj) + (sys.WgT[-1] * mCO2) / sum(mjMj) +
                                                (sys.WgT[-1] * mCH4) / sum(mjMj))

            # WtS
            self.atomicM['s'] = cnst.m['s'] * ((2*sys.WgT[-1] * mS2) / sum(mjMj) + (sys.WgT[-1] * mH2S) / sum(mjMj) +
                                                (sys.WgT[-1] * mSO2) / sum(mjMj))

    def open_system(self, melt, fraction):
        # reduces the gas fraction by removing the amount designated by fraction. I.e. if fraction = 0.9, 10% of the gas fraction remains after the run.
        keep = 1 - fraction
        self.sys.WgT[-1] = self.sys.WgT[-1]*keep
        
        self.sys.atomicM['h'] = cnvs.atomicM_calc(self.sys, melt, self, 'h', -1, WgT=self.sys.WgT[-1])
        self.sys.atomicM['o'] = cnvs.atomicM_calc(self.sys, melt, self, 'o', -1, WgT=self.sys.WgT[-1])
        self.sys.atomicM['o_tot'] = cnvs.atomicM_calc(self.sys, melt, self, 'o_tot', -1, WgT=self.sys.WgT[-1])
        self.sys.atomicM['c'] = cnvs.atomicM_calc(self.sys, melt, self, 'c', -1, WgT=self.sys.WgT[-1])
        self.sys.atomicM['s'] = cnvs.atomicM_calc(self.sys, melt, self, 's', -1, WgT=self.sys.WgT[-1])
        self.sys.atomicM['n'] = cnvs.atomicM_calc(self.sys, melt, self, 'n', -1, WgT=self.sys.WgT[-1])


#  -------------------------------------------------------------------------
class Output:
    """To request the desired outputs through the tmp.out file"""

    n_par = int(6)

    def __init__(self):
        # Graphical outputs
        self.Plt_melt_species = True
        self.Plt_gas_species_wt = True
        self.Plt_gas_species_mol = True
        self.Plt_gas_ratios = []
        self.Plt_gas_fraction = True
        self.Plt_fo2_dFMQ = False

    def set_outputs(self, params):

            for par, val in params.items():
                if self.zTest_Params(par):
                    setattr(self, par, val)

                else:
                    sys.exit("Warning: %s not a valid calculation parameter" % par)

    # helper methods
    def zTest_Params(self, par):

        for item in inspect.getmembers(self):
            if item[0] == par:
                return True
        return 0


    def print_conditions(self):

        l = 0
        for item in inspect.getmembers(self):
            if not inspect.ismethod(self):
                if l < self.n_par:
                    print(item[0], ':', item[1])
                    l += 1
                else:
                    break
        print("\n")
