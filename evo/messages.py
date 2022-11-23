""" Stores warning messages used thought the system to alert the user of an input error or failure to run to completion, and an explaination."""
import numpy as np
import sys

from evo import constants as cnst
from evo import conversions as cnvs
from evo import solubility_laws as sl
from evo import writeout as wt

def open_step():
    exit("Warning: Open system degassing is path dependent. \nFor internal consistency, please set the max (DP_MAX) and min (DP_MIN) pressure steps to be equal (we suggest <1 bar for the run to complete; 0.5 bar is usually sufficent)")

def step_size():
    exit('Error: Please make the minimum stepsize DP_MIN <= maximum stepsize DP_MAX.')

def open_earlyexit(sys, gas, melt): 
    wt.writeout_crash(sys, gas, melt, np.arange(sys.run.P_START, sys.P, sys.run.DP_MAX*-1.0))
    print('Error: This run has finished before the assigned P_STOP, as one or more volatile element is now present in too low an amount to solve the system with the specified pressure step.\nPlease try reducing the stepsize for a complete run.\nThe run so far has been written out in Output/dgs_output.csv')

def closed_earlyexit(sys, gas, melt):
    
    wt.writeout_crash(sys, gas, melt, sys.P_track)

    print('Error: This run has finished before the assigned P_STOP.\n The run so far has been written out in Output/dgs_output.csv')

def earlyexit(sys, gas, melt, msg):
    if sys.run.RUN_TYPE == 'open':
        open_earlyexit(sys, gas, melt)
    elif sys.run.RUN_TYPE == 'closed':
        closed_earlyexit(sys, gas, melt)
    raise AssertionError(msg + ' Try reducing the maximum pressure stepsize (DP_MAX) for model convergence.')
    
def runtime_error(sys, gas, melt, msg):
    if sys.run.RUN_TYPE == 'open':
        open_earlyexit(sys, gas, melt)
    elif sys.run.RUN_TYPE == 'closed':
        closed_earlyexit(sys, gas, melt)
    raise RuntimeError(msg)

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y]/n "
    elif default == "no":
        prompt = " y/[N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

def run_chem_mismatch(mtype, lims, sio2):
    if sio2 < lims[0]:
        print(f'Warning: You specified {mtype.upper()} in the env.yaml file, but the SiO2 content of your melt specified in the chem.yaml file is LESS THAN {lims[0]} wt%.\n')
        answer = query_yes_no('Are you sure you want to continue?', default='no')
    elif sio2> lims[1]:
        print(f'Warning: You specified {mtype.upper()} in the env.yaml file, but the SiO2 content of your melt specified in the chem.yaml file is GREATER THAN {lims[1]} wt%.\n')
        answer = query_yes_no('Are you sure you want to continue?', default='no')

    if answer == False:
        exit('Exiting...')

def solubility_temp(T, limits):
    
    print(f'Warning: The set temperature of {T:3Ng} K is outside the T-dependent solubility limits of {limits[0]} to {limits[1]} K.\n')
    answer = query_yes_no('Would you like to proceed with fixed solubility coefficients?', default='yes')

    if answer == False:
        exit('Exiting...')

def sulfate_warn(s6, sys, gas, melt):
    # Finish off this run
    sys.GvF.append(gas.get_vol_frac(melt))
    sys.rho.append(sys.rho_bulk(melt, gas))
    gas.M.append(cnvs.mean_mol_wt(H2O=gas.mH2O[-1], O2=gas.mO2[-1], H2=gas.mH2[-1], CO=gas.mCO[-1], CO2=gas.mCO2[-1], CH4=gas.mCH4[-1], S2=gas.mS2[-1], SO2=gas.mSO2[-1], H2S=gas.mH2S[-1]))
    
    if not sys.P_track:  # If sys.P hasn't been initiated yet, fill with the starting pressure
        sys.P_track.append(sys.P)
    
    # Writeout results

    closed_earlyexit(sys, gas, melt)

    exit(f"Warning: The sulfate (S6+ speciated as SO4) content of the melt is {s6*100:.1f}% of the total melt S content.\nAs sulfate currently does not exsolve, gas phase sulphur-bearing species will now be significantly underestimated.\nRun is terminated here, results up to this point have been exported to CSV.")

def scss_warn(scss, s, sys, gas, melt):
    # Finish off this run
    sys.GvF.append(gas.get_vol_frac(melt))
    sys.rho.append(sys.rho_bulk(melt, gas))
    gas.M.append(cnvs.mean_mol_wt(H2O=gas.mH2O[-1], O2=gas.mO2[-1], H2=gas.mH2[-1], CO=gas.mCO[-1], CO2=gas.mCO2[-1], CH4=gas.mCH4[-1], S2=gas.mS2[-1], SO2=gas.mSO2[-1], H2S=gas.mH2S[-1]))
    
    if not sys.P_track:  # If sys.P hasn't been initiated yet, fill with the starting pressure
        sys.P_track.append(sys.P)
    
    # Writeout results

    closed_earlyexit(sys, gas, melt)

    exit(f'You have reached sulphide saturation, this code is no longer valid. SCSS: {scss:.1g} ppm, Sulfide in melt: {s:.1g} ppm\nRun is terminated here, results up to this point have been exported to CSV.')

def undersat():
    exit('Error: The system is volatile undersaturated at this pressure. Please lower the starting pressure and try again.')

def vol_setup_standard(run, sys):
    "Checks that the conditions set in Env file match the run type specified."
    
    # check some volatile species aren't being set twice as wt% and fugacity
    if (run.FH2O_SET == True and run.WTH2O_SET == True) or (run.FCO2_SET == True and run.WTCO2_SET == True):
        exit("Warning: Please only use one option to set H2O or CO2. Set either the Wt% or fugacity value to False.")

    if (run.WTCO2_SET == True and run.GRAPHITE_SATURATED == True) or (run.FCO2_SET == True and run.GRAPHITE_SATURATED == True):
        exit("Warning: Please only use one setting for CO2. Set two from fCO2, WTCO2 and graphite_saturated to False.")
    
    if run.ATOMIC_MASS_SET == True:
        lst = {'FH2_SET':run.FH2_SET, 'FH2O_SET':run.FH2O_SET, 'FCO2_SET':run.FCO2_SET, 'WTH2O_SET':run.WTH2O_SET, 
                'WTCO2_SET':run.WTCO2_SET, 'SULFUR_SET':run.SULFUR_SET, 'NITROGEN_SET':run.NITROGEN_SET, 'GRAPHITE_SATURATED':run.GRAPHITE_SATURATED}
        if any(setting == True for setting in lst.values()):
            print(f'Warning: Volatile contents are being set as fugacities/melt contents while ATOMIC_MASS_SET is true.\n')
            
            answer = query_yes_no('Would you like to proceed with the atomic mass set, setting all fugacities and melt contents to False?', default='yes')

            if answer == True:
                for setting in lst.keys():
                    setattr(run, setting, False)
            
            if answer == False:
                answer = query_yes_no('Would you like to set ATOMIC_MASS_SET to False and continue with the fugacities/melt contents already set?', default='yes')
                
                if answer == True:
                    run.ATOMIC_MASS_SET = False

                if answer == False:
                    exit('Exiting...')
    
    elif run.GAS_SYS == 'OH':

        lst = [run.FH2O_SET, run.WTH2O_SET, run.WTCO2_SET, run.FCO2_SET, run.SULFUR_SET]  # PL: Make this a yes/no option to change setup to sat point finder with H2O.
        for opt in lst:
            if opt == True:
                exit(f"Error: This OH system setup cannot run with H2O, CO2 or sulphur contents. To run with H2O and find the saturation point of the system, set run.FIND_SATURATION to True.")
        
        if run.FH2_SET == True and sys.FO2:
            exit("Warning: Please only set one of either fH2, fO2, or the proportions of FeO and Fe2O3 for the OH system.")
        
        elif run.FH2_SET == False and sys.FO2 == None:
            exit("Error: For the OH system, please set either the initial fH2 or fO2 in the env.yaml file, or provide proportions of FeO and Fe2O3 via the chem.yaml file.")

    elif run.GAS_SYS == 'COH' or run.GAS_SYS == 'SOH':
        
        if run.FCO2_SET == True or run.WTCO2_SET == True or run.SULFUR_SET == True:
            exit(f"Warning: CO2 and S cannot be used as an initial condition for the {run.GAS_SYS} system. \nPlease set 2 initial conditions from either fH2, H2O (as fH2O or wt% H2O) or fO2 (either directly or via FeO/Fe2O3 proportions).\nTo run with the complete melt volatile content and find the saturation point of the system, set run.FIND_SATURATION to True.")   
        
        lst = [run.FH2O_SET, run.WTH2O_SET, run.FH2_SET]  # Allowable setup conditions
        true = 0
        if sys.FO2: # final allowable condition, but present as a boolean value
            true += 1
        for ele in lst:
            if ele == True:
                true += 1
        if true > 2 or true <= 1:
            exit("Error: Please set 2 initial conditions from either fH2, H2O (as fH2O or wt% H2O) or fO2 (either directly or via FeO/Fe2O3 proportions)."
                 f"\nThere are currently {true} values set.")

    elif run.GAS_SYS == 'COHS' or run.GAS_SYS == 'COHSN':
        lst = [run.FH2O_SET, run.WTH2O_SET, run.FH2_SET, run.FCO2_SET, run.WTCO2_SET, run.SULFUR_SET, run.GRAPHITE_SATURATED]
        true = 0
        if sys.FO2:
            true += 1
        for ele in lst:
            if ele == True:
                true += 1
        if true > 3 or true <= 2:
            exit("Error: Please set 3 initial conditions from either fH2, H2O (as fH2O or wt% H2O), CO2 (as fCO2, wt% CO2 or graphite), S (as wt% S) or fO2 (either directly or via FeO/Fe2O3 proportions)."
                 f"\nThere are currently {true} values set.")
        
        if run.GAS_SYS == 'COHSN':
            if run.NITROGEN_SET == False:
                exit("Error: Please set a value for the melt nitrogen content using NITROGEN_SET/START.")
        
        elif run.NITROGEN_SET == True: # and system is COHS
                
            print(f'Warning: Nitrogen is being set in the melt, but the system is running with COHS (i.e. no nitrogen).\n')
            
            answer = query_yes_no('Would you like to proceed with the COHS system, ignoring the melt nitrogen content?', default='yes')

            if answer == False:
                exit('Exiting...')        


def vol_setup_saturation(run, sys):
    "Checks that the conditions set in Env file match the run type specified when a finding the saturation point is specified."

    # check some volatile species aren't being set twice as wt% and fugacity
    if (run.FH2O_SET == True and run.WTH2O_SET == True) or (run.FCO2_SET == True and run.WTCO2_SET == True):
        exit("Warning: Please only use one option to set H2O or CO2. Set either the Wt% or fugacity value to False.")

    if (run.WTCO2_SET == True and run.GRAPHITE_SATURATED == True) or (run.FCO2_SET == True and run.GRAPHITE_SATURATED == True):
        exit("Warning: Please only use one setting for CO2. Set two from fCO2, WTCO2 and graphite_saturated to False.")
   
    if run.GAS_SYS == 'OH':
        lst = [run.WTCO2_SET, run.FCO2_SET, run.SULFUR_SET]  # PL: Make this a yes/no option to change setup to sat point finder with H2O.
        for opt in lst:
            if opt == True:
                exit(f"Error: This OH system setup cannot run with CO2 or sulphur contents. Set all CO2 and S values to False.")
        
        if run.FH2_SET == True or run.FH2O_SET== True:
            exit("Error: Please only set fO2 and melt volatile contents to find the volatile saturation point. Turn all fugacities to False in the env.yaml file.")
        
        if not (run.WTH2O_SET == True and sys.FO2):
            exit("Error: To find the saturation point for the OH system, please set both the melt water content, and either fO2 in the env.yaml file, or provide proportions of FeO and Fe2O3 via the chem.yaml file.")

    elif run.GAS_SYS == 'COH':
        lst = [run.FH2O_SET, run.FCO2_SET, run.FH2_SET, run.SULFUR_SET]
        for opt in lst:
            if opt == True:
                exit(f"Error: Please only set fO2 and melt H2O & CO2 contents to find the volatile saturation point. Turn all fugacities and sulfur contents to False.")
        
        if run.WTCO2_SET == True or run.GRAPHITE_SATURATED == True:
            C = True
        else:
            C = False

        if not (run.WTH2O_SET == True and C ==True and sys.FO2):
            exit("Error: To find the saturation point for the COH system, please set both water and CO2 contents in the melt, and either fO2 in the env.yaml file, or provide proportions of FeO and Fe2O3 via the chem.yaml file.")
    
    elif run.GAS_SYS == 'SOH':
        lst = [run.FH2O_SET, run.FCO2_SET, run.WTCO2_SET, run.FH2_SET]
        for opt in lst:
            if opt == True:
                exit(f"Error: Please only set fO2 and melt H2O & sulphur contents to find the volatile saturation point. Turn all fugacities and CO2 contents to False.")
        
        if not (run.WTH2O_SET == True and run.SULFUR_SET==True and sys.FO2):
            exit("Error: To find the saturation point for the SOH system, please set both water and sulphur contents in the melt, and either fO2 in the env.yaml file, or provide proportions of FeO and Fe2O3 via the chem.yaml file.")
    
    elif (run.GAS_SYS == 'COHS' or run.GAS_SYS == 'COHSN'):
        if run.WTCO2_SET == True or run.GRAPHITE_SATURATED == True:
            C = True
        
        lst = [run.WTH2O_SET, C, run.SULFUR_SET, sys.FO2]

        if not all(lst):  # if any are false
            exit("Error: To find the volatile saturation point please make sure the melt water, CO2 and sulfur contents are set and 'True' in env.yaml.")

        if any([run.FH2O_SET, run.FH2_SET, run.FCO2_SET]): # if any are true
            exit('Error: Please only set melt weight fractions, not gas phase fugacities, when searching for the saturation point.')
        
        if run.GAS_SYS == 'COHSN':
            if run.NITROGEN_SET == False:
                exit("Error: Please set a value for the melt nitrogen content using NITROGEN_SET/START.")
                
        elif run.NITROGEN_SET == True: # and system is COHS
                
            print(f'Warning: Nitrogen is being set in the melt, but the system is running with COHS (i.e. no nitrogen).\n')
            
            answer = query_yes_no('Would you like to proceed with the COHS system, ignoring the melt nitrogen content?', default='yes')

            if answer == False:
                exit('Exiting...')
            else:
                run.NITROGEN_SET = False

def graphite_warn(melt):
    """
    Warns of graphite saturation without knowledge of the size of the graphite reservoir if melt_co2 is set but melt is graphite saturated.

    Returns True is the graphite content is set internally, returns false if it still needs to be calculated.
    """
    print(f'Warning: The melt is graphite saturated.\n')
    
    answer = query_yes_no('Does the CO2 content set include the graphite content (melt_co2 + melt_co3 + melt_graphite)?', default='yes')

    if answer == True:
        print('Continuing to find melt carbonate/graphite split.')
        return False

    elif answer == False:
        answer = query_yes_no('Do you know the graphite content of the melt?', default='no')

        if answer == True:
            while True:
                try:
                    graph = float(input("Please enter the graphite content of the melt as a weight fraction: "))
                except ValueError:
                    print("Re-enter value as a float: ")
                    continue

                if graph < 0:
                    print("Sorry, your response must not be negative.")
                    continue
                elif graph> 1:
                    print("Please enter as a weight fraction - your response must be less than 1.")
                else:
                    break
            
            melt.graphite_sat = True
            melt.graph_current = graph/cnst.m['c']

            print('By continuing, the melt CO2 content will be over_written based on the melt CO2 content at graphite saturation.')
            answer = query_yes_no(f'Are you happy to continue with a melt graphite content of {graph*100} wt%?', default='yes')
            
            if answer == True:
                return True
            else:
                exit('Exiting...')

        elif answer == False:
            exit('Error: Cannot continue with this setup. Please either set up the melt so that the total amount of C in the system is known (using atomic mass set),\n\
                or use a different solubility law for CO2 which neglects graphite saturation.')

def graphite_warn_saturation(melt, fCO2, fO2, CO2):
    """
    Warns of graphite saturation in saturation solver & provides options to continue.
    """

    print('Warning: Melt is graphite saturated at the volatile saturation pressure.\nRun cannot proceed without knowing the size of the graphite reservoir at volatile saturation.\n')
    
    answer = query_yes_no('Does the melt CO2 content at saturation include the graphite content (melt_co2 + melt_co3 + melt_graphite)?', default='yes')

    if answer == True:
        co2_melt = sl.co2_melt(fCO2, CO2, fO2, melt.sys.T, melt.sys.P, melt, name = melt.sys.run.C_MODEL)*cnst.m['co2']
        graph_melt = ((melt.sys.run.WTCO2_START - co2_melt)/cnst.m['co2'])*cnst.m['c']   # wt frac graphite in melt
        melt.graph_current = graph_melt/cnst.m['c']

        print(f'The original melt CO2 content was {melt.sys.run.WTCO2_START*100} wt%.')
        print(f'At graphite-present volatile saturation, the melt CO2 content is {co2_melt*100} wt% and the graphite content is {graph_melt*100} wt%.')

        answer = query_yes_no('Do you wish to continue?', default = 'yes')

        if answer == True:
            print('Continuing with calculation')

        elif answer == False:
            exit('Exiting...')
    
    elif answer == False:
        exit('Error: Run cannot proceed without knowing the size of the graphite reservoir at volatile saturation.\nExiting...')

def fo2_model_mismatch(feot, run):

    names = {'r2013':'Righter et al. (2013)', 'kc1991':'Kress & Carmicheal (1991)'}
    
    if feot <= 15.0:
        print(f'Warning: The approximate FeO(t) content of the melt is {feot:.4g} wt%. This is out of the {names[run.FO2_MODEL]} model bounds.')
        answer = query_yes_no('Would you like to switch to the Kress & Carmicheal (1991) fO2 model?', default='yes')

        if answer == True:
            run.FO2_MODEL = 'kc1991'
        elif answer == False:
            print(f"Warning: Continuing with the {names[run.FO2_MODEL]} fO2 model outside of calibration range.")

    elif feot > 15.0:
        print(f'Warning: The approximate FeO(t) content of the melt is {feot:.4g} wt%. The Righter et al. (2013) fO2 model is more appropriate for high FeO melts.')
        answer = query_yes_no('Would you like to switch to the Righter et al. (2013) fO2 model?', default='yes')

        if answer == True:
            print('\n')
            run.FO2_MODEL = 'r2013'
        elif answer == False:
            print(f"Warning: Continuing with the {names[run.FO2_MODEL]} fO2 model outside of calibration range.\n")