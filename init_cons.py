"""
Finds the atomic mass fraction of the system if they haven't already been set (ATOMIC_MASS_SET = True)
and the saturation pressure hasn't been asked for (FIND_SATURATION=TRUE).
"""
import numpy as np
from numpy import sqrt

import constants as cnst
import messages as msgs
import solubility_laws as sl

def get_inputs(sys, run, melt, gas, mols) -> tuple[float, ...]:
    """
    Calculates a subset of gas mole fractions based on the provided
    input parameters.

    Calculates the gas mole fractions corresponding to the initial
    provided parameters. E.g., calculates the mole fraction of H2O if the
    weight fraction of H2O in the melt is provided as an input parameter.

    If the fO2 is calculated at this point, the ferric/ferrous ratio in
    the melt is also set here.

    Parameters
    ----------
    sys : ThermoSystem class
        The active instance of the ThermoSystem class
    run : RunDef class
        The active instance of the RunDef class    
    melt : Melt class
        The active instance of the Melt class
    gas : Gas class
        The active instance of the Gas class
    mols : tuple of Molecule classes
        A tuple of all the Molecule class instances. Setup to be passed
        through all functions, preserving ordering.

    Returns
    -------
    tuple of floats
        mole fractions of H2O, O2, H2, or CO2, H2O, O2 and H2 is the system is COHS(N).
    """

    if run.GAS_SYS == 'COH':
        H2O, O2, H2, CO, CO2, CH4 = mols
    elif run.GAS_SYS == 'SOH':
        H2O, O2, H2, S2, SO2, H2S = mols
    elif run.GAS_SYS == "COHS":
        H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S = mols
    elif run.GAS_SYS == "COHSN":
        H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2 = mols

    if run.FH2_SET == True:
        mH2 = run.FH2_START / (H2.Y * sys.P)
        gas.mH2.append(mH2)

        if sys.FO2 != None:  # if run.FO2 set was true, this is set in the readinenv file. If Fe ratio was set, it's created in the melt init from the melt chemistry
            mO2 = np.exp(sys.FO2) / (O2.Y * sys.P)
            gas.mO2.append(mO2)

            mH2O = (sys.K['K1']*H2.Y*mH2*(O2.Y*mO2*sys.P)**0.5)/H2O.Y
            gas.mH2O.append(mH2O)

        elif run.FH2O_SET == True:
            mH2O = run.FH2O_START / (H2O.Y * sys.P)
            gas.mH2O.append(mH2O)

            mO2 = (((H2O.Y*mH2O*sys.P) / ((H2.Y*mH2*sys.P) * sys.K['K1'])) ** 2) / (O2.Y * sys.P)
            gas.mO2.append(mO2)
            sys.FO2 = np.log(O2.Y * mO2 * sys.P)

        elif run.WTH2O_SET == True:
            mH2O = sl.h2o_fugacity(run.WTH2O_START, H2O, name=run.H2O_MODEL)/(H2O.Y*sys.P)
            assert H2O.Y*mH2O*sys.P > 0, 'H2O fugacity is negative!'
            gas.mH2O.append(mH2O)


            mO2 = (((H2O.Y*mH2O*sys.P) / ((H2.Y*mH2*sys.P) * sys.K['K1'])) ** 2) / (O2.Y * sys.P)
            gas.mO2.append(mO2)
            sys.FO2 = np.log(O2.Y * mO2 * sys.P)

    elif sys.FO2:
        mO2 = np.exp(sys.FO2) / (O2.Y * sys.P)
        gas.mO2.append(mO2)

        if run.FH2O_SET == True:
            mH2O = run.FH2O_START / (H2O.Y * sys.P)
            gas.mH2O.append(mH2O)

        elif run.WTH2O_SET == True:
            mH2O = sl.h2o_fugacity(run.WTH2O_START, H2O, name=run.H2O_MODEL)/(H2O.Y*sys.P)
            assert H2O.Y * mH2O * sys.P > 0, 'fH2O is negative!'
            gas.mH2O.append(mH2O)


        mH2 = (H2O.Y*mH2O*sys.P) / (sys.K['K1']*H2.Y*sys.P*np.exp(sys.FO2)**0.5)
        assert H2.Y*mH2*sys.P > 0, 'fH2 is negative!'
        gas.mH2.append(mH2)

    if not run.FCO2_SET and not run.WTCO2_SET and not run.SULFUR_SET and not run.GRAPHITE_SATURATED:
        
        melt.cm_dry = melt.iron_fraction(sys.FO2)[0] # Uses fO2 to set the Fe2/Fe3 ratio and dry melt chemistry prior to needing sulfide capacity
        
        return mH2O, mO2, mH2

    elif run.SULFUR_SET:  # This comes first
        melt.cm_dry = melt.iron_fraction(sys.FO2)[0] # Uses fO2 to set the Fe2/Fe3 ratio and dry melt chemistry prior to needing sulfide capacity
        
        mS2 = sl.S2_fugacity(run.SULFUR_START, (O2.Y*mO2*sys.P), sys.P, sys.T, melt, run, sulfidename=run.SULFIDE_CAPACITY, sulfatename=run.SULFATE_CAPACITY)/(S2.Y*sys.P)
        gas.mS2.append(mS2)
        
        return mS2, mH2O, mO2, mH2
    
    elif run.FCO2_SET:
        mCO2 = run.FCO2_START / (CO2.Y * sys.P)
        gas.mCO2.append(mCO2)

        melt.cm_dry = melt.iron_fraction(sys.FO2)[0] # Uses fO2 to set the Fe2/Fe3 ratio and dry melt chemistry prior to needing sulfide capacity
    
        return mCO2, mH2O, mO2, mH2        
    
    elif run.WTCO2_SET:
        if run.C_MODEL == 'eguchi2018':
            fCO2 = sl.co2_fugacity(run.WTCO2_START, CO2, np.exp(sys.FO2), sys.T, sys.P, melt, name=run.C_MODEL)
            if melt.graphite_sat == True:
                found = msgs.graphite_warn(melt)                
                if found == False:
                    mCO2 = sl.graphite_fco2(sys.T, sys.P, (O2.Y*mO2*sys.P))/(CO2.Y*sys.P)
                    gas.mCO2.append(mCO2)
                    melt.cm_dry = melt.iron_fraction(sys.FO2)[0]

                    co2_melt = sl.co2_melt((CO2.Y*mCO2*sys.P), CO2, (O2.Y*mO2*sys.P), sys.T, sys.P, melt, name = run.C_MODEL)*cnst.m['co2']
                    graph_melt = ((run.WTCO2_START - co2_melt)/cnst.m['co2'])*cnst.m['c']   # wt frac graphite in melt
                    melt.graphite_sat = True
                    melt.graph_current = graph_melt/cnst.m['c']
                    run.WTCO2_START = co2_melt
                
                elif found == True:
                    mCO2 = sl.graphite_fco2(sys.T, sys.P, (O2.Y*mO2*sys.P))/(CO2.Y*sys.P)
                    gas.mCO2.append(mCO2)
                    melt.cm_dry = melt.iron_fraction(sys.FO2)[0]

                    co2_melt = sl.co2_melt((CO2.Y*mCO2*sys.P), CO2, (O2.Y*mO2*sys.P), sys.T, sys.P, melt, name = run.C_MODEL)*cnst.m['co2']
                    run.WTCO2_START = co2_melt

            elif melt.graphite_sat == False:
                mCO2 = fCO2/(CO2.Y*sys.P)
                gas.mCO2.append(mCO2)

                melt.cm_dry = melt.iron_fraction(sys.FO2)[0]

        else:
            mCO2 = sl.co2_fugacity(run.WTCO2_START, CO2, np.exp(sys.FO2), sys.T, sys.P, melt, name=run.C_MODEL)/(CO2.Y*sys.P)
            gas.mCO2.append(mCO2)

            melt.cm_dry = melt.iron_fraction(sys.FO2)[0] # Uses fO2 to set the Fe2/Fe3 ratio and dry melt chemistry prior to needing sulfide capacity
    
        return mCO2, mH2O, mO2, mH2    

    elif run.GRAPHITE_SATURATED == True:
        melt.graphite_sat = True
        melt.graph_current = run.GRAPHITE_START/cnst.m['c']

        mCO2 = sl.graphite_fco2(sys.T, sys.P, (O2.Y*mO2*sys.P))/(CO2.Y*sys.P)
        gas.mCO2.append(mCO2)
        melt.cm_dry = melt.iron_fraction(sys.FO2)[0]

        return mCO2, mH2O, mO2, mH2

def oh(sys, run, melt, gas, mols):
    """
    Calculates the initial conditions for the OH system.

    Based on the subset of parameters provided at input, where the starting
    pressure is known, the speciation of the gas phase as mole fractions is
    calculated. Based on this data, the volatile content of the melt is
    found and the masses of each volatile element (O, C, H etc) is set,
    ready for use in the equilibrium constant and mass-balance solution.

    Parameters
    ----------
    sys : ThermoSystem class
        The active instance of the ThermoSystem class
    run : RunDef class
        The active instance of the RunDef class    
    melt : Melt class
        The active instance of the Melt class
    gas : Gas class
        The active instance of the Gas class
    mols : tuple of Molecule classes
        A tuple of all the Molecule class instances. Setup to be passed
        through all functions, preserving ordering.

    Returns
    -------
    sys.atomicM : dict
        A dictionary containing the mass fraction of each volatile element
        in the system.
    """
    H2O, O2, H2 = mols

    if run.FH2_SET == True:
        fH2 = run.FH2_START
        assert fH2 > 0, 'fH2 is negative!'
        mH2 = fH2 / (H2.Y * sys.P)
        gas.mH2.append(mH2)

        def quadratic():

            a = 1/(O2.Y*sys.P)
            b = ((sys.K['K1']*H2.Y*mH2)/H2O.Y)**2
            c = -((1-mH2**2)/(O2.Y*sys.P))

            mo2 = [x for x in np.roots([a, b, c]) if x > 0]
            return mo2[0]
            
        mO2 = quadratic()
        
        gas.mO2.append(mO2)
        sys.FO2 = np.log(O2.Y * mO2 * sys.P)

    else:  # assumes FO2 has either been set by the user, or calculated from the iron ratio in an earlier step
        fO2 = np.exp(sys.FO2)
        mO2 = fO2 / (O2.Y * sys.P)
        gas.mO2.append(mO2)

        mH2 = ((1. - mO2) * H2O.Y) / (sys.K['K1'] * H2.Y * (fO2 ** 0.5) + H2O.Y)
        gas.mH2.append(mH2)

    mH2O = 1 - mH2 - mO2
    gas.mH2O.append(mH2O)

    empty_lists = [gas.mCO2, gas.mCO, gas.mCH4, gas.mSO2, gas.mS2, gas.mH2S, gas.mN2]
    for alist in empty_lists:
        alist.append(np.float64(0))

    lst = {'o2':mO2, 'h2':mH2, 'h2o':mH2O}
    mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

    # WtO
    sys.atomicM['o'] = cnst.m['o'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + 2*mO2)
                                     + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL))

    # WtH(total)
    sys.atomicM['h'] = 2 * cnst.m['h'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + mH2) +
                                            sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL) +
                                            sl.h2_melt(mH2, H2, sys.P, melt, name=run.H2_MODEL))
    return sys.atomicM

def coh(sys, run, melt, gas, mols):
    """
    Calculates the initial conditions for the COH system.

    Based on the subset of parameters provided at input, where the starting
    pressure is known, the speciation of the gas phase as mole fractions is
    calculated. Based on this data, the volatile content of the melt is
    found and the masses of each volatile element (O, C, H etc) is set,
    ready for use in the equilibrium constant and mass-balance solution.

    Parameters
    ----------
    sys : ThermoSystem class
        The active instance of the ThermoSystem class
    run : RunDef class
        The active instance of the RunDef class    
    melt : Melt class
        The active instance of the Melt class
    gas : Gas class
        The active instance of the Gas class
    mols : tuple of Molecule classes
        A tuple of all the Molecule class instances. Setup to be passed
        through all functions, preserving ordering.

    Returns
    -------
    sys.atomicM : dict
        A dictionary containing the mass fraction of each volatile element
        in the system.
    """
    H2O, O2, H2, CO, CO2, CH4 = mols

    # uses get_inputs to find which parameters have been set and get the water system.
    mH2O, mO2, mH2 = get_inputs(sys, run, melt, gas, mols)

    # Solve for CH4

    assert 1 > (1-mH2O-mH2-mO2) > 0, 'Some fugacities are negative!'

    mCH4 = (1 - mH2 - mO2 - mH2O)/(1 + ((sys.K['K3']*CH4.Y*(O2.Y*mO2)**2)/(CO2.Y*(H2O.Y*mH2O)**2)) + ((sys.K['K3']*CH4.Y*(O2.Y*mO2)**2)/(sys.K['K2']*CO.Y*(O2.Y*mO2*sys.P)**0.5*(H2O.Y*mH2O)**2)))
    gas.mCH4.append(mCH4)

    mCO2 = (sys.K['K2'] * (O2.Y*mO2*sys.P)**0.5 * CO.Y * (1 - mH2 - mO2 - mH2O - mCH4)) / (CO2.Y + sys.K['K2'] * (O2.Y*mO2*sys.P)**0.5 * CO.Y)
    gas.mCO2.append(mCO2)

    # Check graphite saturation
    if run.C_MODEL == 'eguchi2018':
        graph_fCO2 = sl.graphite_fco2(sys.T, sys.P, (O2.Y*mO2*sys.P))
        if (CO2.Y*mCO2*sys.P) > graph_fCO2:
            exit('Error: Melt is graphite saturated. System is overconstrained; please rerun by setting melt graphite content and finding the saturation pressure, or stop using eguchi2018 in C_MODEL.')

    mCO = 1 - mH2O - mH2 - mO2 - mCH4 - mCO2
    gas.mCO.append(mCO)

    empty_lists = [gas.mSO2, gas.mS2, gas.mH2S, gas.mN2]
    for list in empty_lists:
        list.append(np.float64(0))

    lst = {'o2':mO2, 'h2':mH2, 'h2o':mH2O, 'co':mCO, 'co2':mCO2, 'ch4':mCH4}
    mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

    # WtO
    sys.atomicM['o'] = cnst.m['o'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + 2*mO2 + 2*mCO2 + mCO)
                                    + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
                                    + 2*sl.co2_melt((CO2.Y*mCO2*sys.P), CO2, (O2.Y*mO2*sys.P), sys.T, sys.P, melt, name=run.C_MODEL)
                                    + sl.co_melt((CO.Y*mCO*sys.P), sys.P, name = run.CO_MODEL))

    # WtH(total)
    sys.atomicM['h'] = 2 * cnst.m['h'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + mH2 + 2*mCH4) +
                                        sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL) +
                                        sl.h2_melt(mH2, H2, sys.P, melt, name=run.H2_MODEL) +
                                        2*sl.ch4_melt((CH4.Y*mCH4*sys.P), sys.P, name = run.CH4_MODEL))

    # WtC
    sys.atomicM['c'] = cnst.m['c'] * ((sys.WgT[0] / sum(mjMj)) * (mCO + mCO2 + mCH4) +
                                        sl.co2_melt((CO2.Y*mCO2*sys.P), CO2, (O2.Y*mO2*sys.P), sys.T, sys.P, melt, name=run.C_MODEL) +
                                        sl.ch4_melt((CH4.Y*mCH4*sys.P), sys.P, name = run.CH4_MODEL) +
                                        sl.co_melt((CO.Y*mCO*sys.P), sys.P, name = run.CO_MODEL) +
                                        melt.graph_current)
    return sys.atomicM

def soh(sys, run, melt, gas, mols):
    """
    Calculates the initial conditions for the SOH system.

    Based on the subset of parameters provided at input, where the starting
    pressure is known, the speciation of the gas phase as mole fractions is
    calculated. Based on this data, the volatile content of the melt is
    found and the masses of each volatile element (O, C, H etc) is set,
    ready for use in the equilibrium constant and mass-balance solution.

    Parameters
    ----------
    sys : ThermoSystem class
        The active instance of the ThermoSystem class
    run : RunDef class
        The active instance of the RunDef class    
    melt : Melt class
        The active instance of the Melt class
    gas : Gas class
        The active instance of the Gas class
    mols : tuple of Molecule classes
        A tuple of all the Molecule class instances. Setup to be passed
        through all functions, preserving ordering.

    Returns
    -------
    sys.atomicM : dict
        A dictionary containing the mass fraction of each volatile element
        in the system.
    """
    H2O, O2, H2, S2, SO2, H2S = mols

    # uses get_inputs to find which parameters have been set and get the water system.
    mH2O, mO2, mH2 = get_inputs(sys, run, melt, gas, mols)

    assert 1 > (1 - mH2O - mH2 - mO2) > 0, 'Some fugacities are negative!'

    # Quadratic solve for SO2

    def quadratic():

        a = SO2.Y ** 2 / (((sys.K['K5'] * O2.Y * mO2) ** 2) * S2.Y * sys.P)
        b = 1 + ((sys.K['K4'] * H2O.Y * mH2O * SO2.Y) / (sys.K['K5'] * H2S.Y * ((O2.Y * mO2) ** 1.5) * sqrt(sys.P)))
        c = -(1 - mH2 - mO2 - mH2O)

        mSO2 = [x for x in np.roots([a, b, c]) if x > 0 and x < 1]
        return mSO2[0]

    mSO2 = quadratic()
    gas.mSO2.append(mSO2)

    mS2 = (SO2.Y * mSO2) ** 2 / ((sys.K['K5'] * O2.Y * mO2) ** 2 * S2.Y * sys.P)
    gas.mS2.append(mS2)

    mH2S = (1 - mO2 - mS2 - mSO2) / (1 + (H2S.Y / (sys.K['K1'] * sys.K['K4'] * H2.Y * (S2.Y * mS2 * sys.P) ** 0.5)) + ((H2S.Y * (O2.Y * mO2) ** 0.5) / (H2O.Y * sys.K['K4'] * (S2.Y * mS2) ** 0.5)))
    gas.mH2S.append(mH2S)

    empty_lists = [gas.mCO2, gas.mCO, gas.mCH4, gas.mN2]
    for lst in empty_lists:
        lst.append(np.float64(0))

    lst = {'o2':mO2, 'h2':mH2, 'h2o':mH2O, 'so2':mSO2, 's2':mS2, 'h2s':mH2S}
    mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

    # WtO
    sys.atomicM['o'] = cnst.m['o'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + 2*mO2 + 2*mSO2) + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL))

    # WtH
    sys.atomicM['h'] = 2 * cnst.m['h'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + mH2 + mH2S) +
                                        sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL) +
                                        sl.h2_melt(mH2, H2, sys.P, melt, name=run.H2_MODEL))

    # WtS
    sys.atomicM['s'] = cnst.m['s'] * ((sys.WgT[0] / sum(mjMj)) * (2*mS2 + mH2S + mSO2) + 
                                    sl.sulfide_melt((S2.Y*mS2*sys.P), (O2.Y*mO2*sys.P), sys.P, sys.T, melt, name=run.SULFIDE_CAPACITY) +
                                    sl.sulfate_melt((S2.Y*mS2*sys.P), (O2.Y*mO2*sys.P), sys.P, sys.T, melt, run, name=run.SULFATE_CAPACITY))

    return sys.atomicM

def cohs(sys, run, melt, gas, mols):
    """
    Calculates the initial conditions for the COHS system.

    Based on the subset of parameters provided at input, where the starting
    pressure is known, the speciation of the gas phase as mole fractions is
    calculated. Based on this data, the volatile content of the melt is
    found and the masses of each volatile element (O, C, H etc) is set,
    ready for use in the equilibrium constant and mass-balance solution.

    Parameters
    ----------
    sys : ThermoSystem class
        The active instance of the ThermoSystem class
    run : RunDef class
        The active instance of the RunDef class    
    melt : Melt class
        The active instance of the Melt class
    gas : Gas class
        The active instance of the Gas class
    mols : tuple of Molecule classes
        A tuple of all the Molecule class instances. Setup to be passed
        through all functions, preserving ordering.

    Returns
    -------
    sys.atomicM : dict
        A dictionary containing the mass fraction of each volatile element
        in the system.
    """
    H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S = mols

    # uses get_inputs to find which parameters have been set and get the water system + CO2 for full system.
    if not run.SULFUR_SET:
        mCO2, mH2O, mO2, mH2 = get_inputs(sys, run, melt, gas, mols)
        
        assert 1 > (1 - mH2O - mH2 - mO2) > 0, 'Some fugacities are negative!'
        assert 1 > mCO2 > 0, 'CO2 fraction is too high for specified gas fraction'

        mCO = CO2.Y*mCO2/(sys.K['K2']*CO.Y*(O2.Y*mO2*sys.P)**0.5)
        gas.mCO.append(mCO)

        mCH4 = CO2.Y*mCO2*(H2O.Y*mH2O*sys.P)**2/(sys.K['K3']*CH4.Y*(O2.Y*mO2*sys.P)**2)
        gas.mCH4.append(mCH4)

        # Quadratic solve for SO2

        def quadratic():

            a = (SO2.Y / (sys.K['K5'] * O2.Y * mO2)) ** 2 / (S2.Y * sys.P)
            b = 1+(sys.K['K4']*H2O.Y*mH2O*SO2.Y/(sys.K['K5']*H2S.Y*(O2.Y*mO2)**1.5*sys.P**0.5))
            c = -(1 - mH2 - mO2 - mH2O - mCO - mCO2 - mCH4)

            mSO2 = [x for x in np.roots([a, b, c]) if x > 0 and x < 1]
            if mSO2:
                return mSO2[0]
            else:
                exit('fSO2 is negative. Gas fraction or melt volatile content too high.')

        mSO2 = quadratic()
        gas.mSO2.append(mSO2)

        mS2 = (SO2.Y*mSO2)**2 / ((sys.K['K5']*(O2.Y*mO2))**2 * S2.Y*sys.P)
        gas.mS2.append(mS2)

        mH2S = (sys.K['K4'] * H2O.Y * mH2O * SO2.Y * mSO2) / (sys.K['K5'] * H2S.Y * (O2.Y * mO2) ** 1.5 * sys.P ** 0.5)
        gas.mH2S.append(mH2S)
    
    else: # sulfur has been used as an input
        mS2, mH2O, mO2, mH2 = get_inputs(sys, run, melt, gas, mols)

        assert 1 > (1 - mH2O - mH2 - mO2) > 0, 'Some fugacities are negative!'
        assert 1 > mS2 > 0, 'S2 fraction wrong'
    
        mSO2 = (sys.K['K5'] * O2.Y * mO2 * (S2.Y*mS2*sys.P)**0.5)/SO2.Y
        gas.mSO2.append(mSO2)

        mH2S = (sys.K['K4'] * H2O.Y * mH2O * (S2.Y*mS2)**0.5)/((O2.Y*mO2)**0.5 * H2S.Y)
        gas.mH2S.append(mH2S)

        mCH4 = (1 - mH2 - mO2 - mH2O - mS2 - mSO2 - mH2S)/(1 + ((sys.K['K3']*CH4.Y*(O2.Y*mO2)**2)/(CO2.Y*(H2O.Y*mH2O)**2)) + ((sys.K['K3']*CH4.Y*(O2.Y*mO2)**2)/(sys.K['K2']*CO.Y*(O2.Y*mO2*sys.P)**0.5*(H2O.Y*mH2O)**2)))
        gas.mCH4.append(mCH4)

        mCO2 = (sys.K['K3'] * CH4.Y*mCH4*(O2.Y*mO2)**2)/((H2O.Y*mH2O)**2*CO2.Y)
        gas.mCO2.append(mCO2)

        # Check graphite saturation
        if run.C_MODEL == 'eguchi2018':
            graph_fCO2 = sl.graphite_fco2(sys.T, sys.P, (O2.Y*mO2*sys.P))
            if (CO2.Y*mCO2*sys.P) > graph_fCO2:
                exit('Error: Melt is graphite saturated. System is overconstrained; please rerun by setting melt graphite content rather than sulphur, or stop using eguchi2018 in C_MODEL.')

        mCO = CO2.Y*mCO2/(sys.K['K2']*CO.Y*(O2.Y*mO2*sys.P)**0.5)
        gas.mCO.append(mCO)    


    # Empty list
    gas.mN2.append(0.0)
    
    lst = {'o2':mO2, 'h2':mH2, 'h2o':mH2O, 'co':mCO, 'co2':mCO2, 'ch4':mCH4, 'so2':mSO2, 's2':mS2, 'h2s':mH2S}
    mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

    # WtO
    sys.atomicM['o'] = cnst.m['o'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + 2*mO2 + 2*mCO2 + mCO + 2*mSO2)
                                   + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
                                   + 2*sl.co2_melt((CO2.Y*mCO2*sys.P), CO2, (O2.Y*mO2*sys.P), sys.T, sys.P, melt, name=run.C_MODEL)
                                   + sl.co_melt((CO.Y*mCO*sys.P), sys.P, name = run.CO_MODEL))

    # WtH
    sys.atomicM['h'] = 2 * cnst.m['h'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + mH2 + 2*mCH4 + mH2S) +
                                            sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL) +
                                            sl.h2_melt(mH2, H2, sys.P, melt, name=run.H2_MODEL) +
                                            2*sl.ch4_melt((CH4.Y*mCH4*sys.P), sys.P, name = run.CH4_MODEL))

    # WtC
    sys.atomicM['c'] = cnst.m['c'] * ((sys.WgT[0] / sum(mjMj)) * (mCO + mCO2 + mCH4) +
                                            sl.co2_melt((CO2.Y*mCO2*sys.P), CO2, (O2.Y*mO2*sys.P), sys.T, sys.P, melt, name=run.C_MODEL) +
                                            sl.ch4_melt((CH4.Y*mCH4*sys.P), sys.P, name = run.CH4_MODEL) +
                                            sl.co_melt((CO.Y*mCO*sys.P), sys.P, name = run.CO_MODEL) +
                                            melt.graph_current)

    # WtS
    sys.atomicM['s'] = cnst.m['s'] * ((sys.WgT[0] / sum(mjMj)) * (2*mS2 + mH2S + mSO2) + 
                                    sl.sulfide_melt((S2.Y*mS2*sys.P), (O2.Y*mO2*sys.P), sys.P, sys.T, melt, name=run.SULFIDE_CAPACITY) +
                                    sl.sulfate_melt((S2.Y*mS2*sys.P), (O2.Y*mO2*sys.P), sys.P, sys.T, melt, run, name=run.SULFATE_CAPACITY))

    return sys.atomicM

def cohsn(sys, run, melt, gas, mols):
    """
    Calculates the initial conditions for the COHSN system.

    Based on the subset of parameters provided at input, where the starting
    pressure is known, the speciation of the gas phase as mole fractions is
    calculated. Based on this data, the volatile content of the melt is
    found and the masses of each volatile element (O, C, H etc) is set,
    ready for use in the equilibrium constant and mass-balance solution.

    Parameters
    ----------
    sys : ThermoSystem class
        The active instance of the ThermoSystem class
    run : RunDef class
        The active instance of the RunDef class    
    melt : Melt class
        The active instance of the Melt class
    gas : Gas class
        The active instance of the Gas class
    mols : tuple of Molecule classes
        A tuple of all the Molecule class instances. Setup to be passed
        through all functions, preserving ordering.

    Returns
    -------
    sys.atomicM : dict
        A dictionary containing the mass fraction of each volatile element
        in the system.
    """
    H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2 = mols

    # uses get_inputs to find which parameters have been set and get the water system + CO2 for full system.
    if not run.SULFUR_SET:
        mCO2, mH2O, mO2, mH2 = get_inputs(sys, run, melt, gas, mols)
        
        assert 1 > (1 - mH2O - mH2 - mO2) > 0, 'Some fugacities are negative!'
        assert 1 > mCO2 > 0, 'CO2 fraction is too high for specified gas fraction'

        mCO = CO2.Y*mCO2/(sys.K['K2']*CO.Y*(O2.Y*mO2*sys.P)**0.5)
        gas.mCO.append(mCO)

        mCH4 = CO2.Y*mCO2*(H2O.Y*mH2O*sys.P)**2/(sys.K['K3']*CH4.Y*(O2.Y*mO2*sys.P)**2)
        gas.mCH4.append(mCH4)

        mN2 = sl.n2_fugacity(run.NITROGEN_START, N2.Y, (O2.Y * mO2 * sys.P), sys.P, name=run.N_MODEL)/(N2.Y*sys.P)
        
        if mN2 > 1:
            exit('Warning: Too high a nitrogen content for the pressure.')
        gas.mN2.append(mN2)

        # Quadratic solve for SO2

        def quadratic():

            a = (SO2.Y / (sys.K['K5'] * O2.Y * mO2)) ** 2 / (S2.Y * sys.P)
            b = 1+(sys.K['K4']*H2O.Y*mH2O*SO2.Y/(sys.K['K5']*H2S.Y*(O2.Y*mO2)**1.5*sys.P**0.5))
            c = -(1 - mH2 - mO2 - mH2O - mCO - mCO2 - mCH4 - mN2)

            mSO2 = [x for x in np.roots([a, b, c]) if x > 0 and x < 1]
            if mSO2:
                return mSO2[0]
            else:
                exit('fSO2 is negative. Gas fraction or melt volatile content too high.')

        mSO2 = quadratic()
        gas.mSO2.append(mSO2)

        mS2 = (SO2.Y*mSO2)**2 / ((sys.K['K5']*(O2.Y*mO2))**2 * S2.Y*sys.P)
        gas.mS2.append(mS2)

        mH2S = (sys.K['K4'] * H2O.Y * mH2O * SO2.Y * mSO2) / (sys.K['K5'] * H2S.Y * (O2.Y * mO2) ** 1.5 * sys.P ** 0.5)
        gas.mH2S.append(mH2S)
    
    else: # sulfur has been used as an input
        mS2, mH2O, mO2, mH2 = get_inputs(sys, run, melt, gas, mols)

        assert 1 > (1 - mH2O - mH2 - mO2) > 0, 'Some fugacities are negative!'
        assert 1 > mS2 > 0, 'S2 fraction wrong'
    
        mSO2 = (sys.K['K5'] * O2.Y * mO2 * (S2.Y*mS2*sys.P)**0.5)/SO2.Y
        gas.mSO2.append(mSO2)

        mH2S = (sys.K['K4'] * H2O.Y * mH2O * (S2.Y*mS2)**0.5)/((O2.Y*mO2)**0.5 * H2S.Y)
        gas.mH2S.append(mH2S)

        mN2 = sl.n2_fugacity(run.NITROGEN_START, N2.Y, (O2.Y * mO2 * sys.P), sys.P, name=run.N_MODEL)/(N2.Y*sys.P)
        
        if mN2 > 1:
            exit('Warning: Too high a nitrogen content for the pressure.')
        gas.mN2.append(mN2)

        mCH4 = (1 - mH2 - mO2 - mH2O - mS2 - mSO2 - mH2S - mN2)/(1 + ((sys.K['K3']*CH4.Y*(O2.Y*mO2)**2)/(CO2.Y*(H2O.Y*mH2O)**2)) + ((sys.K['K3']*CH4.Y*(O2.Y*mO2)**2)/(sys.K['K2']*CO.Y*(O2.Y*mO2*sys.P)**0.5*(H2O.Y*mH2O)**2)))
        gas.mCH4.append(mCH4)

        mCO2 = (sys.K['K3'] * CH4.Y*mCH4*(O2.Y*mO2)**2)/((H2O.Y*mH2O)**2*CO2.Y)
        gas.mCO2.append(mCO2)

        # Check graphite saturation
        if run.C_MODEL == 'eguchi2018':
            graph_fCO2 = sl.graphite_fco2(sys.T, sys.P, (O2.Y*mO2*sys.P))
            if (CO2.Y*mCO2*sys.P) > graph_fCO2:
                exit('Error: Melt is graphite saturated. System is overconstrained; please rerun by setting melt graphite content rather than sulphur, or stop using eguchi2018 in C_MODEL.')
                    
        mCO = CO2.Y*mCO2/(sys.K['K2']*CO.Y*(O2.Y*mO2*sys.P)**0.5)
        gas.mCO.append(mCO)    


    lst = {'o2':mO2, 'h2':mH2, 'h2o':mH2O, 'co':mCO, 'co2':mCO2, 'ch4':mCH4, 'so2':mSO2, 's2':mS2, 'h2s':mH2S, 'n2':mN2}
    mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

    # WtO
    sys.atomicM['o'] = cnst.m['o'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + 2*mO2 + 2*mCO2 + mCO + 2*mSO2)
                                            + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
                                        + 2*sl.co2_melt((CO2.Y*mCO2*sys.P), CO2, (O2.Y*mO2*sys.P), sys.T, sys.P, melt, name=run.C_MODEL) +
                                        sl.co_melt((CO.Y*mCO*sys.P), sys.P, name = run.CO_MODEL))

    # WtH
    sys.atomicM['h'] = 2 * cnst.m['h'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + mH2 + 2*mCH4 + mH2S) +
                                            sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL) +
                                            sl.h2_melt(mH2, H2, sys.P, melt, name=run.H2_MODEL) +
                                            2*sl.ch4_melt((CH4.Y*mCH4*sys.P), sys.P, name = run.CH4_MODEL))

    # WtC
    sys.atomicM['c'] = cnst.m['c'] * ((sys.WgT[0] / sum(mjMj)) * (mCO + mCO2 + mCH4) +
                                            sl.co2_melt((CO2.Y*mCO2*sys.P), CO2, (O2.Y*mO2*sys.P), sys.T, sys.P, melt, name=run.C_MODEL) +
                                            sl.ch4_melt((CH4.Y*mCH4*sys.P), sys.P, name = run.CH4_MODEL) +
                                            sl.co_melt((CO.Y*mCO*sys.P), sys.P, name = run.CO_MODEL) + 
                                            melt.graph_current)

    # WtS
    sys.atomicM['s'] = cnst.m['s'] * ((sys.WgT[0] / sum(mjMj)) * (2*mS2 + mH2S + mSO2) + 
                                    sl.sulfide_melt((S2.Y*mS2*sys.P), (O2.Y*mO2*sys.P), sys.P, sys.T, melt, name=run.SULFIDE_CAPACITY) +
                                    sl.sulfate_melt((S2.Y*mS2*sys.P), (O2.Y*mO2*sys.P), sys.P, sys.T, melt, run, name=run.SULFATE_CAPACITY))
    
    # WtN
    sys.atomicM['n'] = cnst.m['n'] * ((sys.WgT[0] / sum(mjMj)) * (2*mN2) + sl.n_melt(mN2, (O2.Y*mO2*sys.P), sys.P, name=run.N_MODEL))
    return sys.atomicM