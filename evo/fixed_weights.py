"""
Code to find the saturation point of a melt based on it's volatile content given
as atomic weights of species in the melt. so not quite total atomic mass in the system..

Sets WgT to 1e-6 wt% and solves to find the speciation of both the gas phase and
the speciation within the melt.
"""

import numpy as np
from scipy.optimize import fsolve, root
import warnings

from evo import constants as cnst
from evo import conversions as cnvs
from evo import messages as msgs
from evo import solubility_laws as sl
from evo import solvgas as sg

def get_molfrac(P, fugacities, gamma):
    "Converts fugacities into mol fractions using pressure and fugacity coefficients."

    fh2o, fo2, fh2, fco, fco2, fch4, fs2, fso2, fh2s, fn2 = fugacities
    h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, n2_y = gamma
    
    if fco2 == 0 and fs2 == 0:
        return fh2o/(h2o_y*P), fo2/(o2_y*P), fh2/(h2_y*P)
    elif fco2 == 0:
        return fh2o/(h2o_y*P), fo2/(o2_y*P), fh2/(h2_y*P), fs2/(s2_y*P), fso2/(so2_y*P), fh2s/(h2s_y*P)
    elif fs2 == 0:
        return fh2o/(h2o_y*P), fo2/(o2_y*P), fh2/(h2_y*P), fco/(co_y*P), fco2/(co2_y*P), fch4/(ch4_y*P)
    elif fn2 == 0:
        return fh2o/(h2o_y*P), fo2/(o2_y*P), fh2/(h2_y*P), fco/(co_y*P), fco2/(co2_y*P), fch4/(ch4_y*P), fs2/(s2_y*P), fso2/(so2_y*P), fh2s/(h2s_y*P)
    else:
        return fh2o/(h2o_y*P), fo2/(o2_y*P), fh2/(h2_y*P), fco/(co_y*P), fco2/(co2_y*P), fch4/(ch4_y*P), fs2/(s2_y*P), fso2/(so2_y*P), fh2s/(h2s_y*P), fn2/(n2_y*P)

def p_tot(P, sys, fugacities, gamma):
    "Returns the sum of the gas phase partial pressures, given the current system pressure"

    fh2o, fo2, fh2, fco, fco2, fch4, fs2, fso2, fh2s, fn2 = fugacities
    h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, n2_y = gamma

    if fco2 == 0 and fs2 == 0:
        return fh2o/h2o_y + fo2/o2_y + fh2/h2_y - P
    elif fco2 == 0:
        return fh2o/h2o_y + fo2/o2_y + fh2/h2_y + fs2/s2_y + fso2/so2_y + fh2s/h2s_y - P
    elif fs2 == 0:
        return fh2o/h2o_y + fo2/o2_y + fh2/h2_y + fco/co_y + fco2/co2_y + fch4/ch4_y - P
    elif fn2 == 0:
        return fh2o/h2o_y + fo2/o2_y + fh2/h2_y + fco/co_y + fco2/co2_y + fch4/ch4_y + fs2/s2_y + fso2/so2_y + fh2s/h2s_y - P
    else:
        return fh2o/h2o_y + fo2/o2_y + fh2/h2_y + fco/co_y + fco2/co2_y + fch4/ch4_y + fs2/s2_y + fso2/so2_y + fh2s/h2s_y + fn2/n2_y - P


def sat_pressure(run, sys, gas, melt, mols):
    """
    Calculate the volatile saturation pressure of a run, set the atomic masses at that pressure then round down to the nearest 1 bar to set the starting pressure
    for the solver loop.
    """

    if run.GAS_SYS == 'OH':
        H2O, O2, H2 = mols
        sys.atomicM['h'] = run.ATOMIC_H/1e6
    
    elif run.GAS_SYS == 'COH':
        H2O, O2, H2, CO, CO2, CH4 = mols
        sys.atomicM['c'] = run.ATOMIC_C/1e6
        sys.atomicM['h'] = run.ATOMIC_H/1e6

    elif run.GAS_SYS == 'SOH':
        H2O, O2, H2, S2, SO2, H2S = mols
        sys.atomicM['h'] = run.ATOMIC_H/1e6
        sys.atomicM['s'] = run.ATOMIC_S/1e6

    elif run.GAS_SYS == "COHS":
        H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S = mols
        sys.atomicM['c'] = run.ATOMIC_C/1e6
        sys.atomicM['h'] = run.ATOMIC_H/1e6
        sys.atomicM['s'] = run.ATOMIC_S/1e6
    
    elif run.GAS_SYS == "COHSN":
        H2O, O2, H2, CO, CO2, CH4, S2, SO2, H2S, N2 = mols
        sys.atomicM['c'] = run.ATOMIC_C/1e6
        sys.atomicM['h'] = run.ATOMIC_H/1e6
        sys.atomicM['n'] = run.ATOMIC_N/1e6
        sys.atomicM['s'] = run.ATOMIC_S/1e6

    def get_f(P, melt_h2o, melt_co2, melt_s, melt_n, sys, melt, gamma):
        """
        Calculate the gas fugacities of species based on the melt composition and total pressure.
        """

        O2.Y = gamma[1]
        
        if run.FO2_buffer_SET:
            fo2 = np.exp(cnvs.generate_fo2(sys, sys.run.FO2_buffer_START, sys.run.FO2_buffer, P))
        
        else:
            fo2 = cnvs.c2fo2(melt.Cm(), sys.T, P*1e5, sys.run.FO2_MODEL)
        
        fh2o = sl.h2o_fugacity(melt_h2o, H2O, name=run.H2O_MODEL)
        fh2 = fh2o / (sys.K['K1']*fo2**0.5)

        if melt_co2 != 0.0:
            fco2 = sl.co2_fugacity(melt_co2, CO2, fo2, sys.T, P, melt, name=run.C_MODEL)
            fco = fco2/(sys.K['K2']*fo2**0.5)
            fch4 = (fco2 * fh2o**2) /(sys.K['K3']*fo2**2)
        else:
            fco2 = 0
            fco = 0
            fch4 = 0

        if melt_s != 0.0:
            fs2 = sl.S2_fugacity(melt_s, fo2, P, sys.T, melt, run, sulfidename=run.SULFIDE_CAPACITY, sulfatename=run.SULFATE_CAPACITY)
            fso2 = sys.K['K5'] * fo2 * fs2**0.5
            fh2s = (sys.K['K4']*fh2o*fs2**0.5)/(fo2**0.5)
        
        else:
            fs2 = 0
            fso2 = 0
            fh2s = 0

        if melt_n!= 0.0:
            n2_y = gamma[-1]
            fn2 = sl.n2_fugacity(melt_n, n2_y, fo2, P, name=run.N_MODEL)
        else:
            fn2 = 0          
        
        return fh2o, fo2, fh2, fco, fco2, fch4, fs2, fso2, fh2s, fn2
    
    def find_p(P, melt_h2o, melt_co2, melt_s, melt_n, sys, melt):

        gamma = sg.find_Y(P, sys.T, sys.SC)[:10]
        fugacity = get_f(P, melt_h2o, melt_co2, melt_s, melt_n, sys, melt, gamma)
        
        return p_tot(P, sys, fugacity, gamma)

    def fixed_weights_oh(guesses=[0.0]):
        melt_h2o = guesses[0]
        melt_co2, melt_n, melt_s = 0.0, 0.0, 0.0

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                P_sat = fsolve(find_p, 1.0, args=(melt_h2o, melt_co2, melt_s, melt_n, sys, melt))
        except:
            raise RuntimeError('Failed to find saturation pressure; solver not converging.')
       
        # Check p_sat hasn't returned 1.0
        if P_sat == 1.0:
            exit('Saturation solver has failed :(. Exiting...')
        
        P_sat = P_sat[0]
        sys.P = P_sat
        
        gamma = sg.find_Y(P_sat, sys.T, sys.SC)[:10]
        fugacities = get_f(P_sat, melt_h2o, melt_co2, melt_s, melt_n, sys, melt, gamma)

        h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, n2_y = gamma
        O2.Y = o2_y
        
        mH2O, mO2, mH2 = get_molfrac(P_sat, fugacities, gamma)
        
        lst = {'o2': mO2, 'h2': mH2, 'h2o': mH2O}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]
        N = 1e-8/sum(mjMj)

        return [(N * (mH2O + mH2) + sl.h2_melt(mH2, H2, P_sat, melt, name=run.H2_MODEL, Y=h2_y) + sl.h2o_melt(mH2O, H2O, P_sat, name=run.H2O_MODEL, Y=h2o_y)) - (sys.atomicM['h']/(2*cnst.m['h']))]
    
    def fixed_weights_coh(guesses=[0.0, 0.0]):
        melt_h2o, melt_co2 = guesses[0], guesses[1]
        melt_n, melt_s = 0.0, 0.0

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                P_sat = fsolve(find_p, 1.0, args=(melt_h2o, melt_co2, melt_s, melt_n, sys, melt))
        except:
            raise RuntimeError('Failed to find saturation pressure; solver not converging.')
       
        # Check p_sat hasn't returned 1.0
        if P_sat == 1.0:
            exit('Saturation solver has failed :(. Exiting...')
        
        P_sat = P_sat[0]
        sys.P = P_sat
        
        gamma = sg.find_Y(P_sat, sys.T, sys.SC)[:10]
        fugacities = get_f(P_sat, melt_h2o, melt_co2, melt_s, melt_n, sys, melt, gamma)

        h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, n2_y = gamma
        O2.Y = o2_y
        
        mH2O, mO2, mH2, mCO, mCO2, mCH4 = get_molfrac(P_sat, fugacities, gamma)
        
        lst = {'o2': mO2, 'h2': mH2, 'h2o': mH2O, 'co': mCO, 'co2': mCO2, 'ch4': mCH4}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]
        N = 1e-8/sum(mjMj)

        if melt.graphite_sat == True:                      
            melt.graph_current = cnvs.get_graphite(sys, melt, P_sat, CO2, mCO, mCO2, mCH4, mO2, co2_y, co_y, ch4_y, o2_y, N)/cnst.m['c']

        elif melt.graphite_sat == False:
            melt.graph_current = 0

        return [(N * (mH2O + mH2 + 2*mCH4) + sl.h2_melt(mH2, H2, P_sat, melt, name=run.H2_MODEL, Y=h2_y) + sl.h2o_melt(mH2O, H2O, P_sat, name=run.H2O_MODEL, Y=h2o_y) + 2*sl.ch4_melt((ch4_y*mCH4*P_sat), P_sat, name = run.CH4_MODEL)) - (sys.atomicM['h']/(2*cnst.m['h'])),
        
        (N * (mCO + mCO2 + mCH4) + sl.co2_melt((co2_y*mCO2*P_sat), CO2, (o2_y*mO2*P_sat), sys.T, P_sat, melt, name=run.C_MODEL) + sl.co_melt((co_y*mCO*P_sat), P_sat, name = run.CO_MODEL) + sl.ch4_melt((ch4_y*mCH4*P_sat), P_sat, name = run.CH4_MODEL) + melt.graph_current) - (sys.atomicM['c']/cnst.m['c'])]
    
    def fixed_weights_soh(guesses=[0.0, 0.0]):
        melt_h2o, melt_s = guesses[0], guesses[1]
        melt_n, melt_co2 = 0.0, 0.0

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                P_sat = fsolve(find_p, 1.0, args=(melt_h2o, melt_co2, melt_s, melt_n, sys, melt))
        except:
            raise RuntimeError('Failed to find saturation pressure; solver not converging.')
       
        # Check p_sat hasn't returned 1.0
        if P_sat == 1.0:
            exit('Saturation solver has failed :(. Exiting...')
        
        P_sat = P_sat[0]
        sys.P = P_sat
        
        gamma = sg.find_Y(P_sat, sys.T, sys.SC)[:10]
        fugacities = get_f(P_sat, melt_h2o, melt_co2, melt_s, melt_n, sys, melt, gamma)

        h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, n2_y = gamma
        O2.Y = o2_y
        
        mH2O, mO2, mH2, mS2, mSO2, mH2S = get_molfrac(P_sat, fugacities, gamma)
        
        lst = {'o2': mO2, 'h2': mH2, 'h2o': mH2O, 's2': mS2, 'so2': mSO2, 'h2s': mH2S}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]
        N = 1e-8/sum(mjMj)

        return [(N * (mH2O + mH2 + mH2S) + sl.h2_melt(mH2, H2, P_sat, melt, name=run.H2_MODEL, Y=h2_y) + sl.h2o_melt(mH2O, H2O, P_sat, name=run.H2O_MODEL, Y=h2o_y)) - (sys.atomicM['h']/(2*cnst.m['h'])),
        
        (N*(mSO2 + mH2S + 2*mS2) + sl.sulfide_melt((s2_y*mS2*P_sat), (o2_y*mO2*P_sat), P_sat, sys.T, melt, name=run.SULFIDE_CAPACITY) + sl.sulfate_melt((s2_y*mS2*P_sat), (o2_y*mO2*P_sat), P_sat, sys.T, melt, run, name=run.SULFATE_CAPACITY)) - (sys.atomicM['s']/cnst.m['s'])]
    
    def fixed_weights_cohs(guesses=[0.0, 0.0, 0.0]):
        melt_h2o, melt_co2, melt_s = guesses[0], guesses[1], guesses[2]
        melt_n = 0.0 

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                P_sat = fsolve(find_p, 1.0, args=(melt_h2o, melt_co2, melt_s, melt_n, sys, melt))
        except:
            raise RuntimeError('Failed to find saturation pressure; solver not converging.')
       
        # Check p_sat hasn't returned 1.0
        if P_sat == 1.0:
            exit('Saturation solver has failed :(. Exiting...')
        
        P_sat = P_sat[0]
        sys.P = P_sat
        
        gamma = sg.find_Y(P_sat, sys.T, sys.SC)[:10]
        fugacities = get_f(P_sat, melt_h2o, melt_co2, melt_s, melt_n, sys, melt, gamma)

        h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, n2_y = gamma
        O2.Y = gamma[1]

        mH2O, mO2, mH2, mCO, mCO2, mCH4, mS2, mSO2, mH2S = get_molfrac(P_sat, fugacities, gamma)
        
        lst = {'o2': mO2, 'h2': mH2, 'h2o': mH2O, 'co': mCO, 'co2': mCO2, 'ch4': mCH4, 's2': mS2, 'so2': mSO2, 'h2s': mH2S}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]
        N = 1e-8/sum(mjMj)

        if melt.graphite_sat == True:                      
            melt.graph_current = cnvs.get_graphite(sys, melt, P_sat, CO2, mCO, mCO2, mCH4, mO2, co2_y, co_y, ch4_y, o2_y, N)/cnst.m['c']

        elif melt.graphite_sat == False:
            melt.graph_current = 0

        return [(N * (mH2O + mH2 + mH2S + 2*mCH4) + sl.h2_melt(mH2, H2, P_sat, melt, name=run.H2_MODEL, Y=h2_y) + sl.h2o_melt(mH2O, H2O, P_sat, name=run.H2O_MODEL, Y=h2o_y) + 2*sl.ch4_melt((ch4_y*mCH4*P_sat), P_sat, name = run.CH4_MODEL)) - (sys.atomicM['h']/(2*cnst.m['h'])),

        (N*(mSO2 + mH2S + 2*mS2) + sl.sulfide_melt((s2_y*mS2*P_sat), (o2_y*mO2*P_sat), P_sat, sys.T, melt, name=run.SULFIDE_CAPACITY) + sl.sulfate_melt((s2_y*mS2*P_sat), (o2_y*mO2*P_sat), P_sat, sys.T, melt, run, name=run.SULFATE_CAPACITY)) - (sys.atomicM['s']/cnst.m['s']),
        
        (N * (mCO + mCO2 + mCH4) + sl.co2_melt((co2_y*mCO2*P_sat), CO2, (o2_y*mO2*P_sat), sys.T, P_sat, melt, name=run.C_MODEL) + sl.co_melt((co_y*mCO*P_sat), P_sat, name = run.CO_MODEL) + sl.ch4_melt((ch4_y*mCH4*P_sat), P_sat, name = run.CH4_MODEL) + melt.graph_current) - (sys.atomicM['c']/cnst.m['c'])]

    def fixed_weights_cohsn(guesses=[0.0, 0.0, 0.0, 0.0]):
        melt_h2o, melt_co2, melt_s, melt_n = guesses[0], guesses[1], guesses[2], guesses[3]      

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                P_sat = fsolve(find_p, 1.0, args=(melt_h2o, melt_co2, melt_s, melt_n, sys, melt))
        except:
            raise RuntimeError('Failed to find saturation pressure; solver not converging.')
       
        # Check p_sat hasn't returned 1.0
        if P_sat == 1.0:
            exit('Saturation solver has failed :(. Exiting...')
        
        P_sat = P_sat[0]
        sys.P = P_sat
        
        gamma = sg.find_Y(P_sat, sys.T, sys.SC)[:10]
        fugacities = get_f(P_sat, melt_h2o, melt_co2, melt_s, melt_n, sys, melt, gamma)

        h2o_y, o2_y, h2_y, co_y, co2_y, ch4_y, s2_y, so2_y, h2s_y, n2_y = gamma
        O2.Y = o2_y

        mH2O, mO2, mH2, mCO, mCO2, mCH4, mS2, mSO2, mH2S, mN2 = get_molfrac(P_sat, fugacities, gamma)

        lst = {'o2': mO2, 'h2': mH2, 'h2o': mH2O, 'co': mCO, 'co2': mCO2, 'ch4': mCH4, 's2': mS2, 'so2': mSO2, 'h2s': mH2S, 'n2':mN2}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]
        N = 1e-8/sum(mjMj)

        if melt.graphite_sat == True:                      
            melt.graph_current = cnvs.get_graphite(sys, melt, P_sat, CO2, mCO, mCO2, mCH4, mO2, co2_y, co_y, ch4_y, o2_y, N)/cnst.m['c']

        elif melt.graphite_sat == False:
            melt.graph_current = 0

        return [(N*(2*mN2) + sl.n_melt(mN2, (o2_y*mO2*P_sat), P_sat, name=run.N_MODEL)) - (sys.atomicM['n']/cnst.m['n']),

        (N * (mH2O + mH2 + mH2S + 2*mCH4) + sl.h2_melt(mH2, H2, P_sat, melt, name=run.H2_MODEL, Y=h2_y) + sl.h2o_melt(mH2O, H2O, P_sat, name=run.H2O_MODEL, Y=h2o_y) + 2*sl.ch4_melt((ch4_y*mCH4*P_sat), P_sat, name = run.CH4_MODEL)) - (sys.atomicM['h']/(2*cnst.m['h'])),

        (N*(mSO2 + mH2S + 2*mS2) + sl.sulfide_melt((s2_y*mS2*P_sat), (o2_y*mO2*P_sat), P_sat, sys.T, melt, name=run.SULFIDE_CAPACITY) + sl.sulfate_melt((s2_y*mS2*P_sat), (o2_y*mO2*P_sat), P_sat, sys.T, melt, run, name=run.SULFATE_CAPACITY)) - (sys.atomicM['s']/cnst.m['s']),
        
        (N * (mCO + mCO2 + mCH4) + sl.co2_melt((co2_y*mCO2*P_sat), CO2, (o2_y*mO2*P_sat), sys.T, P_sat, melt, name=run.C_MODEL) + sl.co_melt((co_y*mCO*P_sat), P_sat, name = run.CO_MODEL) + sl.ch4_melt((ch4_y*mCH4*P_sat), P_sat, name = run.CH4_MODEL) + melt.graph_current) - (sys.atomicM['c']/cnst.m['c'])]
        
        ## start of function ------------------------------------------------------------------

    if run.FO2_buffer_SET == True and cnvs.generate_fo2_buffer(sys, np.exp(cnvs.generate_fo2(sys, run.FO2_buffer_START, run.FO2_buffer, 300)), 300, buffer_choice = 'IW') < -1.5:
        guess_h2o = ((run.ATOMIC_H/1e6)/cnst.m['h'])*0.2*cnst.m['h2o']
    else:
        guess_h2o = (((run.ATOMIC_H/1e6)/cnst.m['h'])*cnst.m['h2o'])/2
    
    if run.GAS_SYS == 'COH' or run.GAS_SYS == 'COHS' or run.GAS_SYS == 'COHSN':
        if run.C_MODEL == 'eguchi2018' and cnvs.generate_fo2_buffer(sys, np.exp(cnvs.generate_fo2(sys, run.FO2_buffer_START, run.FO2_buffer, 300)), 300, buffer_choice = 'IW') < -1.5:
            guess_co2 = 1e-5
            melt.graphite_sat == True
        elif run.CO_MODEL != 'None' and run.CH4_MODEL != 'None':
            guess_co2 = (((run.ATOMIC_C/1e6)/cnst.m['c'])*cnst.m['co2']*0.2)
        else:
            guess_co2 = (((run.ATOMIC_C/1e6)/cnst.m['c'])*cnst.m['co2'])
    else:
        melt_co2 = 0
    
    if run.GAS_SYS == 'SOH' or run.GAS_SYS == 'COHS' or run.GAS_SYS == 'COHSN': 
        guess_s = run.ATOMIC_S/1e6
    else:
        melt_s = 0
    
    if run.GAS_SYS == 'COHSN':
        guess_n = run.ATOMIC_N/1e6
    else:
        melt_n = 0
    
    try:
        if run.GAS_SYS == 'OH':
            melt_h2o = fsolve(fixed_weights_oh, [guess_h2o])[0]     # otherwise melt_h2o is an array!
        elif run.GAS_SYS == 'COH':
            melt_h2o, melt_co2 = fsolve(fixed_weights_coh, [guess_h2o, guess_co2])
        elif run.GAS_SYS == 'SOH':
            melt_h2o, melt_s = fsolve(fixed_weights_soh, [guess_h2o, guess_s])
        elif run.GAS_SYS == 'COHS':
            melt_h2o, melt_co2, melt_s = fsolve(fixed_weights_cohs, [guess_h2o, guess_co2, guess_s])
        elif run.GAS_SYS == 'COHSN':
            try:
                melt_h2o, melt_co2, melt_s, melt_n = fsolve(fixed_weights_cohsn, [guess_h2o, guess_co2, guess_s, guess_n])
            except:
                print('Convergence struggling; trying Levenberg-Marquardt method')
                sol = root(fixed_weights_cohsn, [guess_h2o, guess_co2, guess_s, guess_n], method='lm')
                melt_h2o, melt_co2, melt_s, melt_n = sol['x']
            
    except Exception:
        raise RuntimeError('Failed to find saturation pressure; solver not converging.')
    
    P_sat = sys.P
    gamma = sg.find_Y(sys.P, sys.T, sys.SC)[:10]
    fugacities = get_f(sys.P, melt_h2o, melt_co2, melt_s, melt_n, sys, melt, gamma)

    # Reset the system fO2 to match that at the saturation pressure
    if run.FO2_buffer_SET:
        sys.FO2 = cnvs.generate_fo2(sys, sys.run.FO2_buffer_START, sys.run.FO2_buffer, sys.P)
    
    else:
        sys.FO2 = np.log(cnvs.c2fo2(melt.Cm(), sys.T, sys.Ppa, sys.run.FO2_MODEL))
    
    # Updates the fugacity coefficients stored in the molecule classes
    sg.set_Y(sys, mols)
    
    # set the gas fraction at saturation to 1e-6 wt%
    sys.WgT[0] = 1e-8

    # Find the atomic masses at saturation using method which enforces sum(mol fracs) = 1

    if run.GAS_SYS == 'OH':

        mH2O, mO2, mH2 = get_molfrac(sys.P, fugacities, gamma)

        lst = {'o2':mO2, 'h2':mH2, 'h2o':mH2O}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

        # WtO
        sys.atomicM['o'] = cnst.m['o'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + 2*mO2)
                                            + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL))

        sys.atomicM['c'], sys.atomicM['s'], sys.atomicM['n'] = 0.0, 0.0, 0.0
        
        # add gas speciation to lists to act as initial guess for first 'proper' pressure step
        lists = [gas.mH2O, gas.mO2, gas.mH2]
        empty_list = [gas.mCO, gas.mCO2, gas.mCH4, gas.mS2, gas.mSO2, gas.mH2S, gas.mN2]
        values = [mH2O, mO2, mH2]

    elif run.GAS_SYS == 'COH':
        
        mH2O, mO2, mH2, mCO, mCO2, mCH4 = get_molfrac(sys.P, fugacities, gamma)

        lst = {'o2':mO2, 'h2':mH2, 'h2o':mH2O, 'co':mCO, 'co2':mCO2, 'ch4':mCH4}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

        # WtO
        sys.atomicM['o'] = cnst.m['o'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + 2*mO2 + 2*mCO2 + mCO)
                                         + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
                                        + 2*sl.co2_melt((CO2.Y*mCO2*sys.P), CO2, (O2.Y*mO2*sys.P), sys.T, sys.P, melt, name=run.C_MODEL)
                                        + sl.co_melt((CO.Y*mCO*sys.P), sys.P, name = run.CO_MODEL))
        
        sys.atomicM['s'] = 0.0

        sys.atomicM['n'] = 0.0
        
        lists = [gas.mH2O, gas.mO2, gas.mH2, gas.mCO, gas.mCO2, gas.mCH4]
        empty_list = [gas.mS2, gas.mSO2, gas.mH2S, gas.mN2]
        values = [mH2O, mO2, mH2, mCO, mCO2, mCH4]

    elif run.GAS_SYS == 'SOH':

        mH2O, mO2, mH2, mS2, mSO2, mH2S = get_molfrac(sys.P, fugacities, gamma)

        lst = {'o2':mO2, 'h2':mH2, 'h2o':mH2O, 'so2':mSO2, 's2':mS2, 'h2s':mH2S}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

        # WtO
        sys.atomicM['o'] = cnst.m['o'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + 2*mO2 + 2*mSO2) + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL))

        sys.atomicM['c'], sys.atomicM['n'] = 0.0, 0.0
        
        lists = [gas.mH2O, gas.mO2, gas.mH2, gas.mS2, gas.mSO2, gas.mH2S]
        empty_list = [gas.mCO, gas.mCO2, gas.mCH4, gas.mN2]
        values = [mH2O, mO2, mH2, mS2, mSO2, mH2S]

    elif run.GAS_SYS == 'COHS':

        mH2O, mO2, mH2, mCO, mCO2, mCH4, mS2, mSO2, mH2S = get_molfrac(sys.P, fugacities, gamma)

        lst = {'o2':mO2, 'h2':mH2, 'h2o':mH2O, 'co':mCO, 'co2':mCO2, 'ch4':mCH4, 'so2':mSO2, 's2':mS2, 'h2s':mH2S}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

        # WtO
        sys.atomicM['o'] = cnst.m['o'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + 2*mO2 + 2*mCO2 + mCO + 2*mSO2)
                                                + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
                                            + 2*sl.co2_melt((CO2.Y*mCO2*sys.P), CO2, (O2.Y*mO2*sys.P), sys.T, sys.P, melt, name=run.C_MODEL)
                                            + sl.co_melt((CO.Y*mCO*sys.P), sys.P, name = run.CO_MODEL))
        
        sys.atomicM['n'] = 0.0

        lists = [gas.mH2O, gas.mO2, gas.mH2, gas.mCO, gas.mCO2, gas.mCH4, gas.mS2, gas.mSO2, gas.mH2S]
        empty_list = [gas.mN2]
        values = [mH2O, mO2, mH2, mCO, mCO2, mCH4, mS2, mSO2, mH2S]
    
    elif run.GAS_SYS == 'COHSN':

        mH2O, mO2, mH2, mCO, mCO2, mCH4, mS2, mSO2, mH2S, mN2 = get_molfrac(sys.P, fugacities, gamma)

        lst = {'o2':mO2, 'h2':mH2, 'h2o':mH2O, 'co':mCO, 'co2':mCO2, 'ch4':mCH4, 'so2':mSO2, 's2':mS2, 'h2s':mH2S, 'n2':mN2}
        mjMj = [lst[ele] * cnst.m[ele] for ele in lst]

        # WtO
        sys.atomicM['o'] = cnst.m['o'] * ((sys.WgT[0] / sum(mjMj)) * (mH2O + 2*mO2 + 2*mCO2 + mCO + 2*mSO2)
                                                + sl.h2o_melt(mH2O, H2O, sys.P, name=run.H2O_MODEL)
                                            + 2*sl.co2_melt((CO2.Y*mCO2*sys.P), CO2, (O2.Y*mO2*sys.P), sys.T, sys.P, melt, name=run.C_MODEL)
                                            + sl.co_melt((CO.Y*mCO*sys.P), sys.P, name = run.CO_MODEL))

        lists = [gas.mH2O, gas.mO2, gas.mH2, gas.mCO, gas.mCO2, gas.mCH4, gas.mS2, gas.mSO2, gas.mH2S, gas.mN2]
        empty_list = []
        values = [mH2O, mO2, mH2, mCO, mCO2, mCH4, mS2, mSO2, mH2S, mN2]

    # add gas speciation to lists to act as initial guess for first 'proper' pressure step
    for ls, val in zip(lists, values):
        ls.append(val)
    
    for ls in empty_list:
        ls.append(np.float64(0))
    
    # Round P down to nearest bar ready for solver to find WgT at 'starting' pressure.
    sys.P = cnvs.truncate(sys.P, decimals=0)
    run.P_START = sys.P

    return P_sat, values, gamma, mols, melt.graphite_sat
        
