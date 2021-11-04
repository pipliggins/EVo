# conversions.py

import ferric
import constants as cnst
import solubility_laws as sl
import numpy as np

# ------------------------------------------------------------------------
# FUNCTIONS
# ------------------------------------------------------------------------

def K2C(T):
    'Convert degrees K into degrees C'
    return T - 273.15

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return np.floor(n * multiplier) / multiplier

# ------------------------------------------------------------------------
# FOR EQUILIBRATION OF FO2 TO THE FE CONTENT OF THE MAGMA
# ------------------------------------------------------------------------


def norm(c, f=1., cons=None):
    # sum to f
    # if cons are present, this takes them out and  normalises while holding cons constant.
    
    sm = 0
    tmp = {}
    constant = 0

    if cons:
        for x in c:
            if x in cons:
                constant += c[x]
            else:
                sm += c[x]
    else:
        for x in c:
            sm += c[x]

    if cons:
        for x in c:
            if x in cons:
                tmp[x] = c[x]
            else:
                tmp[x] = (f-constant) * c[x] / sm
    else:
        for x in c:
            tmp[x] = f * c[x] / sm

    return tmp


def wt2mol(c, *args):  # wt% of a system into molar mass (EQ 16 in supplementary material)
    "c is a list of weight fractions we want to convert. Any args specify if we want to return only the specified values."
    
    mol = {}
    sm = 0.0

    for x in c:
        if x == 'o(fe)':
            mol[x] = c[x] / cnst.m['o']
            sm += c[x] / cnst.m['o']
        else:
            mol[x] = c[x] / cnst.m[x]
            sm += c[x] / cnst.m[x]

    for x in c:
        mol[x] = mol[x] / sm

    if args:  # should be able to specify one element from the list to return if necessary.
        for arg in args:
            return mol[arg]
    else:
        return mol


def fo2_2F(cm,t,p,lnfo2,model="kc1991"):
    """
    Generate the fe3/fe2 ratio based on the set fO2 and fO2 model.

    Args:
        cm (dict): Dry silicate melt major element composition as a mole fractions.
        t (float): Temperatue (K)
        p (float): Pressure (Pa)
        lnfo2 (float): Ln(fO2)
        model (string): Name of fO2 model set in 'run'.

    Returns:
        Fe3/Fe2 mole fraction ratio
    """

    if model == "kc1991":
        return ferric.kc91_fo2(cm,t,p,lnfo2)
    elif model == "r2013":
        return ferric.r2013_fo2(cm, t, p, lnfo2)


def c2fo2(cm,t,p, model):
    """
    Calculates fO2 based on Fe3/Fe2 ratio from the melt major element composition and the fO2 model.

    Args:
        cm (dict): Dry silicate melt major element composition as a mole fractions.
        t (float): Temperatue (K)
        p (float): Pressure (Pa)
        model (string): Name of fO2 model set in 'run'.

    Returns:
        fO2 (float)
    """

    if model == "kc1991":
        return ferric.kc91_fe3(cm,t,p)
    elif model == "r2013":
        return ferric.r2013_fe3(cm, t, p)
    

def single_cat(c):
    "Convert dry weight fraction to single cation mol fraction. All iron as FeOt"

    mol = {}
    sm = 0

    for x in c:
        if x in ['al2o3', 'na2o', 'k2o']:
            mol[x] = c[x] / (cnst.m[x]*0.5)
            sm += c[x] / (cnst.m[x]*0.5)
        elif x in ['mgo', 'sio2', 'cao', 'tio2', 'mno', 'feo']:
            mol[x] = c[x] / cnst.m[x]
            sm += c[x] / cnst.m[x]

    for x in c:
        if x not in mol:
            pass
        else:
            mol[x] = mol[x] / sm

    return mol

# ------------------------------------------------------------------------
# CONVERSIONS TO AND FROM THE FMQ and IW BUFFERS FOR EASY COMPARISON TO LITERATURE
# ------------------------------------------------------------------------


def fmq2fo2(dfmq,t,p,name):  # dfmq = the value relative to the fmq buffer
    return(dfmq + fmq(t,p,name))

def fmq_2iw(dfmq, t, p, name):
    fo2 = fmq2fo2(dfmq, t, p, name)
    return fo2_2iw(fo2, t, p, name)

def fo2_2fmq(FO2,t,p,name):  # returns as fo2 relative to the FMQ buffer
    return(FO2 - fmq(t,p,name))

def fmq(t,p,name):
    if name == "frost1991":
        return(-25096.3/t + 8.735 + 0.11*(p-1)/t)  # Input pressure as bar, temp in Kelvin

def iw2fo2(FO2, t, p, name):
    return(FO2 + iw(t, p, name))

def iw_2fmq(diw, t, p, name):
    fo2 = iw2fo2(diw, t, p, name)
    return fo2_2iw(fo2, t, p, name)

def fo2_2iw(FO2, t, p, name):
    return (FO2 - iw(t, p, name))

def iw(t, p, name):
    if name == "frost1991":
        return (-27489/t + 6.702 + 0.055*(p-1)/t)

def nno2fo2(FO2, t, p, name):
    return(FO2 + nno(t, p, name))

def fo2_2nno(FO2, t, p, name):
    return (FO2 - nno(t, p, name))

def nno(t, p, name):
    if name == "frost1991":
        return (-24930/t + 9.36 + 0.046*(p-1)/t)

def generate_fo2(sys, dfo2, buffer, P):
    "Takes an fo2 value relative to any buffer, and returns ln(fo2) for use in sat_pressure."
    # PL: Edit to allow setting from an fe2/fe3 ratio somehow.

    if buffer == 'FMQ':
        fo2 = np.log(10 ** fmq2fo2(dfo2, sys.T, P, sys.run.FMQ_MODEL))  # ln fO2

    elif buffer == 'IW':
        fo2 = np.log(10 ** iw2fo2(dfo2, sys.T, P, sys.run.FMQ_MODEL))
 
    elif buffer == 'NNO':
        fo2 = np.log(10 ** nno2fo2(dfo2, sys.T, P, sys.run.FMQ_MODEL))

    return fo2

def generate_fo2_buffer(sys, fo2, P, buffer_choice = None):
    "Takes an fo2 value and returns the value relative to the chosen rock buffer."

    if buffer_choice is not None:
        buffer = buffer_choice
    elif sys.run.FO2_buffer_SET == True:
        buffer = sys.run.FO2_buffer
    else:
        return fo2_2fmq(np.log10(fo2), sys.T, P, sys.run.FMQ_MODEL)
    
    if buffer == 'FMQ':
        return fo2_2fmq(np.log10(fo2), sys.T, P, sys.run.FMQ_MODEL)

    elif buffer == 'IW':
        return fo2_2iw(np.log10(fo2), sys.T, P, sys.run.FMQ_MODEL)
    
    elif buffer == 'NNO':
        return fo2_2nno(np.log10(fo2), sys.T, P, sys.run.FMQ_MODEL)


# ------------------------------------------------------------------------
# CONVERSIONS FOR THE SOLVER SCRIPT
# ------------------------------------------------------------------------

def mols2wts(**kwargs):  # input as "H2O" = gas.mH2O[-1], "O2"...etc
    x = {}
    sum = 0

    for key, value in kwargs.items():
        sum += value * cnst.m[str(key.lower())]

    for key, value in kwargs.items():
        x[str(key)] = (value * cnst.m[str(key.lower())]) / sum

    return x

def mol2wt(c):
    mol = {}
    sm = 0

    for x in c:
        if x == 'ofe' or x == 'ogas':
            mol[x] = c[x] * cnst.m['o']
            sm += c[x] * cnst.m['o']
        else:
            mol[x] = c[x] * cnst.m[x]
            sm += c[x] * cnst.m[x]

    for x in c:
        mol[x] = mol[x] / sm

    return mol

def mean_mol_wt(**kwargs): # input as "H2O" = gas.mH2O[-1], "O2"...etc
    sum = 0
    
    for key, value in kwargs.items():
        sum += value * (cnst.m[str(key.lower())]/1000)  # in kg/mol

    return sum

# transform a value as a fraction into value as a percentage

def frac2perc(x):
    X = []
    for i in x:
        X.append(i * 100)
    return X

def atomicM_calc(sys, melt, gas, element, i, WgT=None):
    
    if WgT == None:
        WgT = sys.WgT[i+1]
    
    if sys.run.GAS_SYS == 'OH':

        if element == 'o':
                return cnst.m['o'] * (((WgT * gas.Wt['H2O'][i]) + (melt.h2o[i])) / cnst.m['h2o'] + 
                    (2 * WgT * gas.Wt['O2'][i]) / cnst.m['o2'])
        
        elif element == 'h':
            return 2 * cnst.m['h'] * (((WgT * gas.Wt['H2O'][i]) + (melt.h2o[i])) / cnst.m['h2o'] +
                        ((WgT * gas.Wt['H2'][i]) + (melt.h2[i])) / cnst.m['h2'])
        
        elif element == 'c' or element == 'n' or element == 's':
            return 0

        elif element == 'o_tot':
            if sys.run.FE_SYSTEM == False:
                return 0
            else:
                return (cnst.m['o'] * (((WgT * gas.Wt['H2O'][i]) + (melt.h2o[i])) / cnst.m['h2o'] + 
                    (2 * WgT * gas.Wt['O2'][i]) / cnst.m['o2'])) + melt.ofe[i]
    
    elif sys.run.GAS_SYS == 'COH':
    
        if element == 'o':
                return cnst.m['o'] * (((WgT * gas.Wt['H2O'][i]) + (melt.h2o[i])) / cnst.m['h2o'] + 
                    (2 * WgT * gas.Wt['O2'][i]) / cnst.m['o2'] + 
                        (WgT * gas.Wt['CO'][i] + melt.co[i])/cnst.m['co'] + 
                            (2 *((WgT * gas.Wt['CO2'][i]) + (melt.co2[i]))) / cnst.m['co2'])
        
        elif element == 'h':
            return 2 * cnst.m['h'] * (((WgT * gas.Wt['H2O'][i]) + (melt.h2o[i])) / cnst.m['h2o'] +
                        ((WgT * gas.Wt['H2'][i]) + (melt.h2[i])) / cnst.m['h2'] + 
                            2*(WgT * gas.Wt['CH4'][i] + melt.ch4[i])/cnst.m['ch4'])

        elif element == 'c':
            return cnst.m['c']*(((gas.Wt['CO2'][i]*WgT) + melt.co2[i])/cnst.m['co2'] + 
                    (gas.Wt['CO'][i]*WgT + melt.co[i])/cnst.m['co'] + 
                        (gas.Wt['CH4'][i]*WgT + melt.ch4[i])/cnst.m['ch4'] +
                        melt.graphite[i]/cnst.m['c'])

        elif element == 's' or element == 'n':
            return 0

        elif element == 'o_tot':
            if sys.run.FE_SYSTEM == False:
                return 0
            else:
                return (cnst.m['o'] * (((WgT * gas.Wt['H2O'][i]) + (melt.h2o[i])) / cnst.m['h2o'] + 
                    (2 * WgT * gas.Wt['O2'][i]) / cnst.m['o2'] + 
                        (WgT * gas.Wt['CO'][i] + melt.co[i])/cnst.m['co'] + 
                            (2 *((WgT * gas.Wt['CO2'][i]) + (melt.co2[i]))) / cnst.m['co2'])) + melt.ofe[i]
    
    elif sys.run.GAS_SYS == 'SOH':
        if element == 'o':
                return cnst.m['o'] * (((WgT * gas.Wt['H2O'][i]) + (melt.h2o[i])) / cnst.m['h2o'] + 
                    (2 * WgT * gas.Wt['O2'][i]) / cnst.m['o2'] + 
                            (2 *WgT * gas.Wt['SO2'][i]) / cnst.m['so2'])
        
        elif element == 'h':
            return 2 * cnst.m['h'] * (((WgT * gas.Wt['H2O'][i]) + (melt.h2o[i])) / cnst.m['h2o'] +
                        ((WgT * gas.Wt['H2'][i]) + (melt.h2[i])) / cnst.m['h2'] + 
                            (WgT * gas.Wt['H2S'][i]) / cnst.m['h2s'])
        
        elif element == 'c' or element == 'n':
            return 0

        elif element == 's':
            return cnst.m['s'] * ((2*gas.Wt['S2'][i]*WgT)/cnst.m['s2'] +
                            (gas.Wt['H2S'][i]*WgT)/cnst.m['h2s'] +
                            (gas.Wt['SO2'][i]*WgT)/cnst.m['so2'] +
                            melt.s[i]/cnst.m['s'])

        elif element == 'o_tot':
            if sys.run.FE_SYSTEM == False:
                return 0
            else:
                return (cnst.m['o'] * (((WgT * gas.Wt['H2O'][i]) + (melt.h2o[i])) / cnst.m['h2o'] + 
                    (2 * WgT * gas.Wt['O2'][i]) / cnst.m['o2'] + 
                            (2 *WgT * gas.Wt['SO2'][i]) / cnst.m['so2'])) + melt.ofe[i]
    
    elif sys.run.GAS_SYS == 'COHS' or sys.run.GAS_SYS == 'COHSN':

        if element == 'o':
                return cnst.m['o'] * (((WgT * gas.Wt['H2O'][i]) + melt.h2o[i]) / cnst.m['h2o'] + 
                    (2 * WgT * gas.Wt['O2'][i]) / cnst.m['o2'] + 
                        (WgT * gas.Wt['CO'][i] + melt.co[i])/cnst.m['co'] + 
                            (2 *((WgT * gas.Wt['CO2'][i]) + melt.co2[i])) / cnst.m['co2'] +
                            (2 *(WgT * gas.Wt['SO2'][i])) / cnst.m['so2'])
        
        elif element == 'o_tot':
            if sys.run.FE_SYSTEM == False:
                return 0
            else:
                return (cnst.m['o'] * (((WgT * gas.Wt['H2O'][i]) + melt.h2o[i]) / cnst.m['h2o'] + 
                    (2 * WgT * gas.Wt['O2'][i]) / cnst.m['o2'] + 
                        (WgT * gas.Wt['CO'][i] + melt.co[i])/cnst.m['co'] + 
                            (2 *((WgT * gas.Wt['CO2'][i]) + melt.co2[i])) / cnst.m['co2'] +
                            (2 *(WgT * gas.Wt['SO2'][i])) / cnst.m['so2'])) + melt.ofe[i]
        
        elif element == 'h':
            return 2 * cnst.m['h'] * (((WgT * gas.Wt['H2O'][i]) + melt.h2o[i]) / cnst.m['h2o'] +
                        ((WgT * gas.Wt['H2'][i]) + melt.h2[i]) / cnst.m['h2'] + 
                        (gas.Wt['H2S'][i]*WgT)/cnst.m['h2s'] +
                            2*((WgT * gas.Wt['CH4'][i]) + melt.ch4[i])/cnst.m['ch4'])

        elif element == 'c':
            return cnst.m['c']*(((gas.Wt['CO2'][i]*WgT) + melt.co2[i])/cnst.m['co2'] + 
                (gas.Wt['CO'][i]*WgT + melt.co[i])/cnst.m['co'] + 
                    (gas.Wt['CH4'][i]*WgT + melt.ch4[i])/cnst.m['ch4']
                    + melt.graphite[i]/cnst.m['c'])

        elif element == 's':
            return cnst.m['s'] * (2*(gas.Wt['S2'][i]*WgT)/cnst.m['s2'] +
                            (gas.Wt['H2S'][i]*WgT)/cnst.m['h2s'] +
                            (gas.Wt['SO2'][i]*WgT)/cnst.m['so2'] +
                            melt.s[i]/cnst.m['s'])

        elif element == 'n':
            if sys.run.GAS_SYS == 'COHSN':
                return melt.n[i] + cnst.m['n']*((2 * gas.Wt['N2'][i] * WgT)/cnst.m['n2'])
            else:
                return 0

def get_graphite(sys, melt, P, CO2, mCO, mCO2, mCH4, mO2, co2Y, coY, ch4Y, o2Y, N):

    """
    Returns the mass of graphite in the melt.
    """

    return ((sys.atomicM['c']/cnst.m['c']) - (N * (mCO + mCO2 + mCH4) + sl.co2_melt((co2Y*mCO2*P), CO2, (o2Y*mO2*P), sys.T, P, melt, name=sys.run.C_MODEL) +
                 sl.co_melt((coY*mCO*P), P, name = sys.run.CO_MODEL) + sl.ch4_melt((ch4Y*mCH4*P), P, name = sys.run.CH4_MODEL))) * cnst.m['c']

