# ferric.py
# 
# methods for calcualting ferric and ferrous iron from melt composition and fO2
# do one or the other, assumes only one has been specified
# plus methods for calculating the iron fraction during decompression steps, for use in fe_equilibrate

import numpy as np

#------------------------------------------------------------------------
# FUNCTIONS
#------------------------------------------------------------------------

def kc91_fo2(C,T,P,lnfo2):
    """
    Calculates the Fe2O3/FeO mole ratio of a melt where the fO2 is known.

    Parameters
    ----------
    C : dictionary
        Major element composition of the silicate melt as MOLE FRACTIONS
        Required species: Al2O3, FeOt, CaO,  Na2O, K2O    
    T : float
        Temperature in degrees K    
    P : float
        Pressure in pascals (Pa)    
    lnfo2 : float
        ln(fO2)

    Returns
    -------
    F : float
        Fe2O3/FeO mole ratio

    References
    ----------
    Kress and Carmichael (1991) The compressibility of silicate liquids
    containing Fe2O3 and the effect of composition, temperature, oxygen 
    fugacity and pressure on their redox states.
    Contributions to Mineralogy and Petrology.
    """
    
    a = 0.196
    b = 1.1492e4                # K
    c = -6.675

    dal2o3 = -2.243
    dfeo = -1.828
    dcao = 3.201
    dna2o = 5.854
    dk2o = 6.215

    e = -3.36
    f = -7.01e-7                # K/Pa
    g = -1.54e-10               # /Pa
    h = 3.85e-17                # K/Pa^2

    T0 = 1673.0                 # K

    F = np.exp(a*lnfo2 + b/T + c +
               dal2o3*C['al2o3'] + dfeo*C['feo'] + dcao*C['cao'] + dna2o*C['na2o'] + dk2o*C['k2o'] + e*(1.0 - T0/T - np.log(T/T0)) + f*P/T + g*(T-T0)*P/T + h*P**2/T)

    # returns mole fraction Fe2O3/FeO
    return(F)

def kc91_fe3(Cm,T,P):
    """
    Calculates the oxygen fugacity (fO2) of a melt given where the ferric/
    ferrous ratio is known. 

    Parameters
    ----------
    C : dictionary
        Major element composition of the silicate melt as mole fractions
        Required species: Al2O3, FeO, Fe2O3, CaO,  Na2O, K2O    
    T : float
        Temperature in degrees K    
    P : float
        Pressure in pascals (Pa)

    Returns
    -------
    FO2 : float
        Absolute fO2
    
    References
    ----------
    Kress and Carmichael (1991) The compressibility of silicate liquids
    containing Fe2O3 and the effect of composition, temperature, oxygen
    fugacity and pressure on their redox states.
    Contributions to Mineralogy and Petrology.
    """
    a = 0.196
    b = 1.1492e4                # K
    c = -6.675

    dal2o3 = -2.243
    dfeo = -1.828
    dcao = 3.201
    dna2o = 5.854
    dk2o = 6.215

    e = -3.36
    f = -7.01e-7                # K/Pa
    g = -1.54e-10               # /Pa
    h = 3.85e-17                # K/Pa^2

    T0 = 1673.0                 # K

    FO2 = np.exp((np.log(Cm['fe2o3']/Cm['feo']) - b/T - c - dal2o3*Cm['al2o3'] - dfeo*(Cm['feo'] + Cm['fe2o3']*0.8998) - dcao*Cm['cao'] - dna2o*Cm['na2o'] - dk2o*Cm['k2o'] - e*(1.0 - T0/T - np.log(T/T0)) - f*(P/T) - g*(T-T0)*P/T - h*P**2/T)/a)

    return FO2

def r2013_fo2(cm, T, P, lnfo2):
    """
    Calculates the Fe2O3/FeO mole ratio of an FeOt>15 wt% melt where the
    fO2 is known.   

    Parameters
    ----------
    C : dictionary
        Major element composition of the silicate melt as mole fractions
        Required species: Al2O3, FeOt, CaO, Na2O, K2O, P2O5    
    T : float
        Temperature in degrees K    
    P : float
        Pressure in pascals (Pa)    
    lnfo2 : float
        ln(fO2)

    Returns
    -------
    F : float
        Fe2O3/FeO mole ratio
    
    References
    ----------
    Righter et al. (2013) Redox systematics of martian magmas with
    implications for magnetite stability. American Mineralogist.
    """
        
    a = 0.22
    b = 3800
    c = -370
    d = -6.6    # feo
    e = 7.3     # al2o3
    f = 17.3    # cao
    g = 132.3   # na2o
    h = -147.8  # k2o
    i = 0.6     # p2o5
    j = -4.26

    Pgpa = P*1e-9

    F = np.exp(a * lnfo2 + b / T + c * (Pgpa / T) + d * cm['feo'] + e * cm['al2o3'] + f * cm['cao'] + g * cm['na2o'] + h * cm['k2o'] + i * cm['p2o5'] + j)
    
    return F

def r2013_fe3(Cm, T, P):
    """
    Calculates the oxygen fugacity (fO2) of an FeOt>15 wt% melt given
    where the ferric/ferrous ratio is known.

    Parameters
    ----------
    C : dictionary
        Major element composition of the silicate melt as mole fractions
        Required species: Al2O3, FeOt, CaO, Na2O, K2O, P2O5    
    T : float
        Temperature in degrees K    
    P : float
        Pressure in pascals (Pa)

    Returns
    -------
    FO2 : float
        Absolute fO2

    References
    ----------
    Righter et al. (2013) Redox systematics of martian magmas with 
    implications for magnetite stability.
    """
        
    a = 0.22
    b = 3800
    c = -370
    d = -6.6    # feo
    e = 7.3     # al2o3
    f = 17.3    # cao
    g = 132.3   # na2o
    h = -147.8  # k2o
    i = 0.6     # p2o5
    j = -4.26

    Pgpa = P*1e-9   #GPa

    FO2 = np.exp((np.log(Cm['fe2o3']/Cm['feo']) - b/T - c*(Pgpa / T)  - d*(Cm['feo'] + Cm['fe2o3']*0.8998) - e * Cm['al2o3'] - f * Cm['cao'] - g * Cm['na2o'] - h * Cm['k2o'] - i * Cm['p2o5'] - j)/a)
    
    return FO2


