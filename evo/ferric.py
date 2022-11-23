# ferric.py
# 
# methods for calculating ferric and ferrous iron from melt composition and fO2
# do one or the other, assumes only one has been specified
# plus methods for calculating the iron fraction during decompression steps, for use in fe_equilibrate

import numpy as np

#------------------------------------------------------------------------
# FUNCTIONS
#------------------------------------------------------------------------

def kc91_fo2(C,T,P,FO2):  # If fo2 was set
    " Returns the mFe2O3/mFeO ratio; pressure in Pa"
    # Kress and Carmichael 1991, CMP
    # FO2 is imported as ln(fo2)
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

    F = np.exp(a*FO2 + b/T + c +
               dal2o3*C['al2o3'] + dfeo*C['feo'] + dcao*C['cao'] + dna2o*C['na2o'] + dk2o*C['k2o'] + e*(1.0 - T0/T - np.log(T/T0)) + f*P/T + g*(T-T0)*P/T + h*P**2/T)

    # returns mole fraction Fe2O3/FeO
    return(F)

def kc91_fe3(Cm,T,P):  # if the ratio was set
    "Returns fO2 when provided with fe2o3/feo ratio and composition as mol fraction; pressure in Pa"
    # Kress and Carmichael 1991, CMP
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
    Returns the Fe3/Fe2 ratio for a melt where the fO2 is specified, for a high Fe melt.
    
    Takes the melt composition as normalised mole fractions, with all iron as feo(t), pressure as Pa, converted to GPa
    
    Using Righter et al., 2013; valid for FeO(t) values above 15 wt% (Martian shergottities)
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
    Returns the fO2 for a melt where the fe3/fe2 ratio is specified, for a high Fe melt.
    
    Takes the melt composition as normalised mole fractions, P in Pa then converted to GPa.
    
    Using Righter et al., 2013; valid for FeO(t) values above 15 wt% (Martian shergottities)
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


