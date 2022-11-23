# density.py
"""
Calculates density of melt according to changing temp, pressure and compositions.
From Sperra(2000) Physical properties of magmas, Encyclopedia of Volcanology
"""

from evo import constants as cnst

#------------------------------------------------------------------------
# FUNCTIONS 
#------------------------------------------------------------------------

def den_calc_spera2000(Cm,T,P):
    """
    Calculates density of melt calc according to composition, temperature, pressure. 
    Encyclopeida of Volcanoes, Spera 2000

    Args:
        Cm (dict): composition of the melt, optionally including H2O and CO2, as mole fractions.
        T (float): Current temperature, K
        P (float): Current pressure, pascals (Pa)

    Returns:
        melt density, in kg/m3
    """

    T0 = 1673  # Ref temp, Kelvin (1400 celsius)
    P0 = 1.e-4  # Ref pressure, GPa (1 bar)

    # reference volumes: 10^-6 m3/mol, 1673 K (mm^3 per mol) Table S5.1
    vols = {}
    vols['sio2'] = 26.86/10.**6
    vols['tio2'] = 23.16/10.**6
    vols['al2o3'] = 37.42/10.**6
    vols['fe2o3'] = 42.13/10.**6
    vols['feo'] = 13.65/10.**6
    vols['mgo'] = 11.69/10.**6
    vols['cao'] = 16.53/10.**6
    vols['na2o'] = 28.88/10.**6
    vols['k2o'] = 45.07/10.**6
    vols['li2o'] = 16.85/10.**6
    vols['h2o'] = 26.27/10.**6
    vols['co2'] = 33.00/10.**6

    # dV/dT|P: 10^-9 m3/mol.K       #Thermal expansion per K
    dvdt = {}
    dvdt['sio2'] = 0.0/10.**9
    dvdt['tio2'] = 7.24/10.**9
    dvdt['al2o3'] = 0.0/10.**9
    dvdt['fe2o3'] = 9.09/10.**9
    dvdt['feo'] = 2.92/10.**9
    dvdt['mgo'] = 3.27/10.**9
    dvdt['cao'] = 3.74/10.**9
    dvdt['na2o'] = 7.68/10.**9
    dvdt['k2o'] = 12.08/10.**9
    dvdt['li2o'] = 5.25/10.**9
    dvdt['h2o'] = 9.46/10.**9
    dvdt['co2'] = 0.0/10.**9

    # dV/dP: 10^-6 m3/mol.GPa       #Compressability per GPa of pressure
    dvdp = {}
    dvdp['sio2'] = -1.89/10.**6
    dvdp['tio2'] = -2.31/10.**6
    dvdp['al2o3'] = -2.26/10.**6
    dvdp['fe2o3'] = -2.53/10.**6
    dvdp['feo'] = -0.45/10.**6
    dvdp['mgo'] = 0.27/10.**6
    dvdp['cao'] = 0.34/10.**6
    dvdp['na2o'] = -2.40/10.**6
    dvdp['k2o'] = -6.75/10.**6
    dvdp['li2o'] = -1.02/10.**6
    dvdp['h2o'] = -3.15/10.**6
    dvdp['co2'] = 0.0/10.**6
 
    # normalise to available oxides
    sm = 0
    tmp_Cm = {}         #Composition in mole fractions
    for e in vols:
        if e in Cm.keys():
            sm += Cm[e]
    for e in vols:
        if e in Cm.keys():
            tmp_Cm[e] = Cm[e]/sm
        
    V = 0.
    for x in tmp_Cm:
        V += tmp_Cm[x]*(vols[x] + dvdt[x]*(T-T0) + dvdp[x]*(P/10.**9-P0))
    M = sumM(tmp_Cm)
    return((M/1000.)/V)


def sumM(Cm):
    M = 0.
    for e in Cm:
        M += cnst.m[e]*Cm[e]
    return(M)
    
