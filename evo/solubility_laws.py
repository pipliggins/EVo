"""
Module to store the solubility laws of each species present in the melt,
using the method specified in the run.yaml file.
"""

import numpy as np
from scipy.special import erf
import warnings
from scipy.optimize import fsolve
import math

import conversions as cnvs
import constants as cnst

# ------------------------------------------------------------------------------
# MODEL DEFINITIONS (Melt from gas & fugacity from melt)
# ------------------------------------------------------------------------------

def ardia2013(fCH4, P):
    """
    The weight fraction of CH4 in the melt using method of Ardia et.al., 2013

    Args:
        fCH4 (float): The fugacity of CH4 (bars)
        P (float): Pressure (bars)
    
    Returns:
        The weight fraction of CH4 in the melt.
    """
    P = P*1e-4      # Convert to GPa
    fCH4 = fCH4*1e-4
    
    k = np.exp(4.93 - (0.000193*P))

    return (k*fCH4)/1e6

def ardia2013_fugacity(melt_CH4, P):
    """
    The fugacity CH4 using method of Ardia et.al., 2013

    Args:
        melt_CH4 (float): The weight fraction of CH4 in the melt
        P (float): Pressure (bars)
    
    Returns:
        fCH4  Fugacity of CH4 (bar)
    """
    P = P*1e-4      # Convert to GPa
    
    k = np.exp(4.93 - (0.000193*P))

    fugacity = (melt_CH4*1e6)/k

    return fugacity*1e4 # convert fugacity from GPa to bar

def armstrong2015(fCO, P):
    """
    The weight fraction of CARBON in C=O complexes (assume mean stoichiometry is
    a single CO) in a reduced melt from eq. 10 of Armstrong et al., 2015.
    (Technically this is 'non-carbonate C')
    
    Args:
        fCO (float): The fugacity of CO in the gas phase
        P (float): Pressure (bars)
    
    Returns:
        The weight fraction of CO in the melt.
    """
    
    C = 10**(-0.738 + 0.876*np.log10(fCO) - 5.44e-5*P)
    melt_CO = C/1e6       # ppm to weight fraction

    return melt_CO

def armstrong2015_fugacity(melt_CO, P):
    """
    The fugacity of CO in the gas based on non-carbonate C in the melt.
    from eq. 10 of Armstrong et al., 2015.
    
    Args:
        melt_CO (float): weight fraction of CO in the melt
        P (float): Pressure (bars)
    
    Returns:
        fCO (float): FUgacity of CO (bar)
    """
    
    melt_CO = melt_CO*1e6

    fCO = 10**((np.log10(melt_CO) + 0.738 + 5.44e-5*P)/0.876)

    return fCO

def burguisser2015_co2(fCO2, CO2):
    """
    Returns the weight fraction of CO2 in the melt using the method applied in
    D-Compress (Burguisser et. al., 2015).

    Args:
        fCO2 (float): The CO2 fugacity
        CO2 (class): The CO2 instance of the Molecule class
    
    Returns:
        The weight fraction of CO2 in the melt.
    """    
    return (CO2.solCon['a']*(fCO2)**CO2.solCon['b'])

def burguisser2015_co2_fugacity(melt_co2, CO2):
    """
    Returns the fugacity of CO2 in the gas phase according to the solubility law applied in
    D-Compress (Burguisser et. al., 2015).

    Args:
        melt_co2 (float): The weight fraction of CO2 in the melt
        CO2 (class): The CO2 instance of the Molecule class
    
    Returns:
        The fugacity of CO2 in the gas phase
    """     
    return (melt_co2/CO2.solCon['a']) ** (1 / CO2.solCon['b'])

def burguisser2015_h2(mH2, H2, P, H2Y=None):
    """
    Returns the weight fraction of H2 in the melt using the method applied in
    D-Compress (Burguisser et. al., 2015).

    Args:
        mH2 (float): The mol fraction of H2 in the gas phase
        H2 (class): The H2 instance of the Molecule class
        P (float): Pressure (bars)
        H2Y (float): The H2 fugacity coefficient (if different to one stored in H2)
    
    Returns:
        The weight fraction of H2 in the melt.
    """    
    if H2Y == None:
        return (H2.solCon['a']*(H2.Y*mH2*P)**H2.solCon['b'])
    else:
        return (H2.solCon['a']*(H2Y*mH2*P)**H2.solCon['b'])

def burguisser2015_h2_fugacity(melt_h2, H2):
    """
    Returns the fugacity of H2 in the melt using the solubility law applied in
    D-Compress (Burguisser et. al., 2015).

    Args:
        mH2 (float): The weight fraction of H2 in the melt
        H2 (class): The H2 instance of the Molecule class
    
    Returns:
        The fugacity of H2 in the gas phase
    """    
    return (melt_h2/H2.solCon['a'])**(1/H2.solCon['b'])

def gaillard2003_h2(mH2, H2, P, melt, H2Y=None):
    """
    Returns the weight fraction of H2 in the melt using the method derived in
    Gaillard et. al., 2003 (returns g/m3; melt density converts from H2 density to
    weight fraction)
    Applicable for fH2 = 0.02 - 70 bar, and 573 - 1273 K

    Args:
        mH2 (float): The mol fraction of H2 in the gas phase
        H2 (class): The H2 instance of the Molecule class
        P (float): Pressure (bars)
        melt (class): Active instance of Melt class
        H2Y (float): The H2 fugacity coefficient (if different to one stored in H2)
    
    Returns:
        The weight fraction of H2 in the melt.
    """    
    
    melt_density = melt.rho(P = P*1e5)/1e3      # convert to units g/cm3 from kg/m3

    if H2Y == None:
        return (3.4e-7*(H2.Y*mH2*P)**1.28)/melt_density
    else:
        return (3.4e-7*(H2Y*mH2*P)**1.28)/melt_density

def gaillard2003_h2_fugacity(melt_h2, melt):
    """
    Returns the fugacity of H2 in the melt using the solubility law derived in
    Gaillard et. al., 2003.

    Args:
        melt_h2 (float): The weight fraction of H2 in the melt
        melt (class): The active instance of the Melt class
    
    Returns:
        The fugacity of H2 in the gas phase
    """   

    melt_density = melt.rho()/1e3

    return ((melt_h2*melt_density)/3.4e-7)**(1/1.28)

def burguisser2015_h2o(mH2O, H2O, P, H2OY=None):
    """
    Returns the weight fraction of H2O in the melt using the method applied in
    D-Compress (Burguisser et. al., 2015).

    Args:
        mH2O (float): The mol fraction of H2O in the gas phase
        H2O (class): The H2O instance of the Molecule class
        P (float): Pressure (bars)
        H2OY (float): The H2O fugacity coefficient (if different to version in H2O)
    
    Returns:
        The weight fraction of H2O in the melt.
    """    
    if H2OY == None:
        return (H2O.solCon['a']*(H2O.Y*mH2O*P)**H2O.solCon['b'])
    else:
        return (H2O.solCon['a']*(H2OY*mH2O*P)**H2O.solCon['b'])

def burguisser2015_h2o_fugacity(melt_h2o, H2O):
    """
    Returns the fugacity of H2O in the gas phase using the solubility law applied in
    D-Compress (Burguisser et. al., 2015).

    Args:
        melt_h2o (float): The weight fraction of H2O in the melt
        H2O (class): The H2O instance of the Molecule class
    
    Returns:
        The fugacity of H2O in the gas.
    """    
    return (melt_h2o/H2O.solCon['a'])**(1/H2O.solCon['b'])

def eguchi2018(fCO2, fO2, T, P, melt):
    """
    Returns the weight fraction of CO2 dissolved in the melt, as a converted sum of
    both molecular CO2, and carbonate (CO3^2-)
    
    Args:
        fCO2 (float): Fugacity of CO2 (bar)
        fO2 (float): Oxygen fugacity (bar)
        T (float): Temperature (K)
        P (float): Pressure (bar)
        melt (class): active instance of Melt class

    Returns:
        Weight fraction of CO2 dissolved in the melt, as the converted sum for
        CO2(mol) + CO3
    """

    def NBO(mol_fracs):
        nsi = mol_fracs['sio2']
        nti = mol_fracs['tio2']
        nal = mol_fracs['al2o3'] * 2
        nfe2 = mol_fracs['feo']
        nfe3 = mol_fracs['fe2o3'] * 2
        nmn = mol_fracs['mno']
        nmg = mol_fracs['mgo']
        nca = mol_fracs['cao']
        nna = mol_fracs['na2o'] * 2
        nk = mol_fracs['k2o'] * 2
        np = mol_fracs['p2o5'] * 2
        o = mol_fracs['sio2'] * 2 + mol_fracs['tio2'] * 2 + mol_fracs['al2o3'] * 3 + nfe2 + nfe3 * 1.5 + mol_fracs['mno'] + mol_fracs['mgo']\
                    + mol_fracs['cao'] + mol_fracs['na2o'] + mol_fracs['k2o'] + mol_fracs['p2o5'] * 5


        NM = (nmg + nca + nfe2 + nna + nk + nmn)

        Al = nal - NM

        if Al > 0:
            al_tet = NM
        else:
            al_tet = nal

        Fe = nfe3 + Al

        if Al > 0:
            fe_tet = 0
        elif Al <= 0 and Fe > 0:
            fe_tet = Al*-1
        elif Al <= 0 and Fe <= 0:
            fe_tet = nfe3

        Tet = nsi + nti + np + al_tet + fe_tet
        NBO = 2*o - 4*Tet
        return NBO
    
    def xi(fCO2, T, P, nbo, oxides, name='co2'):
        
        cao = oxides['cao']
        na2o = oxides['na2o']
        k2o = oxides['k2o']      
        
        if name =='co2':
            DH = -90212
            DV = 0.000019244
            DS = -43.0815
            B = 1114.9
            yNBO = -7.0937
            A_CaO = 0
            A_Na2O = 0
            A_K2O = 0
            
        elif name == 'co3':
            DH = -164480
            DV = 0.00002384
            DS = -43.6385
            B = 1473.2
            yNBO = 3.291
            A_CaO = 1.68e5
            A_Na2O = 1.759e5
            A_K2O = 2.1085e5
        
        # PL: This isn't great really
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            lnXi = -(DV * (P*1e5) / (8.3144 * T)) + DH / (8.3144 * T) + (np.log(fCO2) * B) / T + DS / 8.3144 + (yNBO * nbo) + \
                    (A_CaO * cao + A_Na2O * na2o + A_K2O * k2o) / (8.3144 * T)
        
        return np.exp(lnXi)

    oxides = melt.iron_fraction(np.log(fO2), ppa = P*1e5)[0]        # dry melt as oxide mol fractions

    nbo = NBO(oxides)
    
    Xco2 = xi(fCO2, T, P, nbo, oxides, name='co2')
    Xco3 = xi(fCO2, T, P, nbo, oxides, name='co3')

    FW_one = melt.formula_weight(fO2, P)

    # This assumes that the only species present in the melt is CO2/CO3. Not completely accurate here.
    CO2_CO2 = ((44.01*Xco2)/(44.01*Xco2+(1-(Xco2+Xco3))*FW_one))
    CO2_CO3 = ((44.01*Xco3)/(44.01*Xco3+(1-(Xco2+Xco3))*FW_one))

    return CO2_CO2 + CO2_CO3

def eguchi2018_fugacity(melt_co2, fO2, T, P, melt):
    """
    Calculates fCO2 from the total melt CO2 content (melt fraction) using 
    Eguchi & Dasgupta 2018.

    Uses fsolve to find the partitioning between molecular CO2 & CO3^2-,
    then the corresponding fCO2.

    Args:
        melt_CO2 (float): Total melt CO2 content (melt frac)
        fO2 (float): Oxygen fugacity
        T (float): Temperature (K)
        P (float): Pressure (bar)
        melt (class): active instance of Melt class

    Returns:
        fCO2 (float): CO2 fugacity
    """

    def f(fco2, fO2, T, P, melt, melt_CO2):
        return eguchi2018(fco2, fO2, T, P, melt) - melt_CO2
    
    fCO2 = fsolve(f, 1.0, args=(fO2, T, P, melt, melt_co2))[0]

    # Check for graphite saturation
    logK1 = 40.07639 - 2.53932 * 10 ** -2 * T + 5.27096 * 10 ** -6 * T ** 2 + 0.0267 * (P - 1) / T
    graph_fco2 = 10 ** logK1 * 10 ** np.log10(fO2)

    if fCO2 > graph_fco2:
        melt.graphite_sat = True
        return graph_fco2
    else:
        melt.graphite_sat = False
        return fCO2

def graphite_fco2(T, P, fO2):
    """
    Calculates the fCO2 in equilibrium with a graphite sturated melt.
    """

    logK1 = 40.07639 - 2.53932 * 10 ** -2 * T + 5.27096 * 10 ** -6 * T ** 2 + 0.0267 * (P - 1) / T
    graph_fco2 = 10 ** logK1 * 10 ** np.log10(fO2)

    return graph_fco2

def libourel2003(mN2, fO2, P):
    """
    Returns the weight fraction of N in the melt.

    Args:
        mN2 (float): The mol fraction of N2 in the gas phase
        mO2 (float): The mol fraction of O2 in the gas phase
        O2 (class): The O2 instance of the Molecule class
        P (float): Pressure (bars)
    
    Returns:
        The weight fraction of N in the melt.
    """
    return (0.0611e-6*P*0.986923*mN2 + 5.97e-16 * fO2**-0.75 * (P*0.986923*mN2)**0.5)

def libourel2003_fugacity(n_melt, nY, fO2, P):
    """
    Returns the fugacity of N2 in the gas phase.

    Args:
        n_melt (float): The weight fraction of N in the melt
        nY (float): The fugacity coefficient of N2
        fO2 (float): The absolute oxygen fugacity
        P (float): Pressure (bars)
    
    Returns:
        The fugacity of N2 in the gas phase
    """
   
    def n_quadratic(n_melt):
        
        a = 0.0611e-6
        b = 5.97e-16

        def f0(x):
            return a*x + b*((fO2)**(-0.75))*(x**0.5) - n_melt
        
        with warnings.catch_warnings():
            warnings.filterwarnings('error')
            try:
                PN2 = fsolve(f0, 1)[0]
            except:
                try:
                    PN2 = fsolve(f0, 1e-10)[0]
                except:
                    raise RuntimeError('Failed to find N2 partial pressure.')

        return PN2

    PN2 = n_quadratic(n_melt)
    mN2 = PN2/(P*0.986923)  # pressure in atmospheres
    fn2 = nY * mN2 * P

    return fn2

def nash2019(fO2, P, T, melt, run):
    """
    Calculates the ratio of sulfate:sulfide (S6+/S2-) in the melt
    after Nash et al., 2019.

    Args:
        fO2 (float): Absolute fO2
        P (float): Pressure (bar)
        T (float): Temperature (K)
        melt (class): active instance of Melt class
        run (class): active instance of the RunDef class
    
    Returns:
        ratio (float): S6+/S2- ratio of sulfur species in the melt
    """
    F = 2*cnvs.fo2_2F(melt.Cm(), T, P*1e5, np.log(fO2), run.FO2_MODEL)  # F is mFe2O3/mFeO. *2 gives fe3/fe2.

    ratio = 10**(8*np.log10(F) + (8.7436e6/T**2) - (27703/T) + 20.273)

    return ratio

def oneill2002(fO2, P, melt):
    """
    Calculates the sulfide capacity of a melt based on the methods of
    O'Neill et al., 2002.

    Args:
        fO2 (float): Absolute fO2
        P (float): Pressure (bar)
        melt (class): active instance of Melt class
    
    Returns:
        capacity (float): Sulfide (S2-) capacity of melt as a weight fraction
    """
    mol_frac = melt.iron_fraction(np.log(fO2), ppa=P*1e5)[0]
    wts = cnvs.mol2wt(mol_frac)
    feo = wts['feo']*100  # wt% of feo in melt
                
    capacity = (0.0003*(100-feo)*np.exp(0.21*feo))/1000000  # Convert ppm -> wt fraction
    return capacity

def oneill2020(T, melt):
    """
    Calculates the sulfide capacity of a melt based on the methods of
    O'Neill et al., 2020.

    Args:
        T (float): temperature, K
        melt (class): active instance of Melt class
    
    Returns:
        capacity (float): Sulfide (S2-) capacity of melt as a weight fraction
    """
    if not bool(melt.cm_dry):
        comp = melt.Cm()
    else:
        comp = melt.cm_dry

    if comp['fe2o3'] > 0:
        comp['feo'] = comp['feo'] + comp['fe2o3']*0.8998
        comp['fe2o3'] = 0.0
    
    comp = cnvs.mol2wt(comp)
    comp = cnvs.single_cat(comp)  # Converts the dry mol fraction to normalized mol fraction on single cation basis.
    
    Na = comp['na2o']
    Ti = comp['tio2']
    K = comp['k2o']
    Ca = comp['cao']
    Mg = comp['mgo']
    Si = comp['sio2']
    Al = comp['al2o3']
    Fe2 = comp['feo'] # FeOt
    Mn = comp['mno']
    FeX = Fe2 + Mn
    
    capacity = (np.exp(-23590/T + 8.77 +(1673/T)*(6.7*(Na+K) + 1.8*Al + 4.9*Mg + 8.1*Ca + 5*Ti + 8.9*FeX - 22.2*FeX*Ti + 7.2*FeX*Si) - 2.06*erf(-7.2*FeX)))/1000000  # Convert ppm -> wt fraction

    return capacity


# ------------------------------------------------------------------------------
# MELT CONTENT MODEL SELECTION
# ------------------------------------------------------------------------------

def ch4_melt(fCH4, P, name='ardia2013'):
    """
    Returns the number of moles of CH4 in the melt.
    Applicable only for low fO2 conditions.

    Args:
        fCH4 (float): the fugacity of CH4 in the gas
        P (float): Pressure (bar)
        name (string): The name of the solubility law to be used, taken from 'run'

    Returns:
    melt_ch4 (float): the number of moles of CH4 in the melt.
    """
    
    if name == 'ardia2013':
        return ardia2013(fCH4, P)/cnst.m['ch4']
    
    elif name == 'None':
        return 0

def co_melt(fCO, P, name='armstrong2015'):
    """
    Returns the number of moles of CO in the melt.
    Applicable only for low fO2 conditions.

    Args:
        fCO (float): the fugacity of CO in the gas
        P (float): Pressure (bar)
        name (string): The name of the solubility law to be used, taken from 'run'

    Returns:
    melt_co (float): the number of moles of CO in the melt.
    """
    
    if name == 'armstrong2015':
        return armstrong2015(fCO, P)/cnst.m['co']
    
    elif name == 'None':
        return 0

def co2_melt(fCO2, CO2, fO2, T, P, melt, name='burguisser2015'):
    """
    Returns the number of moles of CO2 in the melt.

    Args:
        fCO2 (float): CO2 fugacity
        CO2 (class): CO2 instance of the Molecule class
        fO2 (float): Oxygen fugacity
        T (float): Temperature (K)
        P (float): Pressure (bars)
        melt (class): Active instance of the Melt class
        name (string): The name of the solubility law to be used, taken from 'run'
    
    Returns:
        melt_CO2 (float): the number of moles of CO2 in the melt.
    """
    
    if name == 'burguisser2015':
        return burguisser2015_co2(fCO2, CO2)/cnst.m['co2']
    
    elif name == 'eguchi2018':
        return eguchi2018(fCO2, fO2, T, P, melt)/cnst.m['co2']

def h2_melt(mH2, H2, P, melt, name='burguisser2015', Y = None):
    """
    Returns the number of moles of H2 in the melt.

    Args:
        mH2 (float): The mol fraction of H2 in the gas phase
        H2 (class): The H2 instance of the Molecule class
        P (float): Pressure (bars)
        melt (class): Active instance of Melt class
        name (string): The name of the solubility law to be used, taken from 'run'
        Y (float): The H2 fugacity coefficient, if different to version stored in H2
    
    Returns:
        melt_H2 (float): the number of moles of H2 in the melt.
    """
    
    if name == 'burguisser2015':
        return burguisser2015_h2(mH2, H2, P, H2Y=Y)/cnst.m['h2']

    elif name == 'gaillard2003':
        return gaillard2003_h2(mH2, H2, P, melt, H2Y=Y)/cnst.m['h2']

def h2o_melt(mH2O, H2O, P, name='burguisser2015', Y = None):
    """
    Returns the number of moles of H2O in the melt.

    Args:
        mH2O (float): The mol fraction of H2O in the gas phase
        H2O (class): The H2O instance of the Molecule class
        P (float): Pressure (bars)
        name (string): The name of the solubility law to be used, taken from 'run'
        Y (float): The H2O fugacity coefficient, if different to version stored in H2O
    
    Returns:
        melt_H2O (float): the number of moles of H2O in the melt.
    """
    
    if name == 'burguisser2015':
        return burguisser2015_h2o(mH2O, H2O, P, H2OY=Y)/cnst.m['h2o']

def n_melt(mN2, fO2, P, name='libourel2003'):
    """
    Returns the number of moles of N in the melt.

    Args:
        mN2 (float): The mol fraction of N2 in the gas phase
        mO2 (float): The mol fraction of O2 in the gas phase
        O2 (class): The O2 instance of the Molecule class
        P (float): Pressure (bars)
        name (string): The name of the solubility law to be used, taken from 'run'
    
    Returns:
        The number of moles of N in the melt.
    """
    if name == 'libourel2003':
        return libourel2003(mN2, fO2, P)/cnst.m['n']

def sulfate_melt(fS2, fO2, P, T, melt, run, name = 'nash2019'):
    """
    Returns the number of moles of S6+ in the melt.

    Args:
        mS2 (float): The mol fraction of S2 in the gas phase
        S2 (class): The S2 instance of the Molecule class
        mO2 (float): The mol fraction of O2 in the gas phase
        O2 (class): The O2 instance of the Molecule class
        P (float): Pressure (bars)
        T (float): Temperature (K)
        melt (class): active instance of the Melt class
        run (class): active instance of the RunDef class
        name (string): The name of the solubility law to be used, taken from 'run'
    
    Returns:
        melt_S6+ (float): the weight fraction of sulfate (S6+) in the melt.
    """

    if name == 'nash2019':
        
        ratio = nash2019(fO2, P, T, melt, run)
        
        return ratio * (sulfide_melt(fS2, fO2, P, T, melt, name = run.SULFIDE_CAPACITY))

def sulfide_melt(fS2, fO2, P, T, melt, name = 'oneill2020'):
    """
    Returns the number of moles of S2- in the melt.

    Args:
        mS2 (float): The mol fraction of S2 in the gas phase
        S2Y (float): S2 fugacity coefficient
        mO2 (float): The mol fraction of O2 in the gas phase
        O2Y (float): O2 fugacity coefficient
        P (float): Pressure (bars)
        T (float): Temperature (K)
        melt (class): active instance of the Melt class
        name (string): The name of the solubility law to be used, taken from 'run'
    
    Returns:
        melt_S2- (float): the weight fraction of sulfide (S2-) in the melt.
    """
    if name == 'oneill2002':                    
        capacity = oneill2002(fO2, P, melt)
    
    elif name == 'oneill2020':
        capacity = oneill2020(T, melt)
    
    return (capacity * (fS2/fO2) ** 0.5)/cnst.m['s']


# ------------------------------------------------------------------------------
# FUGACITY MODEL SELECTION
# ------------------------------------------------------------------------------

def ch4_fugacity(melt_ch4, P, name='ardia2013'):
    """
    Returns the fugacity of CH4.
    Applicable only for low fO2 conditions.

    Args:
        melt_ch4 (float): the weight fraction of CH4 in the melt
        P (float): Pressure (bar)
        name (string): The name of the solubility law to be used, taken from 'run'

    Returns:
    fCH4 (float): the fugacity of CH4
    """
    if name == 'ardia2013':
        return ardia2013_fugacity(melt_ch4, P)
    
    elif name == 'None':
        return 0

def co_fugacity(melt_co, P, name='armstrong2015'):
    """
    Returns the fugacity of CO.
    Applicable only for low fO2 conditions.

    Args:
        melt_co (float): the weight fraction of CO in the melt
        P (float): Pressure (bar)
        name (string): The name of the solubility law to be used, taken from 'run'

    Returns:
    fCO (float): the fugacity of CO
    """
    if name == 'armstrong2015':
        return armstrong2015_fugacity(melt_co, P)
    
    elif name == 'None':
        return 0

def co2_fugacity(melt_co2, CO2, fO2, T, P, melt, name='burguisser2015'):
    """
    Returns the fugacity of CO2 in the gas.

    Args:
        melt_co2 (float): The weight fraction of CO2 in the melt
        CO2 (class): The CO2 instance of the Molecule class
        fO2 (float): Oxygen fugacity NOT NAT LOGGED
        name (string): The name of the solubility law to be used, taken from 'run'
    
    Returns:
        fCO2 (float): the fugacity of CO2
    """
    
    if name == 'burguisser2015':
        return burguisser2015_co2_fugacity(melt_co2, CO2)
    
    elif name == 'eguchi2018':
        return eguchi2018_fugacity(melt_co2, fO2, T, P, melt)

def h2_fugacity(melt_h2, H2, melt, name='burguisser2015'):
    """
    Returns the fugacity of H2 in the gas.

    Args:
        melt_h2 (float): The weight fraction of H2 in the melt
        H2 (class): The H2 instance of the Molecule class
        melt (class): Active instance of Melt class
        name (string): The name of the solubility law to be used, taken from 'run'
    
    Returns:
        fH2 (float): the fugacity of H2
    """
    
    if name == 'burguisser2015':
        return burguisser2015_h2_fugacity(melt_h2, H2)
    
    elif name == 'gaillard2003':
        return gaillard2003_h2_fugacity(melt_h2, melt)

def h2o_fugacity(melt_h2o, H2O, name='burguisser2015'):
    """
    Returns the fugacity of H2O in the gas.

    Args:
        melt_h2o (float): The weight fraction of H2O in the melt
        H2O (class): The H2O instance of the Molecule class
        name (string): The name of the solubility law to be used, taken from 'run'
    
    Returns:
        fH2O (float): the fugacity of H2O.
    """
    
    if name == 'burguisser2015':
        return burguisser2015_h2o_fugacity(melt_h2o, H2O)

def n2_fugacity(melt_n, N2Y, fO2, P, name='libourel2003'):
    """
    Returns the fugacity of N2 in the gas.

    Args:
        melt_n (float): The weight fraction of N in the melt
        N2Y (float): The fugacity coefficient of N2
        fo2 (float): The absolute fO2
        P (float): Pressure (bars)
        name (string): The name of the solubility law to be used, taken from 'run'
    
    Returns:
        The fugacity of N2
    """
    
    if name == 'libourel2003':
        return libourel2003_fugacity(melt_n, N2Y, fO2, P)

def S2_fugacity(s_melt, fO2, P, T, melt, run, sulfidename = 'oneill2020', sulfatename='nash2019'):
    """
    Returns the fugacity of S2 in the gas.

    Args:
        s_melt (float): The weight fraction of total sulfur (S2- + S6+) in the melt
        fO2 (float): The absolute fO2
        P (float): Pressure (bars)
        T (float): temperature (K)
        melt (class): active instance of the Melt class
        run (class): active instance of the RunDef class
        sulfidename (string): The name of the sulfide capacity law
        sulfatename (string): The name of the sulfate calculation
    
    Returns:
        fS2: the fugacity of S2 in the gas.
    """

    if sulfatename == 'nash2019':
        melt.cm_dry = melt.iron_fraction(np.log(fO2), ppa=P*1e5)[0] # Uses fO2 to set the Fe2/Fe3 ratio and dry melt chemistry prior to needing sulfide capacity
        
        # need to find the ratio of s6+/s2-, then use the amount of s2- to get mS2.
        ratio = nash2019(fO2, P, T, melt, run) # s6+/s2-
        s2_melt = s_melt * (1/(1+ratio))

        if sulfidename == 'oneill2002':                    
            capacity = oneill2002(fO2, P, melt)
    
        elif sulfidename == 'oneill2020':
            capacity = oneill2020(T, melt)

        fS2 = (s2_melt/capacity)**(1/0.5) * fO2

        return fS2

