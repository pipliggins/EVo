# constants.py

# Universal gas constant
R = 8.3144598

m = {}  # atomic and molecular masses  g/mol
m['h'] = 1.00794
m['o'] = 15.9994
m['c'] = 12.0000
m['s'] = 32.066
m['n'] = 14.0067
m['si'] = 28.0855
m['ti'] = 47.867
m['al'] = 26.981539
m['mn'] = 54.938044
m['mg'] = 24.305
m['ca'] = 40.078
m['na'] = 22.990
m['k'] = 39.0983
m['p'] = 30.97376
m['li'] = 6.941
m['h2'] = 2.0164
m['o2'] = 31.9988
m['h2o'] = 18.01528
m['co2'] = 44.0095
m['co'] = 28.0101
m['ch4'] = 16.0425
m['so2'] = 64.0638
m['h2s'] = 34.08088
m['s2'] = 64.13
m['ocs'] = 60.0751
m['n2'] = 28.0134
m['fe'] = 55.845
m['sio2'] = 60.0843
m['tio2'] = 79.8658
m['al2o3'] = 101.961278
m['fe2o3'] = 159.6882
m['feo'] = 71.8444
m['feot'] = 71.8444
m['mno'] = 70.937445
m['mgo'] = 40.3044
m['cao'] = 56.0774
m['na2o'] = 61.978938
m['k2o'] = 94.196
m['p2o5'] = 141.944524
m['li2o'] = 29.8814

#PL: change so reads as a_basalt, a_rhyolite, a_phonolite etc. These are stored in the Molecule instance for each species, so can specify at that stage which is relevant.
basalt = {
    'a' : {  # First solubility constant
        'H2O' : 6.576e-4,
        'H2' : 3.400e-7,
        'S2' : 0,  # this is sulfide capacity, controlled in the molecule class
        'CO2' : 1.729e-6,
        'O2' : 0,  # considered insoluble
        'CO' : 0,
        'CH4' : 0,
        'SO2' : 0,
        'H2S' : 0,
        'OCS' : 0,
        'N2': 0.0611e-6}, # physical dissolution

    'b' : {   # Second solubility constant, as above.
        'H2O' : 0.5698,
        'H2' : 1.2800,
        'S2' : 0.5,
        'CO2' : 0.8540,
        'O2' : 1,  # Considered insoluble
        'CO' : 1,
        'CH4' : 1,
        'SO2' : 1,
        'H2S' : 1,
        'OCS' : 1,
        'N2' : 5.97e-16} # chemical dissolution
        }

def p_h2o_a(T=1473.15):
    return -3.166e-9 * T**2 + 7.48e-6 * T - 3.853e-3

def p_h2o_b(T=1473.15):
    return 2.555e-6 * T**2 - 5.827e-3 * T + 3.918

#T-dependent solubilities can also be replaced with fixed constants outside of the tested temperature limits; this is controlled in the Molecule class.
phonolite = {

    'a' : {  # First solubility constant
        'H2O' : p_h2o_a,  
        'H2' : 3.400e-7,
        'SO2' : 2.019e-4,
        'H2S' : 4.172e-5,
        'CO2' : 4.339e-7,
        'O2' : 0,  # considered insoluble
        'CO' : 0,
        'CH4' : 0,
        'S2' : 0,
        'OCS' : 0,
        'N2': 0.0611e-6},

    'b' : {   # Second solubility constant, as above.
        'H2O' : p_h2o_b,
        'H2' : 1.2800,
        'SO2' : 0.4366,
        'H2S' : 0.5015,
        'CO2' : 0.8006,
        'O2' : 1,  # Considered insoluble
        'CO' : 1,
        'CH4' : 1,
        'S2' : 1,
        'OCS' : 1,
        'N2': 5.97e-16}
        }

def r_h2o_a(T=1473.15):
    return 2.5973e-8 * T**2 - 4.8473e-5 * T + 2.298e-2

def r_h2o_b(T=1473.15):
    return -5.1482e-6 * T**2 + 9.4853e-3 * T - 3.7085

def r_co2_a(T=1473.15):
    return 2.8895e-9 * T - 1.9625e-6

def r_co2_b(T=1473.15):
    return -1.0764e-3 * T + 1.9639

rhyolite = {
    'a' : {  # First solubility constant
        'H2O' : r_h2o_a,
        'H2' : 3.400e-7,
        'SO2' : 5.6322e-8,
        'H2S' : 2.3164e-6,
        'CO2' : r_co2_a,#
        'O2' : 0,  # considered insoluble
        'CO' : 0,
        'CH4' : 0,
        'S2' : 0,
        'OCS' : 0,
        'N2': 0.0611e-6},

    'b' : {   # Second solubility constant, as above.
        'H2O' : r_h2o_b,
        'H2' : 1.2800,
        'SO2' : 1.2937,
        'H2S' : 0.7338,
        'CO2' : r_co2_b,
        'O2' : 1,  # Considered insoluble
        'CO' : 1,
        'CH4' : 1,
        'S2' : 1,
        'OCS' : 1,
        'N2': 5.97e-16}
        }

PTcrit = {}
PTcrit['O2'] = [154.75, 50.7638]
PTcrit['H2O'] = [647.25, 221.1925]
PTcrit['H2'] = [33.25, 12.9696]
PTcrit['CO2'] = [304.15, 73.8659]
PTcrit['CH4'] = [191.05, 46.4069]
PTcrit['CO'] = [133.15, 34.9571]
PTcrit['S2'] = [208.15, 72.954]
PTcrit['SO2'] = [430.95, 78.7295]
PTcrit['OCS'] = [377.55, 65.8612]
PTcrit['H2S'] = [373.55, 90.0779]
PTcrit['N2'] = [126.2, 33.9]  # Roskosz 2006
