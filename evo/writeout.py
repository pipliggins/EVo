"""Writes the results of the run to a file and contains options to produce a graph of the results."""
import conversions as cnvt
import constants as cnst
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sat_pressure as sat
import warnings

def writeout_file(sys, gas, melt, P):
    if not os.path.exists('Output'):
        os.makedirs('Output')
    
    output = open('Output/dgs_output.csv', 'w')
    output.write(f"#Decompressing a {sys.run.COMPOSITION} {sys.run.GAS_SYS} {sys.run.RUN_TYPE} system at {sys.T:.2f} degrees (K). \n"
                 f"#Approx total H2O (wt%) starting fO2: {(gas.Wt['H2O'][0]*sys.WgT[0] + (melt.h2o[0]))*100:.3}, final fO2: {(gas.Wt['H2O'][-1]*sys.WgT[-1] + (melt.h2o[-1]))*100:.3}\n"
                 f"#Approx total CO2 starting fO2: {(gas.Wt['CO2'][0]*sys.WgT[0] + (melt.co2[0]))*100:.3} (wt%) {(gas.Wt['CO2'][0]*sys.WgT[0] + (melt.co2[0]))*1000000:.0} (ppm), final fO2: {((gas.Wt['CO2'][-1]*sys.WgT[-1] + (melt.co2[-1]))*100)*10000:.0} (ppm) \n"
                 f"#Approx total C (ppm) starting fO2: {(((gas.Wt['CO2'][0]*sys.WgT[1]) + (melt.co2[0]))/cnst.m['co2'] + (gas.Wt['CO'][0]*sys.WgT[1])/cnst.m['co'] + (gas.Wt['CH4'][0]*sys.WgT[1])/cnst.m['ch4'])*1000000*cnst.m['c']}, final fO2: {(((gas.Wt['CO2'][-1]*sys.WgT[-1]) + (melt.co2[-1]))/cnst.m['co2'] + (gas.Wt['CO'][-1]*sys.WgT[-1])/cnst.m['co'] + (gas.Wt['CH4'][-1]*sys.WgT[-1])/cnst.m['ch4'])*1000000*cnst.m['c']}\n"
                 f"#Total H (ppm): {sys.atomicM['h']*1000000} wt%: {sys.atomicM['h']*100}\n"
                 f"#fTH2 = {gas.mH2[-1]+gas.mH2S[-1]+(2*gas.mCH4[-1])}\n")

    output.write("#\n")
    output.write(f"#{'':8} \t {'':8} \t {'':8} \t {'':8} \t {'':12} \t {'':12} \t {'Gas:':12} \t {'':10} \t {'':14} \t {'Gas by mol:':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'Gas by mass:':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'Melt (wt%):':8} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'Ratios:':8}\n")
    
    output.write(f"{'P'}\t{sys.run.FO2_buffer}\t{'fo2'}\t{'F'}\t{'rho_bulk'}\t{'rho_melt'}\t{'Exsol_vol%'}\t{'Gas_wt'}\t{'mol_mass'}\t{'mH2O'}\t{'mH2'}\t{'mO2'}\t{'mCO2'}\t{'mCO'}\t{'mCH4'}\t{'mSO2'}\t{'mH2S'}\t{'mS2'}\t{'mN2'}\t{'wH2O'}\t{'wH2'}\t{'wO2'}\t{'wCO2'}\t{'wCO'}\t{'wCH4'}\t{'wSO2'}\t{'wH2S'}\t{'wS2'}\t{'wN2'}\t{'H2O_melt'}\t{'H2_melt'}\t{'CO2_melt'}\t{'CO_melt'}\t{'CH4_melt'}\t{'graph_melt'}\t{'S2-_melt'}\t{'S6+_melt'}\t{'Stot_melt'}\t{'N_melt'}\t{'mCO2/CO'}\t{'mCO2/H2O'}\t{'mCO2/SO2'}\t{'mH2S/SO2'}\t{'fH2'}\t{'fH2O'}\t{'fCO2'}\t{'fCO'}\t{'fCH4'}\t{'fSO2'}\t{'fH2S'}\t{'fS2'}\t{'fN2'}\t{'tot_H'}\t{'tot_C'}\t{'tot_O_gas'}\t{'tot_O'}\t{'tot_S'}\t{'tot_N'}\n")

    i = 0

    if sys.run.FIND_SATURATION == True or sys.run.ATOMIC_MASS_SET == True:
        output.write(sat.satp_writeout(sys, melt, gas, sys.sat_conditions[0], sys.sat_conditions[1], sys.sat_conditions[2], sys.sat_conditions[3], graph_sat=sys.sat_conditions[4]))
    
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        if sys.run.SINGLE_STEP == False:

            while i < len(P):
                output.write(f"{P[i]} \t {cnvt.generate_fo2_buffer(sys, np.exp(gas.fo2[i]), P[i]):+.3} \t {np.exp(gas.fo2[i]):.3e} \t {melt.F[i+1]:.3} \t {sys.rho[i]:.3} \t {melt.rho_store[i]:.3} \t {sys.GvF[i] * 100:} \t {sys.WgT[i+1] * 100:.5g} \t {gas.M[i]} \t {gas.mH2O[i+1]:.8g} \t {gas.mH2[i+1]:.8g} \t {gas.mO2[i+1]:.6g} \t {gas.mCO2[i+1]:.6g} \t {gas.mCO[i+1]:.6g} \t {gas.mCH4[i+1]:.6g} \t {gas.mSO2[i+1]:.6g} \t {gas.mH2S[i+1]:.6g} \t {gas.mS2[i+1]:.6g} \t {gas.mN2[i+1]:.6g} \t {gas.Wt['H2O'][i]:.6g} \t {gas.Wt['H2'][i]:.6g} \t {gas.Wt['O2'][i]:.6g} \t {gas.Wt['CO2'][i]:.6g} \t {gas.Wt['CO'][i]:.6g} \t {gas.Wt['CH4'][i]:.6g} \t {gas.Wt['SO2'][i]:.6g} \t {gas.Wt['H2S'][i]:.6g} \t {gas.Wt['S2'][i]:.6g} \t {gas.Wt['N2'][i]:.6g} \t {melt.h2o[i]*100:.6g} \t {melt.h2[i]*100:.6g} \t {melt.co2[i]*100:.6g} \t {melt.co[i]*100:.6g} \t {melt.ch4[i]*100:.6g} \t {melt.graphite[i]*100:.6g} \t {melt.sulfide[i]*100:.6g} \t {melt.sulfate[i]*100:.6g} \t {melt.s[i]*100:.6g} \t {melt.n[i]*100:.6g} \t {gas.mCO2[i+1]/gas.mCO[i+1]:.6g} \t {gas.mCO2[i+1]/gas.mH2O[i+1]:.6g} \t {gas.mCO2[i+1]/gas.mSO2[i+1]:.6g} \t {gas.mH2S[i+1]/gas.mSO2[i+1]:.6g} \t {gas.f['H2'][i]:.8g} \t {gas.f['H2O'][i]:.6g} \t {gas.f['CO2'][i]:.6g} \t {gas.f['CO'][i]:.6g} \t {gas.f['CH4'][i]:.6g} \t {gas.f['SO2'][i]:.8g} \t {gas.f['H2S'][i]:.8g} \t {gas.f['S2'][i]:.8g} \t {gas.f['N2'][i]:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'h', i)*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'c', i)*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'o', i)*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'o_tot', i)*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 's', i)*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'n', i)*1000000:.8g} \n")
                i += 1
        else:
            output.write(
                    f"{P[-1]} \t {cnvt.generate_fo2_buffer(sys, np.exp(gas.fo2[-1]), P[-1]):+.3} \t {np.exp(gas.fo2[-1]):.3e} \t {melt.F[-1]:.3} \t {sys.rho[-1]:.3} \t {melt.rho_store[-1]:.3} \t {sys.GvF[-1] * 100} \t {sys.WgT[-1] * 100:.5g} \t {gas.M[-1]} \t {gas.mH2O[-1]:.8g} \t {gas.mH2[-1]:.8g} \t {gas.mO2[-1]:.6g} \t {gas.mCO2[-1]:.6g} \t {gas.mCO[-1]:.6g} \t {gas.mCH4[-1]:.6g} \t {gas.mSO2[-1]:.6g} \t {gas.mH2S[-1]:.6g} \t {gas.mS2[-1]:.6g} \t {gas.mN2[-1]:.6g} \t {gas.Wt['H2O'][-1]:.6g} \t {gas.Wt['H2'][-1]:.6g} \t {gas.Wt['O2'][-1]:.6g} \t {gas.Wt['CO2'][-1]:.6g} \t {gas.Wt['CO'][-1]:.6g} \t {gas.Wt['CH4'][-1]:.6g} \t {gas.Wt['SO2'][-1]:.6g} \t {gas.Wt['H2S'][-1]:.6g} \t {gas.Wt['S2'][-1]:.6g} \t {gas.Wt['N2'][-1]:.6g} \t {melt.h2o[-1]*100:.6g} \t {melt.h2[-1]*100:.6g} \t {melt.co2[-1]*100:.6g} \t {melt.co[-1]*100:.6g} \t {melt.ch4[-1]*100:.6g} \t {melt.graphite[-1]*100:.6g} \t {melt.sulfide[-1]*100:.6g} \t {melt.sulfate[-1]*100:.6g} \t {melt.s[-1]*100:.6g} \t {melt.n[-1]*100:.6g} \t {gas.mCO2[-1]/gas.mCO[-1]:.6g} \t {gas.mCO2[-1]/gas.mH2O[-1]:.6g} \t {gas.mCO2[-1]/gas.mSO2[-1]:.6g} \t {gas.mH2S[-1]/gas.mSO2[-1]:.6g} \t {gas.f['H2'][-1]:.8g} \t {gas.f['H2O'][-1]:.6g} \t {gas.f['CO2'][-1]:.8g} \t {gas.f['CO'][-1]:.8g} \t {gas.f['CH4'][-1]:.6g} \t {gas.f['SO2'][-1]:.8g} \t {gas.f['H2S'][-1]:.8g} \t {gas.f['S2'][-1]:.8g} \t {gas.f['N2'][-1]:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'h', -1, WgT=sys.WgT[-1])*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'c', -1, WgT=sys.WgT[-1])*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'o', -1, WgT=sys.WgT[-1])*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'o_tot', -1, WgT=sys.WgT[-1])*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 's', -1, WgT=sys.WgT[-1])*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'n', -1, WgT=sys.WgT[-1])*1000000:.8g} \n")
        output.close()

def writeout_crash(sys, gas, melt, P):
    if not os.path.exists('Output'):
        os.makedirs('Output')
    
    output = open('Output/dgs_output.csv', 'w')
    output.write(f"#Decompressing a {sys.run.RUN_TYPE} {sys.run.COMPOSITION} {sys.run.GAS_SYS} system at {sys.T:.2f} degrees (K). \n#Crashed data file.\n")
    output.write("#\n")
    output.write(f"#{'':8} \t {'':8} \t {'':8} \t {'':8} \t {'':12} \t {'':12} \t {'Gas:':12} \t {'':10} \t {'':14} \t {'Gas by mol:':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'Gas by mass:':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'Melt (wt%):':8} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'Ratios:':8}\n")
    
    output.write(f"{'P'}\t{sys.run.FO2_buffer}\t{'fo2'}\t{'F'}\t{'rho_bulk'}\t{'rho_melt'}\t{'Exsol_vol%'}\t{'Gas_wt'}\t{'mol_mass'}\t{'mH2O'}\t{'mH2'}\t{'mO2'}\t{'mCO2'}\t{'mCO'}\t{'mCH4'}\t{'mSO2'}\t{'mH2S'}\t{'mS2'}\t{'mN2'}\t{'wH2O'}\t{'wH2'}\t{'wO2'}\t{'wCO2'}\t{'wCO'}\t{'wCH4'}\t{'wSO2'}\t{'wH2S'}\t{'wS2'}\t{'wN2'}\t{'H2O_melt'}\t{'H2_melt'}\t{'CO2_melt'}\t{'CO_melt'}\t{'CH4_melt'}\t{'graph_melt'}\t{'S2-_melt'}\t{'S6+_melt'}\t{'Stot_melt'}\t{'N_melt'}\t{'mCO2/CO'}\t{'mCO2/H2O'}\t{'mCO2/SO2'}\t{'mH2S/SO2'}\t{'fH2'}\t{'fH2O'}\t{'fCO2'}\t{'fCO'}\t{'fCH4'}\t{'fSO2'}\t{'fH2S'}\t{'fS2'}\t{'fN2'}\t{'tot_H'}\t{'tot_C'}\t{'tot_O_gas'}\t{'tot_O'}\t{'tot_S'}\t{'tot_N'}\n")

    i = 0

    if sys.run.FIND_SATURATION == True or sys.run.ATOMIC_MASS_SET == True:
        output.write(sat.satp_writeout(sys, melt, gas, sys.sat_conditions[0], sys.sat_conditions[1], sys.sat_conditions[2], sys.sat_conditions[3], graph_sat=sys.sat_conditions[4]))
    
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')    
        if len(P) > 1:
            while (i < len(P)):
                output.write(f"{P[i]} \t {cnvt.generate_fo2_buffer(sys, np.exp(gas.fo2[i]), P[i]):+.3} \t {np.exp(gas.fo2[i]):.3e} \t {melt.F[i+1]:.3} \t {sys.rho[i]:.3} \t {melt.rho_store[i]:.3} \t {sys.GvF[i] * 100:} \t {sys.WgT[i+1] * 100:.6g} \t {gas.M[i]} \t {gas.mH2O[i+1]:.6g} \t {gas.mH2[i+1]:.6g} \t {gas.mO2[i+1]:.6g} \t {gas.mCO2[i+1]:.6g} \t {gas.mCO[i+1]:.6g} \t {gas.mCH4[i+1]:.6g} \t {gas.mSO2[i+1]:.6g} \t {gas.mH2S[i+1]:.6g} \t {gas.mS2[i+1]:.6g} \t {gas.mN2[i+1]:.6g} \t {gas.Wt['H2O'][i]:.6g} \t {gas.Wt['H2'][i]:.6g} \t {gas.Wt['O2'][i]:.6g} \t {gas.Wt['CO2'][i]:.6g} \t {gas.Wt['CO'][i]:.6g} \t {gas.Wt['CH4'][i]:.6g} \t {gas.Wt['SO2'][i]:.6g} \t {gas.Wt['H2S'][i]:.6g} \t {gas.Wt['S2'][i]:.6g} \t {gas.Wt['N2'][i]:.6g} \t {melt.h2o[i]*100:.6g} \t {melt.h2[i]*100:.6g} \t {melt.co2[i]*100:.6g} \t {melt.co[i]*100:.6g} \t {melt.ch4[i]*100:.6g} \t {melt.graphite[i]*100:.6g} \t {melt.sulfide[i]*100:.6g} \t {melt.sulfate[i]*100:.6g} \t {melt.s[i]*100:.6g} \t {melt.n[i]*100:.6g} \t {gas.mCO2[i+1]/gas.mCO[i+1]:.6g} \t {gas.mCO2[i+1]/gas.mH2O[i+1]:.6g} \t {gas.mCO2[i+1]/gas.mSO2[i+1]:.6g} \t {gas.mH2S[i+1]/gas.mSO2[i+1]:.6e} \t {gas.f['H2'][i]:.8g} \t {gas.f['H2O'][i]:.6g} \t {gas.f['CO2'][i]:.6g} \t {gas.f['CO'][i]:.6g} \t {gas.f['CH4'][i]:.6g} \t {gas.f['SO2'][i]:.8g} \t {gas.f['H2S'][i]:.8g} \t {gas.f['S2'][i]:.8g} \t {gas.f['N2'][i]:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'h', i)*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'c', i)*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'o', i)*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'o_tot', i)*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 's', i)*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'n', i)*1000000:.8g} \n")
                i += 1
        else:
            output.write(
                f"{P[- 1]} \t {cnvt.generate_fo2_buffer(sys, np.exp(gas.fo2[-1]), P[-1]):+.3} \t {np.exp(gas.fo2[-1]):.3e} \t {melt.F[-1]:.3} \t {sys.rho[-1]:.3} \t {melt.rho_store[-1]:.3} \t {sys.GvF[-1] * 100} \t {sys.WgT[-1] * 100:.6g} \t {gas.M[-1]} \t {gas.mH2O[-1]:.6g} \t {gas.mH2[-1]:.6g} \t {gas.mO2[-1]:.6g} \t {gas.mCO2[-1]:.6g} \t {gas.mCO[-1]:.6g} \t {gas.mCH4[-1]:.6g} \t {gas.mSO2[-1]:.6g} \t {gas.mH2S[-1]:.6g} \t {gas.mS2[-1]:.6g} \t {gas.mN2[-1]:.6g} \t {gas.Wt['H2O'][-1]:.6g} \t {gas.Wt['H2'][-1]:.6g} \t {gas.Wt['O2'][-1]:.6g} \t {gas.Wt['CO2'][-1]:.6g} \t {gas.Wt['CO'][-1]:.6g} \t {gas.Wt['CH4'][-1]:.6g} \t {gas.Wt['SO2'][-1]:.6g} \t {gas.Wt['H2S'][-1]:.6g} \t {gas.Wt['S2'][-1]:.6g} \t {gas.Wt['N2'][-1]:.6g} \t {melt.h2o[-1]*100:.6g} \t {melt.h2[-1]*100:.6g} \t {melt.co2[-1]*100:.6g} \t {melt.co[-1]*100:.6g} \t {melt.ch4[-1]*100:.6g} \t {melt.graphite[-1]*100:.6g} \t {melt.sulfide[-1]*100:.6g} \t {melt.sulfate[-1]*100:.6g} \t {melt.s[-1]*100:.6g} \t {melt.n[-1]*100:.6g} \t {gas.mCO2[-1]/gas.mCO[-1]:.6g} \t {gas.mCO2[-1]/gas.mH2O[-1]:.6g} \t {gas.mCO2[-1]/gas.mSO2[-1]:.6g} \t {gas.mH2S[-1]/gas.mSO2[-1]:.6g} \t {gas.f['H2'][-1]:.8g} \t {gas.f['H2O'][-1]:.6g} \t {gas.f['CO2'][-1]:.8g} \t {gas.f['CO'][-1]:.8g} \t {gas.f['CH4'][-1]:.6g} \t {gas.f['SO2'][-1]:.8g} \t {gas.f['H2S'][-1]:.8g} \t {gas.f['S2'][-1]:.8g} \t {gas.f['N2'][-1]:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'h', -1, WgT=sys.WgT[-1])*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'c', -1, WgT=sys.WgT[-1])*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'o', -1, WgT=sys.WgT[-1])*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'o_tot', -1, WgT=sys.WgT[-1])*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 's', -1, WgT=sys.WgT[-1])*1000000:.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'n', -1, WgT=sys.WgT[-1])*1000000:.8g} \n")
    
    output.close()

def writeout_ocs(sys, gas, melt, P):
    output = open('Output/dgs_output.csv', 'w')
    output.write(f"#Decompressing a {sys.run.COMPOSITION} {sys.run.GAS_SYS} system at {sys.T:.2f} degrees (K) with a final conversion to OCS production. \n"
                 f"#Approx total H2O (wt%) starting fO2: {(gas.Wt['H2O'][0]*sys.WgT[0] + (melt.h2o[0]/100))*100:.3}, final fO2: {(gas.Wt['H2O'][-1]*sys.WgT[-1] + (melt.h2o[-1]/100))*100:.3}\n"
                 f"#Approx total CO2 starting fO2: {(gas.Wt['CO2'][0]*sys.WgT[0] + (melt.co2[0]/100))*100:.3} (wt%) {(gas.Wt['CO2'][0]*sys.WgT[0] + (melt.co2[0]/100))*1000000:.0} (ppm), final fO2: {((gas.Wt['CO2'][-1]*sys.WgT[-1] + (melt.co2[-1]/100))*100)*10000:.0} (ppm) \n"
                 f"#Approx total C (ppm) starting fO2: {(((gas.Wt['CO2'][0]*sys.WgT[1]) + (melt.co2[0]/100))/cnst.m['co2'] + (gas.Wt['CO'][0]*sys.WgT[1])/cnst.m['co'] + (gas.Wt['CH4'][0]*sys.WgT[1])/cnst.m['ch4'])*1000000*cnst.m['c']}, final fO2: {(((gas.Wt['CO2'][-1]*sys.WgT[-1]) + (melt.co2[-1]/100))/cnst.m['co2'] + (gas.Wt['CO'][-1]*sys.WgT[-1])/cnst.m['co'] + (gas.Wt['CH4'][-1]*sys.WgT[-1])/cnst.m['ch4'])*1000000*cnst.m['c']}\n"
                 f"#Total H (ppm): {sys.atomicM['h']*1000000} wt%: {sys.atomicM['h']*100}\n")

    output.write("#\n")
    output.write(f"#{'':8} \t {'':8} \t {'':8} \t {'':12} \t {'':12} \t {'Gas:':12} \t {'':10} \t {'by mol:':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'by mass:':10} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'Melt (wt%):':8} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'Ratios:':8} \t {'':8} \t {'':8} \t {'':8} \n")
    
    output.write(f"#{'P (bar)':8} \t {'FMQ':8} \t {'fo2':8} \t {'rho (bulk):':12} \t {'rho (melt):':12} \t {'Exsol. vol%':12} \t {'Gas wt%':10} \t {'mH2O':12} \t {'mH2':12} \t {'mO2':12} \t {'mCO2':12} \t {'mCO':12} \t {'mCH4':12} \t {'mSO2':12} \t {'mH2S':12} \t {'mS2':12} \t {'mOCS':12} \t {'wH2O':10} \t {'wH2':12} \t {'wO2':12} \t {'wCO2':8} \t {'wCO':12} \t {'wCH4':12} \t {'wSO2':12} \t {'wH2S':12} \t {'wS2':12} \t {'H2Omelt:':10} \t {'H2melt:':12} \t {'CO2melt:':10} \t {'S2-_melt:':10} \t {'S6+_melt:':10} \t {'mCO2/CO':12} \t {'mCO2/H2O':12} \t {'mCO2/SO2':12} \t {'mH2S/SO2':12} \t {'fH2':10} \t {'fH2O':10} \t {'fCO2':10} \t {'fCO':10} \t {'fCH4':10} \t {'fSO2':10}\t {'fH2S':10} \t {'fS2':10} \t {'fN2':10} \t {'tot H (ppm)':12} \t {'tot C (ppm)':12} \t {'tot O (ppm)':12} \t {'tot S (ppm)':12} \t {'tot N (ppm)':12} \n")

    i = 0

    if sys.run.SINGLE_STEP == False:

        while i < len(P):
            output.write(f"{P[i]:8.1f} \t {cnvt.fo2_2fmq(np.log10(np.exp(gas.fo2[i])),sys.T,P[i],sys.run.FMQ_MODEL):+8.3} \t {np.exp(gas.fo2[i]):8.3e} \t {sys.rho[i]:12.3} \t {melt.rho(P=(P[i])*1e5):12.3} \t {sys.GvF[i] * 100:12} \t {sys.WgT[i+1] * 100:10} \t {gas.mH2O[i+1]:12.g} \t {gas.mH2[i+1]:12.g} \t {gas.mO2[i+1]:12.g} \t {gas.mCO2[i+1]:12.g} \t {gas.mCO[i+1]:12.g} \t {gas.mCH4[i+1]:12.g} \t {gas.mSO2[i+1]:12.g} \t {gas.mH2S[i+1]:12.g} \t {gas.mS2[i+1]:12.g} \t {0:12} \t {gas.Wt['H2O'][i]:10.g} \t {gas.Wt['H2'][i]:12.g} \t {gas.Wt['O2'][i]:12.g} \t {gas.Wt['CO2'][i]:8.g} \t {gas.Wt['CO'][i]:12.g} \t {gas.Wt['CH4'][i]:12.g} \t {gas.Wt['SO2'][i]:12.g} \t {gas.Wt['H2S'][i]:12.g} \t {gas.Wt['S2'][i]:12.g} \t {melt.h2o[i]:10.g} \t {melt.h2[i]:12.g} \t {melt.co2[i]:10.g} \t {melt.sulfide[i]:10.g} \t {melt.sulfate[i]:10.g} \t {gas.mCO2[i+1]/gas.mCO[i+1]:8.6e} \t {gas.mCO2[i+1]/gas.mH2O[i+1]:8.6e} \t {gas.mCO2[i+1]/gas.mSO2[i+1]:8.6e} \t {gas.mH2S[i+1]/gas.mSO2[i+1]:8.6e} \t {gas.f['H2'][i]:10.8g} \t {gas.f['H2O'][i]:10.6g} \t {gas.f['CO2'][i]:10.6e} \t {gas.f['CO'][i]:10.6e} \t {gas.f['CH4'][i]:10.6g} \t {gas.f['SO2'][i]:10.8g} \t {gas.f['H2S'][i]:10.8g} \t {gas.f['S2'][i]:10.8g} \t {gas.f['N2'][i]:10.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'h', i)*1000000:12.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'c', i)*1000000:12.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'o', i)*1000000:12.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 's', i)*1000000:12.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'n', i)*1000000:12.8g} \n")
            i += 1
            
    else:
        output.write(
                f"{P[- 1]:8.1f} \t {cnvt.fo2_2fmq(np.log10(np.exp(gas.fo2[-2])), sys.T, P[-1], sys.run.FMQ_MODEL):+8.3} \t {np.exp(gas.fo2[-2]):8.3e} \t {sys.rho[-2]:12.3} \t {melt.rho(P=(P[-1])*1e5):12.3} \t {sys.GvF[-2] * 100:12} \t {sys.WgT[-2] * 100:10} \t {gas.mH2O[-2]:12.g} \t {gas.mH2[-2]:12.g} \t {gas.mO2[-2]:12.g} \t {gas.mCO2[-2]:12.g} \t {gas.mCO[-2]:12.g} \t {gas.mCH4[-2]:12.g} \t {gas.mSO2[-2]:12.g} \t {gas.mH2S[-2]:12.g} \t {gas.mS2[-2]:12.g} \t {0:12} \t {gas.Wt['H2O'][-1]:10.g} \t {gas.Wt['H2'][-1]:12.g} \t {gas.Wt['O2'][-1]:12.g} \t {gas.Wt['CO2'][-1]:8.g} \t {gas.Wt['CO'][-1]:12.g} \t {gas.Wt['CH4'][-1]:12.g} \t {gas.Wt['SO2'][-1]:12.g} \t {gas.Wt['H2S'][-1]:12.g} \t {gas.Wt['S2'][-1]:12.g} \t {melt.h2o[-2]:10.g} \t {melt.h2[-2]:12.g} \t {melt.co2[-2]:10.g} \t {melt.sulfide[-2]:10.g} \t {melt.sulfate[-2]:10.g} \t {gas.mCO2[-2]/gas.mCO[-2]:8.g} \t {gas.mCO2[-2]/gas.mH2O[-2]:8.6g} \t {gas.mCO2[-2]/gas.mSO2[-2]:8.6g} \t {gas.mH2S[-2]/gas.mSO2[-2]:8.6g} \t {gas.f['H2'][-2]:10.8g} \t {gas.f['H2O'][-2]:10.6g} \t {gas.f['CO2'][-2]:10.8g} \t {gas.f['CO'][-2]:10.8g} \t {gas.f['CH4'][-2]:10.6g} \t {gas.f['SO2'][-2]:10.8g} \t {gas.f['H2S'][-2]:10.8g} \t {gas.f['S2'][-2]:10.8g} \t {gas.f['N2'][-2]:10.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'h', -1)*1000000:12.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'c', -1)*1000000:12.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 'o', -1)*1000000:12.8g} \t {cnvt.atomicM_calc(sys, melt, gas, 's', -1)*1000000:12.8g}  \t {cnvt.atomicM_calc(sys, melt, gas, 'n', -1)*1000000:12.8g} \n")
    
    # OCS line last
    output.write(f"{P[- 1]:8.1f} \t {cnvt.fo2_2fmq(np.log10(np.exp(gas.fo2[-1])), sys.T, P[-1], sys.run.FMQ_MODEL):+8.3} \t {np.exp(gas.fo2[-1]):8.3e} \t {sys.rho[-1]:12.3} \t {melt.rho(P=(P[-1])*1e5):12.3} \t {sys.GvF[-1] * 100:12} \t {sys.WgT[-1] * 100:10} \t {gas.mH2O[-1]:12.g} \t {gas.mH2[-1]:12.g} \t {gas.mO2[-1]:12.g} \t {gas.mCO2[-1]:12.g} \t {gas.mCO[-1]:12.g} \t {gas.mCH4[-1]:12.g} \t {gas.mSO2[-1]:12.g} \t {gas.mH2S[-1]:12.g} \t {gas.mS2[-1]:12.g} \t {gas.mOCS[-1]:12.g} \t {'N/A':10} \t {'N/A':12} \t {'N/A':12} \t {'N/A':8} \t {'N/A':12} \t {'N/A':12} \t {'N/A':12} \t {'N/A':12} \t {'N/A':12} \t {melt.h2o[-1]:10.g} \t {melt.h2[-1]:12.g} \t {melt.co2[-1]:10.g} \t {melt.sulfide[-1]:10.g} \t {melt.sulfate[-1]:10.g} \t {gas.mCO2[-1]/gas.mCO[-1]:8.6e} \t {gas.mCO2[-1]/gas.mH2O[-1]:8.6e} \t {gas.mCO2[-1]/gas.mSO2[-1]:8.6e} \t {gas.mH2S[-1]/gas.mSO2[-1]:8.6e} \t {gas.f['H2'][-1]:10.8g} \t {gas.f['H2O'][-1]:10.6g} \t {gas.f['CO2'][-1]:10.8g} \t {gas.f['CO'][-1]:10.8g} \t {gas.f['CH4'][-1]:10.6g} \t {gas.f['SO2'][-1]:10.8g} \t {gas.f['H2S'][-1]:10.8g} \t {gas.f['S2'][-1]:10.8g} \t {gas.f['N2'][-1]:10.8g} \t {'N/A':12} \t {'N/A':12} \t {'N/A':12} \t {'N/A':12} \t {'N/A':12} \n")
    
    output.close()

def writeout_figs(sys, melt, gas, out, P):

    filelist = glob.glob("Output/*.png")
    for file in filelist:
        os.remove(file)  # Removes previous files so if output specification is changed there is no confusion as to up to date files.

    if out is not None:  # If an output file listing requested figures has been included:
        if out.Plt_melt_species:
            plot_meltspecies(melt, P)

        if out.Plt_gas_species_wt:
            plot_gasspecies_wt(gas, P)

        if out.Plt_gas_species_mol:
            plot_gasspecies_mol(gas, P)

        if out.Plt_gas_fraction:
            plot_gasfraction(sys, P)

        if out.Plt_fo2_dFMQ:
            plot_fo2FMQ(melt, gas, P)

        if out.Plt_gas_fraction != None:
            pass
            #plot

    else:  # else plot every option
        plot_meltspecies(melt, P)
        plot_gasspecies_wt(gas, P)
        plot_gasspecies_mol(gas, P)
        plot_gasfraction(sys, P)
        plot_fo2FMQ(melt, gas, P)


# fO2 and dFMQ
def plot_fo2FMQ(melt, gas, P):
    plt.plot(P, melt.fmq)
    plt.xscale('log')
    plt.xlabel('Pressure (bars)')
    plt.ylabel(r'$\Delta$ FMQ')
    plt.savefig('Output/FMQ.png')
    plt.close()

    fo2 = []
    for x in gas.fo2:
        fo2.append(np.log10(np.exp(x)))

    plt.plot(P, fo2)
    #plt.xscale('log')
    plt.xlabel('Pressure (bars)')
    plt.ylabel('log(10) fO2')
    plt.savefig('Output/fO2.png')
    plt.close()

# Gas speciation mol
def plot_gasspecies_mol(gas, P):

    plt.plot(gas.mH2O[1:], P, c='darkblue')
    plt.plot(gas.mH2[1:], P, c='steelblue')
    plt.plot(gas.mO2[1:], P, c='firebrick')
    plt.annotate('H2O', xy=[gas.mH2O[-1], P[-1]], color='darkblue')
    plt.annotate('H2', xy=[gas.mH2[-1], P[-1]], color='steelblue')
    plt.annotate('O2', xy=[gas.mO2[-1], P[-1]], color='firebrick')
    if gas.sys.run.GAS_SYS == 'COH' or gas.sys.run.GAS_SYS == 'COHS' or gas.sys.run.GAS_SYS == 'COHSN':
        plt.plot(gas.mCO2[1:], P, c='saddlebrown')
        plt.plot(gas.mCO[1:], P, c='darkgreen')
        plt.plot(gas.mCH4[1:], P, c='forestgreen')
        plt.annotate('CO2', xy=[gas.mCO2[-1], P[-1]], color='saddlebrown')
        plt.annotate('CO', xy=[gas.mCO[-1], P[-1]], color='darkgreen')
        plt.annotate('CH4', xy=[gas.mCH4[-1], P[-1]], color='forestgreen')
    if gas.sys.run.GAS_SYS == 'SOH' or gas.sys.run.GAS_SYS == 'COHS' or gas.sys.run.GAS_SYS == 'COHSN':
        plt.plot(gas.mSO2[1:], P, c='maroon')
        plt.plot(gas.mS2[1:], P, c='goldenrod')
        plt.plot(gas.mH2S[1:], P, c='pink')
        plt.annotate('SO2', xy=[gas.mSO2[-1], P[-1]], color='maroon')
        plt.annotate('S2', xy=[gas.mS2[-1], P[-1]], color='goldenrod')
        plt.annotate('H2S', xy=[gas.mH2S[-1], P[-1]], color='pink')
    if gas.sys.run.GAS_SYS == 'COHSN':
        plt.plot(gas.mN2[1:], P, c='grey')
        plt.annotate('N2', xy=[gas.mN2[-1], P[-1]], color='grey')
    plt.xscale('log')
    plt.gca().invert_yaxis()
    plt.xlabel(f'Speciation in a {gas.sys.run.GAS_SYS} gas (mol frac)')
    plt.ylabel('Pressure (bar)')
    plt.savefig('Output/speciation(mol).png')
    plt.close()

# Gas speciation wt%

def plot_gasspecies_wt(gas, P):

    plt.plot(gas.Wt['H2O'], P, c='darkblue')
    plt.plot(gas.Wt['H2'], P, c='steelblue')
    plt.plot(gas.Wt['O2'], P, c='firebrick')
    plt.annotate('H2O', xy=[gas.Wt['H2O'][-1], P[-1]], color='darkblue')
    plt.annotate('H2', xy=[gas.Wt['H2'][-1], P[-1]], color='steelblue')
    plt.annotate('O2', xy=[gas.Wt['O2'][-1], P[-1]], color='firebrick')
    if gas.sys.run.GAS_SYS == 'COH' or gas.sys.run.GAS_SYS == 'COHS' or gas.sys.run.GAS_SYS == 'COHSN':
        plt.plot(gas.Wt['CO2'], P, c='saddlebrown')
        plt.plot(gas.Wt['CO'], P, c='darkgreen')
        #plt.plot(gas.Wt['CH4'], P, c='mediumseagreen')
        plt.annotate('CO2', xy=[gas.Wt['CO2'][-1], P[-1]], color='saddlebrown')
        plt.annotate('CO', xy=[gas.Wt['CO'][-1], P[-1]], color='darkgreen')
        #plt.annotate('CH4', xy=[gas.mCH4[-1], 0.9], color='mediumseagreen')
    if gas.sys.run.GAS_SYS == 'SOH' or gas.sys.run.GAS_SYS == 'COHS' or gas.sys.run.GAS_SYS == 'COHSN':
        plt.plot(gas.Wt['SO2'], P, c='maroon')
        plt.plot(gas.Wt['S2'], P, c='goldenrod')
        plt.plot(gas.Wt['H2S'], P, c='pink')
        plt.annotate('SO2', xy=[gas.Wt['SO2'][-1], P[-1]], color='maroon')
        plt.annotate('S2', xy=[gas.Wt['S2'][-1], P[-1]], color='goldenrod')
        plt.annotate('H2S', xy=[gas.Wt['H2S'][-1], P[-1]], color='pink')
    if gas.sys.run.GAS_SYS == 'COHSN':
        plt.plot(gas.Wt['N2'], P, c='grey')
        plt.annotate('N2', xy=[gas.Wt['N2'][-1], P[-1]], color='grey')
    plt.xscale('log')
    plt.gca().invert_yaxis()
    plt.xlabel(f'Gas phase speciation of a {gas.sys.run.GAS_SYS} system (wt %)')
    plt.ylabel('Pressure (bar)')
    plt.savefig('Output/speciation(wt).png')
    plt.close()

# Melt volatile wt%

def plot_meltspecies(melt, P):

    plt.plot(melt.h2o, P, c='darkblue')
    plt.plot(melt.h2, P, c='steelblue')
    plt.annotate('H2O', xy=[melt.h2o[-1], P[-1]], color='darkblue')
    plt.annotate('H2', xy=[melt.h2[-1], P[-1]], color='steelblue')
    if melt.sys.run.GAS_SYS == 'COH' or melt.sys.run.GAS_SYS == 'COHS' or melt.sys.run.GAS_SYS == 'COHSN':
        plt.plot(melt.co2, P, c='saddlebrown')
        plt.annotate('CO2', xy=[melt.co2[-1], P[-1]], color='saddlebrown')
    if melt.sys.run.GAS_SYS == 'SOH' or melt.sys.run.GAS_SYS == 'COHS' or melt.sys.run.GAS_SYS == 'COHSN':
        plt.plot(melt.s, P, c='maroon')
        plt.annotate('S', xy=[melt.s[-1], P[-1]], color='maroon')
    if melt.sys.run.GAS_SYS == 'COHSN':
        plt.plot(melt.n, P, c='grey')
        plt.annotate('N', xy=[melt.n[-1], P[-1]], color='grey')
    plt.xscale('log')
    plt.gca().invert_yaxis()
    plt.xlabel('Melt volatile content (wt%)')
    plt.ylabel('Pressure (bar)')
    plt.savefig('Output/meltspecies.png')
    plt.close()

# Exsolved gas mass and volume

def plot_gasfraction(sys, P):

    fig, (ax1, ax2) = plt.subplots(1,2, sharey=True)
    ax1.plot(cnvt.frac2perc(sys.WgT[1:]), P)
    ax1.invert_yaxis()
    ax1.set_xlabel('Exsolved gas wt%')
    ax1.set_ylabel('Pressure (bar)')
    ax2.plot(cnvt.frac2perc(sys.GvF), P)
    ax2.set_xlabel('Exsolved gas volume %')
    fig.savefig('Output/exsolved_gas.png')



