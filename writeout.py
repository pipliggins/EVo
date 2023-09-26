"""
Writes the results of the run to a file and
contains options to produce a graph of the results.
"""
import conversions as cnvt
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sat_pressure as sat
import warnings


def writeout_file(sys, gas, melt, P):
    """
    At the end of a run, writes data to a single output CSV file

    Parameters
    ----------
    sys : ThermoSystem class
        Active instance of the ThermoSystem class
    gas : Gas class
        Active instance of the Gas class
    melt : Melt class
        Active instance of the Melt class
    P : float
        Final pressure (bar)
    """

    if not os.path.exists("Output"):
        os.makedirs("Output")

    output = open("Output/dgs_output.csv", "w")
    # write file header
    output.write(
        f"#Decompressing a {sys.run.COMPOSITION} {sys.run.GAS_SYS} {sys.run.RUN_TYPE} "
        f"system at {sys.T:.2f} degrees (K). \n"
    )

    output.write("#\n")
    output.write(
        f"#{'':8} \t {'':8} \t {'':8} \t {'':8} \t {'':12} \t {'':12} \t {'Gas:':12} \t"
        f" {'':10} \t {'':14} \t {'Gas by mol:':12} \t {'':12} \t {'':12} \t {'':12} \t"
        f" {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t "
        f"{'Gas by mass:':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t"
        f" {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'Melt (wt%):':8} \t {'':10} \t "
        f"{'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t "
        f"{'':10} \t {'Ratios:':8}\n"
    )
    # write column headings
    output.write(
        f"{'P'}\t{sys.run.FO2_buffer}\t{'fo2'}\t{'F'}\t{'rho_bulk'}\t{'rho_melt'}\t"
        f"{'Exsol_vol%'}\t{'Gas_wt'}\t{'mol_mass'}\t{'mH2O'}\t{'mH2'}\t{'mO2'}\t"
        f"{'mCO2'}\t{'mCO'}\t{'mCH4'}\t{'mSO2'}\t{'mH2S'}\t{'mS2'}\t{'mN2'}\t{'wH2O'}\t"
        f"{'wH2'}\t{'wO2'}\t{'wCO2'}\t{'wCO'}\t{'wCH4'}\t{'wSO2'}\t{'wH2S'}\t{'wS2'}\t"
        f"{'wN2'}\t{'H2O_melt'}\t{'H2_melt'}\t{'CO2_melt'}\t{'CO_melt'}\t{'CH4_melt'}\t"
        f"{'graph_melt'}\t{'S2-_melt'}\t{'S6+_melt'}\t{'Stot_melt'}\t{'N_melt'}\t"
        f"{'mCO2/CO'}\t{'mCO2/H2O'}\t{'mCO2/SO2'}\t{'mH2S/SO2'}\t{'fH2'}\t{'fH2O'}\t"
        f"{'fCO2'}\t{'fCO'}\t{'fCH4'}\t{'fSO2'}\t{'fH2S'}\t{'fS2'}\t{'fN2'}\t{'tot_H'}"
        f"\t{'tot_C'}\t{'tot_O_gas'}\t{'tot_O'}\t{'tot_S'}\t{'tot_N'}\n"
    )

    i = 0

    if sys.run.FIND_SATURATION is True or sys.run.ATOMIC_MASS_SET is True:
        output.write(
            sat.satp_writeout(
                sys,
                melt,
                gas,
                sys.sat_conditions[0],
                sys.sat_conditions[1],
                sys.sat_conditions[2],
                sys.sat_conditions[3],
                graph_sat=sys.sat_conditions[4],
            )
        )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        if sys.run.SINGLE_STEP is False:
            while i < len(P):
                output.write(
                    f"{P[i]} \t "
                    f"{cnvt.generate_fo2_buffer(sys, np.exp(gas.fo2[i]), P[i]):+.3} \t "
                    f"{np.exp(gas.fo2[i]):.3e} \t {melt.F[i+1]:.3} \t {sys.rho[i]:.3} "
                    f"\t {melt.rho_store[i]:.3} \t {sys.GvF[i] * 100:} \t "
                    f"{sys.WgT[i+1] * 100:.5g} \t {gas.M[i]} \t {gas.mH2O[i+1]:.8g} \t "
                    f"{gas.mH2[i+1]:.8g} \t {gas.mO2[i+1]:.6g} \t {gas.mCO2[i+1]:.6g} "
                    f"\t {gas.mCO[i+1]:.6g} \t {gas.mCH4[i+1]:.6g} \t "
                    f"{gas.mSO2[i+1]:.6g} \t {gas.mH2S[i+1]:.6g} \t {gas.mS2[i+1]:.6g} "
                    f"\t {gas.mN2[i+1]:.6g} \t {gas.Wt['H2O'][i]:.6g} \t "
                    f"{gas.Wt['H2'][i]:.6g} \t {gas.Wt['O2'][i]:.6g} \t "
                    f"{gas.Wt['CO2'][i]:.6g} \t {gas.Wt['CO'][i]:.6g} \t "
                    f"{gas.Wt['CH4'][i]:.6g} \t {gas.Wt['SO2'][i]:.6g} \t "
                    f"{gas.Wt['H2S'][i]:.6g} \t {gas.Wt['S2'][i]:.6g} \t "
                    f"{gas.Wt['N2'][i]:.6g} \t {melt.h2o[i]*100:.6g} \t "
                    f"{melt.h2[i]*100:.6g} \t {melt.co2[i]*100:.6g} \t "
                    f"{melt.co[i]*100:.6g} \t {melt.ch4[i]*100:.6g} \t "
                    f"{melt.graphite[i]*100:.6g} \t {melt.sulfide[i]*100:.6g} \t "
                    f"{melt.sulfate[i]*100:.6g} \t {melt.s[i]*100:.6g} \t "
                    f"{melt.n[i]*100:.6g} \t {gas.mCO2[i+1]/gas.mCO[i+1]:.6g} \t "
                    f"{gas.mCO2[i+1]/gas.mH2O[i+1]:.6g} \t "
                    f"{gas.mCO2[i+1]/gas.mSO2[i+1]:.6g} \t "
                    f"{gas.mH2S[i+1]/gas.mSO2[i+1]:.6g} \t {gas.f['H2'][i]:.8g} \t "
                    f"{gas.f['H2O'][i]:.6g} \t {gas.f['CO2'][i]:.6g} \t "
                    f"{gas.f['CO'][i]:.6g} \t {gas.f['CH4'][i]:.6g} \t "
                    f"{gas.f['SO2'][i]:.8g} \t {gas.f['H2S'][i]:.8g} \t "
                    f"{gas.f['S2'][i]:.8g} \t {gas.f['N2'][i]:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 'h', i)*1000000:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 'c', i)*1000000:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 'o', i)*1000000:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 'o_tot', i)*1000000:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 's', i)*1000000:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 'n', i)*1000000:.8g} \n"
                )
                i += 1
        else:
            output.write(
                f"{P[-1]} \t "
                f"{cnvt.generate_fo2_buffer(sys, np.exp(gas.fo2[-1]), P[-1]):+.3} \t "
                f"{np.exp(gas.fo2[-1]):.3e} \t {melt.F[-1]:.3} \t {sys.rho[-1]:.3} \t "
                f"{melt.rho_store[-1]:.3} \t {sys.GvF[-1] * 100} \t "
                f"{sys.WgT[-1] * 100:.5g} \t {gas.M[-1]} \t {gas.mH2O[-1]:.8g} \t "
                f"{gas.mH2[-1]:.8g} \t {gas.mO2[-1]:.6g} \t {gas.mCO2[-1]:.6g} \t "
                f"{gas.mCO[-1]:.6g} \t {gas.mCH4[-1]:.6g} \t {gas.mSO2[-1]:.6g} \t "
                f"{gas.mH2S[-1]:.6g} \t {gas.mS2[-1]:.6g} \t {gas.mN2[-1]:.6g} \t "
                f"{gas.Wt['H2O'][-1]:.6g} \t {gas.Wt['H2'][-1]:.6g} \t "
                f"{gas.Wt['O2'][-1]:.6g} \t {gas.Wt['CO2'][-1]:.6g} \t "
                f"{gas.Wt['CO'][-1]:.6g} \t {gas.Wt['CH4'][-1]:.6g} \t "
                f"{gas.Wt['SO2'][-1]:.6g} \t {gas.Wt['H2S'][-1]:.6g} \t "
                f"{gas.Wt['S2'][-1]:.6g} \t {gas.Wt['N2'][-1]:.6g} \t "
                f"{melt.h2o[-1]*100:.6g} \t {melt.h2[-1]*100:.6g} \t "
                f"{melt.co2[-1]*100:.6g} \t {melt.co[-1]*100:.6g} \t "
                f"{melt.ch4[-1]*100:.6g} \t {melt.graphite[-1]*100:.6g} \t "
                f"{melt.sulfide[-1]*100:.6g} \t {melt.sulfate[-1]*100:.6g} \t "
                f"{melt.s[-1]*100:.6g} \t {melt.n[-1]*100:.6g} \t "
                f"{gas.mCO2[-1]/gas.mCO[-1]:.6g} \t {gas.mCO2[-1]/gas.mH2O[-1]:.6g} \t "
                f"{gas.mCO2[-1]/gas.mSO2[-1]:.6g} \t {gas.mH2S[-1]/gas.mSO2[-1]:.6g} \t"
                f" {gas.f['H2'][-1]:.8g} \t {gas.f['H2O'][-1]:.6g} \t "
                f"{gas.f['CO2'][-1]:.8g} \t {gas.f['CO'][-1]:.8g} \t "
                f"{gas.f['CH4'][-1]:.6g} \t {gas.f['SO2'][-1]:.8g} \t "
                f"{gas.f['H2S'][-1]:.8g} \t {gas.f['S2'][-1]:.8g} \t "
                f"{gas.f['N2'][-1]:.8g} \t "
                f"{cnvt.atomicM_calc(sys, melt, gas, 'h', -1, WgT=sys.WgT[-1])*1000000:.8g} \t "  # noqa:E501
                f"{cnvt.atomicM_calc(sys, melt, gas, 'c', -1, WgT=sys.WgT[-1])*1000000:.8g} \t "  # noqa:E501
                f"{cnvt.atomicM_calc(sys, melt, gas, 'o', -1, WgT=sys.WgT[-1])*1000000:.8g} \t "  # noqa:E501
                f"{cnvt.atomicM_calc(sys, melt, gas, 'o_tot', -1, WgT=sys.WgT[-1])*1000000:.8g} \t "  # noqa:E501
                f"{cnvt.atomicM_calc(sys, melt, gas, 's', -1, WgT=sys.WgT[-1])*1000000:.8g} \t "  # noqa:E501
                f"{cnvt.atomicM_calc(sys, melt, gas, 'n', -1, WgT=sys.WgT[-1])*1000000:.8g} \n"  # noqa:E501
            )
        output.close()


def writeout_crash(sys, gas, melt, P):
    """
    For a run that's ended early, writes data to a single output CSV file

    Parameters
    ----------
    sys : ThermoSystem class
        Active instance of the ThermoSystem class
    gas : Gas class
        Active instance of the Gas class
    melt : Melt class
        Active instance of the Melt class
    P : float
        Current pressure (bar)
    """

    if not os.path.exists("Output"):
        os.makedirs("Output")

    output = open("Output/dgs_output.csv", "w")
    output.write(
        f"#Decompressing a {sys.run.RUN_TYPE} {sys.run.COMPOSITION} {sys.run.GAS_SYS} "
        f"system at {sys.T:.2f} degrees (K). \n#Crashed data file.\n"
    )
    output.write("#\n")
    output.write(
        f"#{'':8} \t {'':8} \t {'':8} \t {'':8} \t {'':12} \t {'':12} \t {'Gas:':12} \t"
        f" {'':10} \t {'':14} \t {'Gas by mol:':12} \t {'':12} \t {'':12} \t {'':12} \t"
        f" {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t "
        f"{'Gas by mass:':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'':12} \t"
        f" {'':12} \t {'':12} \t {'':12} \t {'':12} \t {'Melt (wt%):':8} \t {'':10} \t "
        f"{'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t {'':10} \t "
        f"{'':10} \t {'Ratios:':8}\n"
    )
    # write column headings
    output.write(
        f"{'P'}\t{sys.run.FO2_buffer}\t{'fo2'}\t{'F'}\t{'rho_bulk'}\t{'rho_melt'}\t"
        f"{'Exsol_vol%'}\t{'Gas_wt'}\t{'mol_mass'}\t{'mH2O'}\t{'mH2'}\t{'mO2'}\t"
        f"{'mCO2'}\t{'mCO'}\t{'mCH4'}\t{'mSO2'}\t{'mH2S'}\t{'mS2'}\t{'mN2'}\t{'wH2O'}\t"
        f"{'wH2'}\t{'wO2'}\t{'wCO2'}\t{'wCO'}\t{'wCH4'}\t{'wSO2'}\t{'wH2S'}\t{'wS2'}\t"
        f"{'wN2'}\t{'H2O_melt'}\t{'H2_melt'}\t{'CO2_melt'}\t{'CO_melt'}\t{'CH4_melt'}\t"
        f"{'graph_melt'}\t{'S2-_melt'}\t{'S6+_melt'}\t{'Stot_melt'}\t{'N_melt'}\t"
        f"{'mCO2/CO'}\t{'mCO2/H2O'}\t{'mCO2/SO2'}\t{'mH2S/SO2'}\t{'fH2'}\t{'fH2O'}\t"
        f"{'fCO2'}\t{'fCO'}\t{'fCH4'}\t{'fSO2'}\t{'fH2S'}\t{'fS2'}\t{'fN2'}\t{'tot_H'}"
        f"\t{'tot_C'}\t{'tot_O_gas'}\t{'tot_O'}\t{'tot_S'}\t{'tot_N'}\n"
    )

    i = 0

    if sys.run.FIND_SATURATION is True or sys.run.ATOMIC_MASS_SET is True:
        output.write(
            sat.satp_writeout(
                sys,
                melt,
                gas,
                sys.sat_conditions[0],
                sys.sat_conditions[1],
                sys.sat_conditions[2],
                sys.sat_conditions[3],
                graph_sat=sys.sat_conditions[4],
            )
        )

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        if len(P) > 1:
            while i < len(P):
                output.write(
                    f"{P[i]} \t "
                    f"{cnvt.generate_fo2_buffer(sys, np.exp(gas.fo2[i]), P[i]):+.3} \t "
                    f"{np.exp(gas.fo2[i]):.3e} \t {melt.F[i+1]:.3} \t {sys.rho[i]:.3} "
                    f"\t {melt.rho_store[i]:.3} \t {sys.GvF[i] * 100:} \t "
                    f"{sys.WgT[i+1] * 100:.5g} \t {gas.M[i]} \t {gas.mH2O[i+1]:.8g} \t "
                    f"{gas.mH2[i+1]:.8g} \t {gas.mO2[i+1]:.6g} \t {gas.mCO2[i+1]:.6g} "
                    f"\t {gas.mCO[i+1]:.6g} \t {gas.mCH4[i+1]:.6g} \t "
                    f"{gas.mSO2[i+1]:.6g} \t {gas.mH2S[i+1]:.6g} \t {gas.mS2[i+1]:.6g} "
                    f"\t {gas.mN2[i+1]:.6g} \t {gas.Wt['H2O'][i]:.6g} \t "
                    f"{gas.Wt['H2'][i]:.6g} \t {gas.Wt['O2'][i]:.6g} \t "
                    f"{gas.Wt['CO2'][i]:.6g} \t {gas.Wt['CO'][i]:.6g} \t "
                    f"{gas.Wt['CH4'][i]:.6g} \t {gas.Wt['SO2'][i]:.6g} \t "
                    f"{gas.Wt['H2S'][i]:.6g} \t {gas.Wt['S2'][i]:.6g} \t "
                    f"{gas.Wt['N2'][i]:.6g} \t {melt.h2o[i]*100:.6g} \t "
                    f"{melt.h2[i]*100:.6g} \t {melt.co2[i]*100:.6g} \t "
                    f"{melt.co[i]*100:.6g} \t {melt.ch4[i]*100:.6g} \t "
                    f"{melt.graphite[i]*100:.6g} \t {melt.sulfide[i]*100:.6g} \t "
                    f"{melt.sulfate[i]*100:.6g} \t {melt.s[i]*100:.6g} \t "
                    f"{melt.n[i]*100:.6g} \t {gas.mCO2[i+1]/gas.mCO[i+1]:.6g} \t "
                    f"{gas.mCO2[i+1]/gas.mH2O[i+1]:.6g} \t "
                    f"{gas.mCO2[i+1]/gas.mSO2[i+1]:.6g} \t "
                    f"{gas.mH2S[i+1]/gas.mSO2[i+1]:.6g} \t {gas.f['H2'][i]:.8g} \t "
                    f"{gas.f['H2O'][i]:.6g} \t {gas.f['CO2'][i]:.6g} \t "
                    f"{gas.f['CO'][i]:.6g} \t {gas.f['CH4'][i]:.6g} \t "
                    f"{gas.f['SO2'][i]:.8g} \t {gas.f['H2S'][i]:.8g} \t "
                    f"{gas.f['S2'][i]:.8g} \t {gas.f['N2'][i]:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 'h', i)*1000000:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 'c', i)*1000000:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 'o', i)*1000000:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 'o_tot', i)*1000000:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 's', i)*1000000:.8g} \t "
                    f"{cnvt.atomicM_calc(sys, melt, gas, 'n', i)*1000000:.8g} \n"
                )
                i += 1
        else:
            output.write(
                f"{P[-1]} \t "
                f"{cnvt.generate_fo2_buffer(sys, np.exp(gas.fo2[-1]), P[-1]):+.3} \t "
                f"{np.exp(gas.fo2[-1]):.3e} \t {melt.F[-1]:.3} \t {sys.rho[-1]:.3} \t "
                f"{melt.rho_store[-1]:.3} \t {sys.GvF[-1] * 100} \t "
                f"{sys.WgT[-1] * 100:.5g} \t {gas.M[-1]} \t {gas.mH2O[-1]:.8g} \t "
                f"{gas.mH2[-1]:.8g} \t {gas.mO2[-1]:.6g} \t {gas.mCO2[-1]:.6g} \t "
                f"{gas.mCO[-1]:.6g} \t {gas.mCH4[-1]:.6g} \t {gas.mSO2[-1]:.6g} \t "
                f"{gas.mH2S[-1]:.6g} \t {gas.mS2[-1]:.6g} \t {gas.mN2[-1]:.6g} \t "
                f"{gas.Wt['H2O'][-1]:.6g} \t {gas.Wt['H2'][-1]:.6g} \t "
                f"{gas.Wt['O2'][-1]:.6g} \t {gas.Wt['CO2'][-1]:.6g} \t "
                f"{gas.Wt['CO'][-1]:.6g} \t {gas.Wt['CH4'][-1]:.6g} \t "
                f"{gas.Wt['SO2'][-1]:.6g} \t {gas.Wt['H2S'][-1]:.6g} \t "
                f"{gas.Wt['S2'][-1]:.6g} \t {gas.Wt['N2'][-1]:.6g} \t "
                f"{melt.h2o[-1]*100:.6g} \t {melt.h2[-1]*100:.6g} \t "
                f"{melt.co2[-1]*100:.6g} \t {melt.co[-1]*100:.6g} \t "
                f"{melt.ch4[-1]*100:.6g} \t {melt.graphite[-1]*100:.6g} \t "
                f"{melt.sulfide[-1]*100:.6g} \t {melt.sulfate[-1]*100:.6g} \t "
                f"{melt.s[-1]*100:.6g} \t {melt.n[-1]*100:.6g} \t "
                f"{gas.mCO2[-1]/gas.mCO[-1]:.6g} \t {gas.mCO2[-1]/gas.mH2O[-1]:.6g} \t "
                f"{gas.mCO2[-1]/gas.mSO2[-1]:.6g} \t {gas.mH2S[-1]/gas.mSO2[-1]:.6g} \t"
                f" {gas.f['H2'][-1]:.8g} \t {gas.f['H2O'][-1]:.6g} \t "
                f"{gas.f['CO2'][-1]:.8g} \t {gas.f['CO'][-1]:.8g} \t "
                f"{gas.f['CH4'][-1]:.6g} \t {gas.f['SO2'][-1]:.8g} \t "
                f"{gas.f['H2S'][-1]:.8g} \t {gas.f['S2'][-1]:.8g} \t "
                f"{gas.f['N2'][-1]:.8g} \t "
                f"{cnvt.atomicM_calc(sys, melt, gas, 'h', -1, WgT=sys.WgT[-1])*1000000:.8g} \t "  # noqa:E501
                f"{cnvt.atomicM_calc(sys, melt, gas, 'c', -1, WgT=sys.WgT[-1])*1000000:.8g} \t "  # noqa:E501
                f"{cnvt.atomicM_calc(sys, melt, gas, 'o', -1, WgT=sys.WgT[-1])*1000000:.8g} \t "  # noqa:E501
                f"{cnvt.atomicM_calc(sys, melt, gas, 'o_tot', -1, WgT=sys.WgT[-1])*1000000:.8g} \t "  # noqa:E501
                f"{cnvt.atomicM_calc(sys, melt, gas, 's', -1, WgT=sys.WgT[-1])*1000000:.8g} \t "  # noqa:E501
                f"{cnvt.atomicM_calc(sys, melt, gas, 'n', -1, WgT=sys.WgT[-1])*1000000:.8g} \n"  # noqa:E501
            )

    output.close()


def writeout_figs(sys, melt, gas, out, P):
    """
    Selects the figures to produce at the end of a successful run.

    If an output file is not used to select the required figures, all
    options are plotted.

    Parameters
    ----------
    sys : ThermoSystem class
        Active instance of the ThermoSystem class
    melt : Melt class
        Active instance of the Melt class
    gas : Gas class
        Active instance of the Gas class
    out : Output class or None
        Active instance of the Output class, or None
    P : list of float
        List of the pressures each step was calculated at (bar)
    """

    filelist = glob.glob("Output/*.png")
    # Removes previous files so if output specification is changed
    # there is no confusion as to up to date files.

    for file in filelist:
        os.remove(file)
    if (
        out is not None
    ):  # If an output file listing requested figures has been included:
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

        if out.Plt_gas_fraction is not None:
            pass
            # plot

    else:  # else plot every option
        plot_meltspecies(melt, P)
        plot_gasspecies_wt(gas, P)
        plot_gasspecies_mol(gas, P)
        plot_gasfraction(sys, P)
        plot_fo2FMQ(melt, gas, P)


# fO2 and dFMQ
def plot_fo2FMQ(melt, gas, P):
    """Plots the fO2 relative to FMQ against pressure"""
    plt.plot(P, melt.fmq)
    plt.xscale("log")
    plt.xlabel("Pressure (bars)")
    plt.ylabel(r"$\Delta$ FMQ")
    plt.savefig("Output/FMQ.png")
    plt.close()

    fo2 = []
    for x in gas.fo2:
        fo2.append(np.log10(np.exp(x)))

    plt.plot(P, fo2)
    # plt.xscale('log')
    plt.xlabel("Pressure (bars)")
    plt.ylabel("log(10) fO2")
    plt.savefig("Output/fO2.png")
    plt.close()


# Gas speciation mol
def plot_gasspecies_mol(gas, P):
    """Plots the gas speciation (mole fraction) vs pressure"""

    plt.plot(gas.mH2O[1:], P, c="darkblue")
    plt.plot(gas.mH2[1:], P, c="steelblue")
    plt.plot(gas.mO2[1:], P, c="firebrick")
    plt.annotate("H2O", xy=[gas.mH2O[-1], P[-1]], color="darkblue")
    plt.annotate("H2", xy=[gas.mH2[-1], P[-1]], color="steelblue")
    plt.annotate("O2", xy=[gas.mO2[-1], P[-1]], color="firebrick")
    if (
        gas.sys.run.GAS_SYS == "COH"
        or gas.sys.run.GAS_SYS == "COHS"
        or gas.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(gas.mCO2[1:], P, c="saddlebrown")
        plt.plot(gas.mCO[1:], P, c="darkgreen")
        plt.plot(gas.mCH4[1:], P, c="forestgreen")
        plt.annotate("CO2", xy=[gas.mCO2[-1], P[-1]], color="saddlebrown")
        plt.annotate("CO", xy=[gas.mCO[-1], P[-1]], color="darkgreen")
        plt.annotate("CH4", xy=[gas.mCH4[-1], P[-1]], color="forestgreen")
    if (
        gas.sys.run.GAS_SYS == "SOH"
        or gas.sys.run.GAS_SYS == "COHS"
        or gas.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(gas.mSO2[1:], P, c="maroon")
        plt.plot(gas.mS2[1:], P, c="goldenrod")
        plt.plot(gas.mH2S[1:], P, c="pink")
        plt.annotate("SO2", xy=[gas.mSO2[-1], P[-1]], color="maroon")
        plt.annotate("S2", xy=[gas.mS2[-1], P[-1]], color="goldenrod")
        plt.annotate("H2S", xy=[gas.mH2S[-1], P[-1]], color="pink")
    if gas.sys.run.GAS_SYS == "COHSN":
        plt.plot(gas.mN2[1:], P, c="grey")
        plt.annotate("N2", xy=[gas.mN2[-1], P[-1]], color="grey")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel(f"Speciation in a {gas.sys.run.GAS_SYS} gas (mol frac)")
    plt.ylabel("Pressure (bar)")
    plt.savefig("Output/speciation(mol).png")
    plt.close()


# Gas speciation wt%


def plot_gasspecies_wt(gas, P):
    """Plots the gas speciation (weight fraction) vs pressure"""

    plt.plot(gas.Wt["H2O"], P, c="darkblue")
    plt.plot(gas.Wt["H2"], P, c="steelblue")
    plt.plot(gas.Wt["O2"], P, c="firebrick")
    plt.annotate("H2O", xy=[gas.Wt["H2O"][-1], P[-1]], color="darkblue")
    plt.annotate("H2", xy=[gas.Wt["H2"][-1], P[-1]], color="steelblue")
    plt.annotate("O2", xy=[gas.Wt["O2"][-1], P[-1]], color="firebrick")
    if (
        gas.sys.run.GAS_SYS == "COH"
        or gas.sys.run.GAS_SYS == "COHS"
        or gas.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(gas.Wt["CO2"], P, c="saddlebrown")
        plt.plot(gas.Wt["CO"], P, c="darkgreen")
        # plt.plot(gas.Wt['CH4'], P, c='mediumseagreen')
        plt.annotate("CO2", xy=[gas.Wt["CO2"][-1], P[-1]], color="saddlebrown")
        plt.annotate("CO", xy=[gas.Wt["CO"][-1], P[-1]], color="darkgreen")
        # plt.annotate('CH4', xy=[gas.mCH4[-1], 0.9], color='mediumseagreen')
    if (
        gas.sys.run.GAS_SYS == "SOH"
        or gas.sys.run.GAS_SYS == "COHS"
        or gas.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(gas.Wt["SO2"], P, c="maroon")
        plt.plot(gas.Wt["S2"], P, c="goldenrod")
        plt.plot(gas.Wt["H2S"], P, c="pink")
        plt.annotate("SO2", xy=[gas.Wt["SO2"][-1], P[-1]], color="maroon")
        plt.annotate("S2", xy=[gas.Wt["S2"][-1], P[-1]], color="goldenrod")
        plt.annotate("H2S", xy=[gas.Wt["H2S"][-1], P[-1]], color="pink")
    if gas.sys.run.GAS_SYS == "COHSN":
        plt.plot(gas.Wt["N2"], P, c="grey")
        plt.annotate("N2", xy=[gas.Wt["N2"][-1], P[-1]], color="grey")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel(f"Gas phase speciation of a {gas.sys.run.GAS_SYS} system (wt %)")
    plt.ylabel("Pressure (bar)")
    plt.savefig("Output/speciation(wt).png")
    plt.close()


# Melt volatile wt%


def plot_meltspecies(melt, P):
    """Plots the volatile content of the melt vs pressure"""

    plt.plot(melt.h2o, P, c="darkblue")
    plt.plot(melt.h2, P, c="steelblue")
    plt.annotate("H2O", xy=[melt.h2o[-1], P[-1]], color="darkblue")
    plt.annotate("H2", xy=[melt.h2[-1], P[-1]], color="steelblue")
    if (
        melt.sys.run.GAS_SYS == "COH"
        or melt.sys.run.GAS_SYS == "COHS"
        or melt.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(melt.co2, P, c="saddlebrown")
        plt.annotate("CO2", xy=[melt.co2[-1], P[-1]], color="saddlebrown")
    if (
        melt.sys.run.GAS_SYS == "SOH"
        or melt.sys.run.GAS_SYS == "COHS"
        or melt.sys.run.GAS_SYS == "COHSN"
    ):
        plt.plot(melt.s, P, c="maroon")
        plt.annotate("S", xy=[melt.s[-1], P[-1]], color="maroon")
    if melt.sys.run.GAS_SYS == "COHSN":
        plt.plot(melt.n, P, c="grey")
        plt.annotate("N", xy=[melt.n[-1], P[-1]], color="grey")
    plt.xscale("log")
    plt.gca().invert_yaxis()
    plt.xlabel("Melt volatile content (wt%)")
    plt.ylabel("Pressure (bar)")
    plt.savefig("Output/meltspecies.png")
    plt.close()


# Exsolved gas mass and volume


def plot_gasfraction(sys, P):
    """Plots the gas volume fraction vs pressure"""

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.plot(cnvt.frac2perc(sys.WgT[1:]), P)
    ax1.invert_yaxis()
    ax1.set_xlabel("Exsolved gas wt%")
    ax1.set_ylabel("Pressure (bar)")
    ax2.plot(cnvt.frac2perc(sys.GvF), P)
    ax2.set_xlabel("Exsolved gas volume %")
    fig.savefig("Output/exsolved_gas.png")
