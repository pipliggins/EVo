import matplotlib.pyplot as plt
import numpy as np


# fO2 and dFMQ
def plot_fo2FMQ(df):
    """Plots the fO2 relative to FMQ against pressure"""
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.plot(df["FMQ"], df["P"])
    ax1.set_yscale("log")
    ax1.set_ylabel("Pressure (bars)")
    ax1.invert_yaxis()
    ax1.set_xlabel(r"$\Delta$ FMQ")

    fo2 = []
    for x in df["fo2"]:
        fo2.append(np.log10(np.exp(x)))

    ax2.plot(fo2, df["P"])
    ax2.set_xlabel("log(10) fO2")

    return fig


# Gas speciation mol
def plot_gasspecies_mol(df):
    """Plots the gas speciation (mole fraction) vs pressure"""

    fig, ax = plt.subplots()
    ax.plot(df["mH2O"], df["P"], c="darkblue", label="H2O")
    ax.plot(df["mH2"], df["P"], c="steelblue", label="H2")
    ax.plot(df["mO2"], df["P"], c="firebrick", label="O2")
    ax.annotate("H2O", xy=[df["mH2O"].iloc[-1], df["P"].iloc[-1]], color="darkblue")
    ax.annotate("H2", xy=[df["mH2"].iloc[-1], df["P"].iloc[-1]], color="steelblue")
    ax.annotate("O2", xy=[df["mO2"].iloc[-1], df["P"].iloc[-1]], color="firebrick")

    ax.plot(df["mCO2"], df["P"], c="saddlebrown", label="CO2")
    ax.plot(df["mCO2"], df["P"], c="darkgreen", label="CO")
    ax.plot(df["mCH4"], df["P"], c="forestgreen", label="CH4")
    ax.annotate("CO2", xy=[df["mCO2"].iloc[-1], df["P"].iloc[-1]], color="saddlebrown")
    ax.annotate("CO", xy=[df["mCO2"].iloc[-1], df["P"].iloc[-1]], color="darkgreen")
    ax.annotate("CH4", xy=[df["mCH4"].iloc[-1], df["P"].iloc[-1]], color="forestgreen")

    ax.plot(df["mSO2"], df["P"], c="maroon", label="SO2")
    ax.plot(df["mS2"], df["P"], c="goldenrod", label="S2")
    ax.plot(df["mH2S"], df["P"], c="pink", label="H2S")
    ax.annotate("SO2", xy=[df["mSO2"].iloc[-1], df["P"].iloc[-1]], color="maroon")
    ax.annotate("S2", xy=[df["mS2"].iloc[-1], df["P"].iloc[-1]], color="goldenrod")
    ax.annotate("H2S", xy=[df["mH2S"].iloc[-1], df["P"].iloc[-1]], color="pink")

    ax.plot(df["mN2"], df["P"], c="grey", label="N2")
    ax.annotate("N2", xy=[df["mN2"].iloc[-1], df["P"].iloc[-1]], color="grey")

    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_xlabel("Gas speciation (mol frac)")
    ax.set_ylabel("Pressure (bar)")
    return fig


# Gas speciation wt%


def plot_gasspecies_wt(df):
    """Plots the gas speciation (weight fraction) vs pressure"""

    fig, ax = plt.subplots()

    for g, color in zip(
        ["wH2O", "wH2", "wO2", "wCO2", "wCO", "wCH4", "wSO2", "wS2", "wH2S", "wN2"],
        [
            "darkblue",
            "steelblue",
            "firebrick",
            "saddlebrown",
            "darkgreen",
            "mediumseagreen",
            "maroon",
            "goldenrod",
            "pink",
            "grey",
        ],
    ):
        ax.plot(df[g], df["P"], c=color, label=g[1:])
        ax.annotate(g[1:], xy=[df[g].iloc[-1], df["P"].iloc[-1]], color=color)

    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_xlabel("Gas phase (wt %)")
    ax.set_ylabel("Pressure (bar)")
    return fig


# Melt volatile wt%


def plot_meltspecies(df):
    """Plots the volatile content of the melt vs pressure"""

    fig, ax = plt.subplots()
    for s, color in zip(
        [
            "H2O_melt",
            "H2_melt",
            "CO2_melt",
            "CO_melt",
            "CH4_melt",
            "S2-_melt",
            "S6+_melt",
            "N_melt",
        ],
        [
            "darkblue",
            "steelblue",
            "saddlebrown",
            "darkgreen",
            "mediumseagreen",
            "maroon",
            "pink",
            "grey",
        ],
    ):
        ax.plot(df[s], df["P"], c=color, label=s.split("_")[0])
        ax.annotate(s.split("_")[0], xy=[df[s].iloc[-1], df["P"].iloc[-1]], color=color)

    ax.set_xscale("log")
    ax.invert_yaxis()
    ax.set_xlabel("Melt volatile content (wt%)")
    ax.set_ylabel("Pressure (bar)")
    return fig


# Exsolved gas mass and volume


def plot_gasfraction(df):
    """Plots the gas volume fraction vs pressure"""

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    ax1.plot(df["Gas_wt"], df["P"])
    ax1.invert_yaxis()
    ax1.set_xlabel("Exsolved gas wt%")
    ax1.set_ylabel("Pressure (bar)")
    ax2.plot(df["Exsol_vol%"], df["P"])
    ax2.set_xlabel("Exsolved gas volume %")
    return fig
