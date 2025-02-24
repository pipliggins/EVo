"""
EVo

A Python model for volcanic degassing, using the equilibrium constants and mass balance
method.
"""

# ----------------- IMPORTS ----------------- #
from evo.dgs import run_evo
from evo.multirun import multirun

__version__ = "1.0.2"
__author__ = "Philippa Liggins"
__all__ = ["run_evo", "multirun"]
