"""
EVo

A Python model for volcanic degassing, using the equilibrium constants and mass balance
method.
"""

__version__ = "1.0.1"
__author__ = "Philippa Liggins"
__all__ = ["main", "multirun"]

# ----------------- IMPORTS ----------------- #
from evo.dgs import main
from evo.multirun import multirun
