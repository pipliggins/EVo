"""
EVo

A Python model for volcanic degassing, using the equilibrium constants and mass balance
method.
"""

from evo.dgs import run_evo
from evo.multirun import multirun

__all__ = ["multirun", "run_evo"]
