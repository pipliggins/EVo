"""
EVo
A Python model for volcanic degassing, using the equilibrium constants and mass balance
method.
"""

from .dgs import run_evo as run_evo

__all__ = ["run_evo"]
