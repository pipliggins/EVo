"""
Tests the functions in ferric.py
Kress & Carmichael (1991) is tested using values calculated from
acovino, Kayla. (2021). Ferric/Ferrous, Fe3+/FeT, fO2 Converter (Kress and Carmichael,
1991) (3.2). Zenodo. https://doi.org/10.5281/zenodo.5907844
"""

import evo
import pytest
import numpy as np


def test_kc91_fo2_where_correctInput_gives_appropriateValue():
    test_composition = {
        "sio2": 0.3931,  # mole fractions
        "tio2": 0.0009,
        "al2o3": 0.0069,
        "feo": 0.0593,
        "mno": 0.001,
        "mgo": 0.5078,
        "cao": 0.0299,
        "na2o": 0.0011,
        "k2o": 0.000034,
        "p2o5": 0.0001,
    }
    assert evo.ferric.kc91_fo2(
        test_composition, 1473.15, 1, -21.415901
    ) == pytest.approx(0.0467, 0.001)


def test_kc91_fe3_where_correctInput_gives_appropriateValue():
    test_composition = {
        "sio2": 0.825741,  # mole fractions
        "tio2": 0.001943,
        "al2o3": 0.060749,
        "fe2o3": 0.016706,
        "feo": 0.014331,
        "mno": 0.002377,
        "mgo": 0.0,
        "cao": 0.002406,
        "na2o": 0.045508,
        "k2o": 0.030239,
        "p2o5": 0.0,
    }
    assert evo.ferric.kc91_fe3(test_composition, 1473.15, 100000) == pytest.approx(
        0.0015106, 0.001
    )


def test_r2013_fo2_where_correctInput_gives_appropriateValue():
    test_composition = {
        "sio2": 0.48121,  # mole fractions
        "tio2": 0.003642,
        "al2o3": 0.060635,
        "feo": 0.181667,
        "mno": 0.002792,
        "mgo": 0.159075,
        "cao": 0.082879,
        "na2o": 0.025262,
        "k2o": 0.00092,
        "p2o5": 0.001918,
    }
    assert evo.ferric.r2013_fo2(
        test_composition, 1679.15, 1.2e9, -17.09620888
    ) == pytest.approx(0.117923975, 0.0001)


def test_r2013_fe3_where_correctInput_gives_appropriateValue():
    test_composition = {
        "sio2": 0.5771,  # mole fractions
        "tio2": 0.0102,
        "al2o3": 0.0262,
        "fe2o3": 0.04,
        "feo": 0.1285,
        "mno": 0.0023,
        "mgo": 0.0785,
        "cao": 0.1155,
        "na2o": 0.0183,
        "k2o": 0.00173,
        "p2o5": 0.0016,
    }
    assert np.log(
        evo.ferric.r2013_fe3(test_composition, 1473.15, 1e9)
    ) == pytest.approx(-11.3887, 0.001)
