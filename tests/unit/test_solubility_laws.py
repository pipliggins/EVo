import evo.dgs_classes
import evo.solubility_laws as sl
import evo
import pytest
import numpy as np


@pytest.mark.parametrize(
    "fCO, P, expected",
    [
        (10**2.5, 0, 30e-6),
        (10**2.5, 10000, 8e-6),
        (10**2.5, 20000, 2e-6),
        (10**3.5, 10000, 60e-6),
        (10**3.5, 20000, 18e-6),
        (10**3.5, 30000, 5e-6),
    ],
)
def test_armstrong2015(fCO, P, expected):
    """
    Expected values inferred from Fig.7 of Armstrong et al. (2015)
    fCO as an absolute value (i.e. not log(fCO) as in plot)
    P in bar, not GPa
    Output in weight fraction, not ppm.
    """
    assert sl.armstrong2015(fCO, P) == pytest.approx(expected, 0.16)


@pytest.mark.parametrize(
    "composition, T, expected",
    [
        (
            {
                "na2o": 2.63,
                "mgo": 5.88,
                "al2o3": 13.36,
                "sio2": 50.93,
                "k2o": 0.17,
                "cao": 10.81,
                "tio2": 2.46,
                "mno": 0.25,
                "feo": 14.21,
                "fe2o3": 0,
            },
            1432,
            -2.25,
        ),
        (
            {
                "na2o": 2.27,
                "mgo": 9.18,
                "al2o3": 16.3,
                "sio2": 49.83,
                "k2o": 0.05,
                "cao": 12.11,
                "tio2": 0.87,
                "mno": 0.17,
                "feo": 9.37,
                "fe2o3": 0,
            },
            1494,
            -2.44,
        ),
    ],
)
def test_oneill2020(composition, T, expected):
    """
    Tested against calculation spreadsheet provided by O'Neill & Mavrogenes (2022),
    sheet "Table S6 S redox calculator", expected value = CS_2-
    """
    run = evo.dgs_classes.RunDef()
    sys = evo.dgs_classes.ThermoSystem(run)
    sys.T = T
    melt = evo.dgs_classes.Melt(run, sys)

    melt.cm_dry = evo.conversions.wt2mol(composition)

    assert sl.oneill2020(T, melt) == pytest.approx(np.exp(expected) / 1e6, 0.01)
