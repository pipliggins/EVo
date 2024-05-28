"""
Tests the warnings provided when a system is being set up with an invalid
set of volatile inputs.
"""

import pytest
import evo
from evo.messages import vol_setup_standard, InvalidSetupError


@pytest.mark.parametrize(
    "options, error_msg",
    # fh2, fh2o, fco2, wth2o, wtco2, wts, wtn, graph_sat
    [
        (
            [True, False, False, False, False, False, False, False],
            "Set only one of either fH2, fO2",
        ),
        (
            [False, True, False, False, False, False, False, False],
            "The OH system cannot be set up with",
        ),
        (
            [False, False, True, False, False, False, False, False],
            "The OH system cannot be set up with",
        ),
        (
            [False, False, False, True, False, False, False, False],
            "The OH system cannot be set up with",
        ),
        (
            [False, False, False, False, True, False, False, False],
            "The OH system cannot be set up with",
        ),
        (
            [False, False, False, False, False, True, False, False],
            "The OH system cannot be set up with",
        ),
        (
            [False, False, False, False, False, False, True, False],
            "The OH system cannot be set up with",
        ),
        (
            [False, False, False, False, False, False, False, True],
            "The OH system cannot be set up with",
        ),
        (
            [False, True, False, True, False, False, False, False],
            "Only use one option to set H2O or CO2",
        ),
        (
            [False, False, True, False, True, False, False, False],
            "Only use one option to set H2O or CO2",
        ),
        (
            [False, False, True, False, False, False, False, True],
            "Only use one setting for CO2.",
        ),
        (
            [False, False, False, False, True, False, False, True],
            "Only use one setting for CO2.",
        ),
    ],
)
def test_oh_standard_setup_errors_fO2isSet(options, error_msg):
    fh2, fh2o, fco2, wth2o, wtco2, wts, wtn, graph_sat = options

    run = evo.dgs_classes.RunDef()
    run.GAS_SYS = "OH"
    run.FO2_buffer_SET = True  # sets fO2 to FMQ
    run.FH2_SET = fh2
    run.FH2O_SET = fh2o
    run.FCO2_SET = fco2
    run.WTH2O_SET = wth2o
    run.WTCO2_SET = wtco2
    run.SULFUR_SET = wts
    run.NITROGEN_SET = wtn
    run.GRAPHITE_SATURATED = graph_sat
    sys = evo.dgs_classes.ThermoSystem(run)
    sys.FO2 = 1  # value unimportant

    with pytest.raises(InvalidSetupError, match=error_msg):
        vol_setup_standard(run, sys)


@pytest.mark.parametrize(
    "options, error_msg",
    # fh2, fh2o, fco2, wth2o, wtco2, wts, wtn, graph_sat
    [
        (
            [True, True, False, True, False, False, False, False],
            "Only use one option to set H2O or CO2",
        ),
        (
            [True, True, False, False, False, False, False, False],
            "Set exactly 2 initial conditions from either fH2, H2O ",
        ),
        (
            [False, False, True, False, False, False, False, False],
            "C, S and N species cannot be used as an initial condition for",
        ),
        (
            [True, False, False, True, False, False, False, False],
            "Set exactly 2 initial conditions from either fH2, H2O ",
        ),
        (
            [False, False, False, False, True, False, False, False],
            "C, S and N species cannot be used as an initial condition for",
        ),
        (
            [False, False, False, False, False, True, False, False],
            "C, S and N species cannot be used as an initial condition for",
        ),
        (
            [False, False, False, False, False, False, True, False],
            "C, S and N species cannot be used as an initial condition for",
        ),
        (
            [False, False, False, False, False, False, False, True],
            "C, S and N species cannot be used as an initial condition for",
        ),
    ],
)
def test_coh_soh_standard_setup_errors_fO2isSet(options, error_msg):
    fh2, fh2o, fco2, wth2o, wtco2, wts, wtn, graph_sat = options

    run = evo.dgs_classes.RunDef()
    run.GAS_SYS = "COH"
    run.FO2_buffer_SET = True  # sets fO2 to FMQ
    run.FH2_SET = fh2
    run.FH2O_SET = fh2o
    run.FCO2_SET = fco2
    run.WTH2O_SET = wth2o
    run.WTCO2_SET = wtco2
    run.SULFUR_SET = wts
    run.NITROGEN_SET = wtn
    run.GRAPHITE_SATURATED = graph_sat
    sys = evo.dgs_classes.ThermoSystem(run)
    sys.FO2 = 1  # value unimportant

    with pytest.raises(InvalidSetupError, match=error_msg):
        vol_setup_standard(run, sys)

    sys.run.GAS_SYS = "SOH"
    with pytest.raises(InvalidSetupError, match=error_msg):
        vol_setup_standard(run, sys)


@pytest.mark.parametrize(
    "options, error_msg",
    # fh2, fh2o, fco2, wth2o, wtco2, wts, wtn, graph_sat
    [
        (
            [True, True, True, False, False, False, False, False],
            "Set exactly 3 initial conditions from either fH2, H2O ",
        ),
        (
            [False, False, False, False, True, True, False, True],
            "Only use one setting for CO2.",
        ),
        # (
        #     [False, False, False, False, True, True, True, False],
        #     "Set exactly 3 initial conditions from either fH2, H2O",
        # ), # output printed to terminal
    ],
)
def test_cohs_standard_setup_errors_fO2isSet(options, error_msg):
    fh2, fh2o, fco2, wth2o, wtco2, wts, wtn, graph_sat = options

    run = evo.dgs_classes.RunDef()
    run.GAS_SYS = "COHS"
    run.FO2_buffer_SET = True  # sets fO2 to FMQ
    run.FH2_SET = fh2
    run.FH2O_SET = fh2o
    run.FCO2_SET = fco2
    run.WTH2O_SET = wth2o
    run.WTCO2_SET = wtco2
    run.SULFUR_SET = wts
    run.NITROGEN_SET = wtn
    run.GRAPHITE_SATURATED = graph_sat
    sys = evo.dgs_classes.ThermoSystem(run)
    sys.FO2 = 1  # value unimportant

    with pytest.raises(InvalidSetupError, match=error_msg):
        vol_setup_standard(run, sys)
