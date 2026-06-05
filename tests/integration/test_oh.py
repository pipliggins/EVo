import yaml
import evo
import pandas as pd

_CHEM = "tests/integration/input_files/chem.yaml"
_ENV = "tests/integration/input_files/env_oh.yaml"


def _load_env(path, **overrides):
    with open(path) as f:
        d = yaml.full_load(f)
    d.update(overrides)
    return pd.Series(d)


def test_oh_default(tmp_path):
    df = evo.run_evo(_CHEM, _ENV, None, folder=tmp_path)
    assert isinstance(df, pd.DataFrame)


def test_oh_saturation(tmp_path):
    env = _load_env(
        _ENV,
        FIND_SATURATION=True,
        FO2_buffer_SET=True,
        FO2_buffer="FMQ",
        FO2_buffer_START=0,
        WTH2O_SET=True,
        WTH2O_START=0.03,
        FH2_SET=False,
    )
    df = evo.run_evo(_CHEM, env, None, folder=tmp_path)
    assert isinstance(df, pd.DataFrame)


def test_oh_atomic_mass(tmp_path):
    env = _load_env(
        _ENV,
        ATOMIC_MASS_SET=True,
        FO2_buffer_SET=True,
        FO2_buffer="FMQ",
        FO2_buffer_START=0,
        FH2_SET=False,
    )
    df = evo.run_evo(_CHEM, env, None, folder=tmp_path)
    assert isinstance(df, pd.DataFrame)
