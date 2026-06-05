import yaml
import evo
import pandas as pd

_CHEM = "tests/integration/input_files/chem.yaml"
_ENV = "tests/integration/input_files/env_cohs.yaml"


def _load_env(path, **overrides):
    with open(path) as f:
        d = yaml.full_load(f)
    d.update(overrides)
    return pd.Series(d)


def test_cohs_saturation(tmp_path):
    df = evo.run_evo(_CHEM, _ENV, None, folder=tmp_path)
    assert isinstance(df, pd.DataFrame)


def test_cohs_atomic_mass(tmp_path):
    env = _load_env(
        _ENV,
        ATOMIC_MASS_SET=True,
        FIND_SATURATION=False,
        WTH2O_SET=False,
        WTCO2_SET=False,
        SULFUR_SET=False,
    )
    df = evo.run_evo(_CHEM, env, None, folder=tmp_path)
    assert isinstance(df, pd.DataFrame)
