import evo
import pandas as pd


def test_cohs_default():
    df = evo.main(
        "evo/input/chem.yaml", "tests/integration/input_files/env_cohs.yaml", None
    )

    assert isinstance(df, pd.DataFrame)
