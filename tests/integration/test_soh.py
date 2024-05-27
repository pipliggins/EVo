import evo
import pandas as pd


def test_soh_default():
    df = evo.main(
        "evo/input/chem.yaml", "tests/integration/input_files/env_soh.yaml", None
    )

    assert isinstance(df, pd.DataFrame)
