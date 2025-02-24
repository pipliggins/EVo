import evo
import pandas as pd


def test_cohs_default(tmp_path):
    df = evo.main(
        "evo/input/chem.yaml",
        "tests/integration/input_files/env_cohs.yaml",
        None,
        folder=tmp_path,
    )

    assert isinstance(df, pd.DataFrame)
