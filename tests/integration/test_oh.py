import evo
import pandas as pd


def test_oh_default(tmp_path):
    df = evo.run_evo(
        "tests/integration/input_files/chem.yaml",
        "tests/integration/input_files/env_oh.yaml",
        None,
        folder=tmp_path,
    )

    assert isinstance(df, pd.DataFrame)
