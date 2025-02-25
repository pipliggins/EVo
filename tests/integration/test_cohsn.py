import evo
import pandas as pd


def test_cohsn_default(tmp_path):
    df = evo.run_evo(
        "tests/integration/input_files/chem.yaml",
        "tests/integration/input_files/env_cohsn.yaml",
        None,
        folder=tmp_path,
    )

    assert isinstance(df, pd.DataFrame)
