import evo
import pandas as pd


def test_cohsn_default():
    df = evo.main(
        "evo/input/chem.yaml", "tests/integration/input_files/env_cohsn.yaml", None
    )

    assert isinstance(df, pd.DataFrame)
