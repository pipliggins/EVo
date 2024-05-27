import evo
import pandas as pd


def test_coh_default():
    df = evo.main(
        "evo/input/chem.yaml", "tests/integration/input_files/env_coh.yaml", None
    )

    assert isinstance(df, pd.DataFrame)
