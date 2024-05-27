import evo
import pandas as pd


def test_oh_default():
    df = evo.main(
        "evo/input/chem.yaml", "tests/integration/input_files/env_oh.yaml", None
    )

    assert isinstance(df, pd.DataFrame)
