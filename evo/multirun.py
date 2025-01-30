import os
from itertools import product
from shutil import copyfile

import numpy as np
import ruamel.yaml

from evo.dgs import main

ryaml = ruamel.yaml.YAML()


def amend_env(file, **kwargs):
    """
    Duplicates the environment input file to use in `multirun()`

    Duplicates the environment file, then edits the file with the relevant
    values for performing a run within the gridsearch. Saves the duplicated
    file as 'multirun.yaml'.

    Parameters
    ----------
    file : str
        Path to the standard environment file
    **kwargs :
        str:any pairs corresponding to items in the env.yaml file
    """

    with open(file) as f:
        env_doc = ryaml.load(f)

    for param, val in kwargs.items():
        if isinstance(val, np.float64):
            env_doc[param] = float(val)
        else:
            env_doc[param] = val

    with open("multirun.yaml", "w") as f:
        ryaml.dump(env_doc, f)


def multirun(**kwargs):
    """
    Performs a gridsearch, running EVo over all vars given as parameters

    Gridsearch over all the kwargs given (as lists or arrays), having
    setup variables you don't want to search over within the env.yaml
    already.
    Each run result will be saved as a separate output file in Outputs, named using the
    variable values used for that run.

    Parameters
    ----------
    **kwargs :
        str:list of floats pairs corresponding to items in the env.yaml file
    """

    # delete any pre-existing multirun setup files
    if os.path.exists("multirun.yaml"):
        os.remove("multirun.yaml")

    # creates an iterable object where each 'column' is in the order given here.
    options = product(*kwargs.values())
    keys = kwargs.keys()
    onerun = {}

    for run in options:
        for a, b in zip(keys, run):
            onerun[a] = b

        run_name = "_".join([str(x) for x in run])

        amend_env("env.yaml", **onerun)
        onerun = {}

        # Runs EVo
        main("chem.yaml", "multirun.yaml", None)

        # Copies the dgs_output file into a separate file ready to be run again.
        copyfile("Output/dgs_output.csv", f"Output/output_{run_name}.csv")


# multirun(FO2_buffer_START=[1, -2], ATOMIC_C=[150, 550], ATOMIC_H=[200, 500])
