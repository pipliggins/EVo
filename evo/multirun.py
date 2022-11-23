
from itertools import product
import numpy as np
import os
from shutil import copyfile
import yaml

from evo.dgs import main

def amend_env(file, **kwargs):
    """
    Duplicates the env.yaml file, then edits with the gridsearch values relevant
    to the specific run.
    """

    with open(file, 'r') as f:
        env_doc = yaml.full_load(f)

    for param, val in kwargs.items():
        if type(val) == np.float64:
            env_doc[param] = float(val)
        else:
            env_doc[param] = val

    with open('multirun.yaml', "w") as f:
        yaml.dump(env_doc, f)

def multirun(**kwargs):
    """
    Gridsearch over all the kwargs given (as lists or arrays), having setup variables you don't want to search over
    within the env.yaml already.

    Each run result will be saved as a separate output file in Outputs, named using the
    values used for that gridsearch.
    """
    
    # delete any pre-existing multirun setup files
    if os.path.exists('multirun.yaml'):
        os.remove('multirun.yaml')
    
    #creates an iterable object where each 'column' is in the order given here.
    options = product(*kwargs.values())
    keys=kwargs.keys()
    onerun={}    
    
    for run in options:
        for a, b in zip(keys, run):
            onerun[a] = b
        
        run_name = '_'.join([str(x) for x in run])
        
        amend_env('env.yaml', **onerun)
        onerun = {}

        # Runs EVo
        main('chem.yaml', 'multirun.yaml', None)

        # Copies the dgs_output file into a separate file ready to be run again.
        copyfile('Output/dgs_output.csv', f'Output/output_{run_name}.csv')

# multirun(FO2_buffer_START=[1, -2], ATOMIC_C=[150, 550], ATOMIC_H = [200, 500])





