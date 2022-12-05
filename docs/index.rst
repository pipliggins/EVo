.. evo-outgas documentation master file, created by
   sphinx-quickstart on Mon Dec  5 10:50:53 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to EVo's documentation!
=======================================

EVo is a volcanic outgassing model for COHSN elements.

The current implementation has been developed in Python 3 and tested on Windows and Linux.

Installation/Usage:
*******************

EVo can be used either through local installation, or using a web-app at `<https://evo.pipliggins.co.uk>`_.

To install locally, EVo must be downloaded from GitHub using
::
   
   git clone -b streamlit-restructure --single-branch git@github.com:pipliggins/EVo.git

into the project directory where you wish to use EVo. EVo must then be locally pip-installed:
::

   cd EVO
   python -m pip install -e evo/

From this point, EVo can either be imported into your python scripts as a regular module using
::

   install evo

and run using
::
   
   evo.main('chem_file', 'env_file', 'output_options_file')

Or EVo can be run directly from the terminal from inside the `evo` directory:
::
   
   cd EVO/evo
   python dgs.py input/chem.yaml input/env.yaml --output input/output.yaml



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   evo_doc

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
