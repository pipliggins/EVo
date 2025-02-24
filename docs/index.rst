Welcome to EVo's documentation!
=======================================

EVo is a volcanic outgassing model for COHSN elements.

The current implementation has been developed in Python 3 and tested on Windows, Linux and macOS.

Installation/Usage:
*******************

EVo can be used through local installation, either through a CLI or a webapp.

To install locally, EVo must be downloaded from GitHub using
::

   git clone git@github.com:pipliggins/EVo.git

into the project directory where you wish to use EVo. EVo must then be locally pip-installed:
::

   cd EVO
   python -m pip install ".[streamlit]"

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

Alternatively you can interact with EVo using the webapp interface, by running:
::

   cd EVO/webapp
   streamlit run streamlit-app.py


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   evo_doc

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
