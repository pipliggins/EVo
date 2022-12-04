"""Contains code to enable a web interface for EVo, using StreamLit."""

import evo
import streamlit as st
from pathlib import Path
import pandas as pd
import numpy as np
import yaml
import ruamel.yaml

ryaml = ruamel.yaml.YAML()

def amend_chem(file, chem_dict):
    """
    Duplicates the chem.yaml file, then edits with the gridsearch values relevant
    to the specific run.
    """

    with open(file, 'r') as f:
        chem_doc = ryaml.load(f)

    for param, val in chem_dict.items():
        if type(val) == np.float64:
            chem_doc[param] = float(val)
        elif val == 0.0:
            pass
        else:
            chem_doc[param] = val

    with open('chem.yaml', "w") as f:
        ryaml.dump(chem_doc, f)

def amend_env(file, env_dict):
    """
    Duplicates the env.yaml file, then edits with the gridsearch values relevant
    to the specific run.
    """

    with open(file, 'r') as f:
        env_doc = ryaml.load(f)

    for param, val in env_dict.items():
        print(param, type(val), val)
        if type(val) == np.float64: # This forces the scientific notation to be output correctly to be read as a float
            env_doc[param] = float(val)
        else:
            env_doc[param] = val

    with open('env.yaml', "w") as f:
        ryaml.dump(env_doc, f)

oxides = {}
env_options = {}
default_laws = {'FO2_MODEL':'kc1991', 'H2_MODEL':'gaillard2003', 'C_MODEL':'burguisser2015', 'CO_MODEL':'armstrong2015', 'CH4_MODEL':'ardia2013', 'SULFIDE_CAPACITY':'oneill2020'}
volatile_params = {}
models = {}

st. set_page_config(layout="wide")

# # Use this to control which things are visible.
## i.e. use 'disabled' to grey out options which are not valid based on the run options.
# https://docs.streamlit.io/library/api-reference/widgets/st.radio
# if "startType" not in st.session_state:
#     # st.session_state.visibility = "visible"
#     # st.session_state.disabled = False
#     # st.session_state.horizontal = False

if "default_sollaws" not in st.session_state:
    st.session_state.default_sollaws = True

st.title('EVo: A volcanic degassing model')

st.text("This is a web interface for the EVo volcanic degassing model.\nUse the input options below to select setup options and load your input options by clicking 'Load data'.\nThen press 'Run EVo' to generate a dataframe of results.")

col1, col2, col3 = st.columns(3)

with col1:
   st.header("Major Oxide Composition")
   sio2 = st.number_input('SiO2', min_value=0.0, max_value=100.0, value=47.95)
   tio2 = st.number_input('TiO2', min_value=0.0, max_value=100.0, value=1.67)
   al2o3 = st.number_input('Al2O3', min_value=0.0, max_value=100.0, value=17.32)
   feo = st.number_input('FeO', min_value=0.0, max_value=100.0, value=10.24)
   # set this to zero and greyed out if fo2 is set
   fe2o3 = st.number_input('Fe2O3', min_value=0.0, max_value=100.0, value=0.0)
   mno = st.number_input('MnO', min_value=0.0, max_value=100.0, value=0.17)
   mgo = st.number_input('MgO', min_value=0.0, max_value=100.0, value=5.76)
   cao = st.number_input('CaO', min_value=0.0, max_value=100.0, value=10.93)
   na2o = st.number_input('Na2O', min_value=0.0, max_value=100.0, value=3.45)
   k2o = st.number_input('K2O', min_value=0.0, max_value=100.0, value=1.99)
   p2o5 = st.number_input('P2O5', min_value=0.0, max_value=100.0, value=0.51)
   nio = st.number_input('NiO', min_value=0.0, max_value=100.0, value=0.0)
   cr2o3 = st.number_input('Cr2O3', min_value=0.0, max_value=100.0, value=0.0)

   oxides = {'SIO2':sio2, 'TIO2':tio2, 'AL2O3':al2o3, 'FEO':feo, 'FE2O3':fe2o3,
   'MNO':mno, 'MGO':mgo, 'CAO':cao, 'NA2O':na2o, 'K2O':k2o, 'P2O5':p2o5, 'NIO':nio, 'CR2O3':cr2o3}

with col2:
   st.header("Model setup options")
   rock_type = st.selectbox("Magma type", ('basalt', 'phonolite', 'rhyolite'))
   run_type = st.selectbox("Run Type - Closed System or Open System", ('closed', 'open'))
   if run_type == 'Open system':
    loss_frac = st.number_input('Fraction of gas phase lost at each P step during open system degassing', min_value=0.0, max_value=0.9999, value=0.99)
   else:
    loss_frac = 0.0

   single_step = st.checkbox('Run for a single pressure?')
   
   env_options.update({'COMPOSITION':rock_type,
                  'RUN_TYPE':run_type,
                  'SINGLE_STEP':single_step})

   saturation = st.radio("Choose how to find the starting pressure:", ('Set  manually', 'Find volatile saturation pressure', 'Set volatile element mass, and find volatile saturation pressure'))
   if saturation == 'Set  manually':
    env_options['FIND_SATURATION'] = False
    env_options['ATOMIC_MASS_SET'] = False
   elif saturation == 'Find volatile saturation pressure':
    env_options['FIND_SATURATION'] = True
    env_options['ATOMIC_MASS_SET'] = False
   else:
    env_options['FIND_SATURATION'] = False
    env_options['ATOMIC_MASS_SET'] = True

   # Eventually use this to control which volatile input options are shown
   gas_sys = st.selectbox('Select the volatile element system:', ('OH', 'COH', 'SOH', 'COHS', 'COHSN'))   
   fe_equil = st.checkbox("Allow redox equilibration between gas and melt?")
   s_sat_warn = st.checkbox("Stop the run if sulfide saturation is reached? (EVo is not suitable for use under sulfide-saturated conditions)")

   st.subheader("Starting conditions")
   t_start = st.number_input('Temperature (K)', min_value=600.0, max_value=2000.0, value=1473.15)
   # work out how to grey this out if you don't need to select a starting pressure
   p_start = st.number_input('Starting Pressure (bar)', min_value=1.0, max_value=10000.0, value=3000.0)
   p_stop = st.number_input('Final Pressure (bar)', min_value=1e-5, max_value=10000.0, value=1.0)
   dp_max = st.number_input('Maximum pressure step size (bar)', min_value=1.0, max_value=100.0, value=100.0)
   dp_min = st.number_input('Minimum pressure step size (bar)', min_value=1e-5, max_value=10.0, value=0.1)
   # grey this one out too
   wgt = st.number_input('Gas weight fraction at starting pressure', min_value=0.000000001, max_value=1.0, value=0.00001)
   mass = st.number_input('Mass of the magma packet (g)', min_value=0.0, value=100.0)

   env_options.update({'GAS_SYS':gas_sys,
                       'FE_SYSTEM':fe_equil,
                       'S_SAT_WARN':s_sat_warn,
                       'T_START':t_start,
                       'P_START':p_start,
                       'P_STOP':p_stop,
                       'DP_MAX':dp_max,
                       'DP_MIN':dp_min,
                       'WgT':wgt,
                       'MASS':mass})
   
   sol_laws = st.checkbox("Use default solubility laws? (if unchecked, dropdown options will be offered)", value = True, key='default_sollaws')
   # These will be greyed out if default solubility laws are chosen
   
   fo2_model = st.selectbox('Ferric/Ferrous -> fO2 conversions model', ('kc1991', 'righter2015'), disabled=st.session_state.default_sollaws)
   h2o_model = st.selectbox('H2O solubility law', ('burguisser2015',), disabled=st.session_state.default_sollaws)
   h2_model = st.selectbox('H2 solubility law', ('gaillard2003', 'burguisser2015'), disabled=st.session_state.default_sollaws)
   co2_model = st.selectbox('CO2 solubility law', ('burguisser2015', 'eguchi2018'), disabled=st.session_state.default_sollaws)
   co_model = st.selectbox('CO solubility law', ('armstrong2015', 'None'), disabled=st.session_state.default_sollaws)
   ch4_model = st.selectbox('CH4 solubility law', ('ardia2013', 'None'), disabled=st.session_state.default_sollaws)
   sulfide_model = st.selectbox('sulfide capacity law', ('oneill2020', 'oneill2002'), disabled=st.session_state.default_sollaws)
   n2_model = st.selectbox('N2 solubility law', ('libourel2003',), disabled=st.session_state.default_sollaws)

   if sol_laws == True:
    models = default_laws
   else:
    models = {'FO2_MODEL':fo2_model,
              'H2_MODEL':h2_model,
              'C_MODEL':co2_model,
              'CO_MODEL':co_model,
              'CH4_MODEL': ch4_model,
              'SULFIDE_CAPACITY': sulfide_model}

with col3:
    st.header("Volatile content")
    buffer_type = st.selectbox('Select the mineral buffer for the oxygen fugacity:', ('FMQ', 'IW', 'NNO'))
    fo2 = st.number_input('Starting fO2 relative to the mineral buffer:', value = 0.0)
    
    volatile_params['FO2_buffer'] = buffer_type
    volatile_params['FO2_buffer_START'] = fo2

    vol_setup = st.radio("Choose a way to set the system volatile content:", ("Set gas phase", "Set melt content", "Set atomic weight fractions"))
    
    if vol_setup == 'Set gas phase':
        st.subheader('Set gas phase fugacities')
        fh2_set = st.checkbox('Set fH2 (bar)', value=False)
        fh2_start = st.number_input('fh2', value = 0.0, label_visibility="collapsed")

        fh2o_set = st.checkbox('Set fH2O (bar)', value=False)
        fh2o_start = st.number_input('fh2o', value = 0.0, label_visibility="collapsed")

        fco2_set = st.checkbox('Set fCO2 (bar)', value=False)
        fco2_start = st.number_input('fco2', value = 0.0, label_visibility="collapsed")

        volatile_params.update({'FH2_SET':fh2_set,
                           'FH2_START':fh2_start,
                           'FH2O_SET':fh2_set,
                           'FH2O_START':fh2o_start,
                           'FCO2_SET':fco2_set,
                           'FCO2_START':fco2_start})
        
        for param in ['wth2o_set', 'wtco2_set', 'sulfur_set', 'nitrogen_set', 'graphite_saturated']:
            volatile_params[param.upper()] = False

        for param in ['wth2o_start', 'wtco2_start', 'sulfur_start', 'nitrogen_start', 'graphite_start', 'atomic_h', 'atomic_c', 'atomic_s', 'atomic_n']:
            volatile_params[param.upper()] = 0.0       

    elif vol_setup == 'Set melt content':
        st.subheader('Set melt volatile content')
        wth2o_set = st.checkbox('Set melt H2O content (wt fraction)', value=False)
        wth2o_start = st.number_input('wth2o', value = 0.0, label_visibility="collapsed")

        wtco2_set = st.checkbox('Set melt CO2 content (wt fraction)', value=False)
        wtco2_start = st.number_input('wtco2', value = 0.0, label_visibility="collapsed")

        sulfur_set = st.checkbox('Set melt sulfur content (wt fraction)', value=False)
        sulfur_start = st.number_input('sulfur', value = 0.0, label_visibility="collapsed")

        nitrogen_set = st.checkbox('Set melt N2 content (wt fraction)', value=False)
        nitrogen_start = st.number_input('nitrogen', value = 0.0, label_visibility="collapsed")

        graphite_saturated = st.checkbox('Is the melt graphite saturated?', value=False)
        graphite_start = st.number_input('Weight fraction of graphite in the melt', value = 0.0)

        volatile_params.update({'WTH2O_SET':wth2o_set,
                           'WTH2O_START':wth2o_start,
                           'WTCO2_SET':wtco2_set,
                           'WTCO2_START':wtco2_start,
                           'SULFUR_SET':sulfur_set,
                           'SULFUR_START':sulfur_start,
                           'NITROGEN_SET':nitrogen_set,
                           'NITROGEN_START':nitrogen_start,
                           'GRAPHITE_SATURATED':graphite_saturated,
                           'GRAPHITE_START':graphite_start})
        
        for param in ['fh2_set', 'fh2o_set', 'fco2_set']:
            volatile_params[param.upper()] = False

        for param in ['fh2_start', 'fh2o_start', 'fco2_start', 'atomic_h', 'atomic_c', 'atomic_s', 'atomic_n']:
            volatile_params[param.upper()] = 0.0

    elif vol_setup == 'Set atomic weight fractions':
        st.subheader('Set the mass fraction of each volatile element in the system (ppm)')

        atomic_h = st.number_input('Weight fraction of H in the system', value = 0.0)
        atomic_c = st.number_input('Weight fraction of C in the system', value = 0.0)
        atomic_s = st.number_input('Weight fraction of S in the system', value = 0.0)
        atomic_n = st.number_input('Weight fraction of N in the system', value = 0.0)

        volatile_params.update({'ATOMIC_H':atomic_h,
                           'ATOMIC_C':atomic_c,
                           'ATOMIC_S':atomic_s,
                           'ATOMIC_N':atomic_n})
        
        for param in ['fh2_set', 'fh2o_set', 'fco2_set', 'wth2o_set', 'wtco2_set', 'sulfur_set', 'nitrogen_set', 'graphite_saturated']:
            volatile_params[param.upper()] = False

        for param in ['fh2_start', 'fh2o_start', 'fco2_start', 'wth2o_start', 'wtco2_start', 'sulfur_start', 'nitrogen_start', 'graphite_start']:
            volatile_params[param.upper()] = 0.0

col4, col5 = st.columns(2)

with col4:
    load_data = st.button("Load input data")
    if load_data:
        amend_chem('chem.yaml', oxides)
        env_params = env_options | models | volatile_params
        amend_env('env.yaml', env_params) # update to get all options

with col5:
    # Want to grey this out as an option until 'load data' has been pressed...
    run_evo = st.button("Run EVo")
    if run_evo:
        filepath = filepath = Path(__file__).parent
        df = evo.main(filepath / 'chem.yaml', filepath / 'env.yaml', filepath / 'output.yaml')

st.header("Results:")
if run_evo:
    st.dataframe(df)


