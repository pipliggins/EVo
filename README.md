# EVo
This code calculates the speciation and volume of a volcanic gas phase erupting in equilibrium with it's parent magma.
Models can be run to calculate the gas phase in equilibrium with a melt at a single pressure, or the melt can be decompressed from depth rising to the surface as a closed-system case.
Single pressure and decompression can be run for OH, COH, SOH, COHS and COHSN systems, and outputs include:
* Gas phase weight and volume fraction within the system
* Gas phase speciation as mole fraction or weight fraction across H2O, H2, O2, CO2, CO, CH4, SO2, S2, H2S and N2.
* Volatile content of the melt at each pressure
* Melt density
* fO2 of the system
* Total volatile contents as ppm atomic H, O, C, S and N

EVo can be set up using either melt volatile contents, or for a set amount of atomic volatile which is preferable for conducting experiments over a wide range of fO2 values.

**Note: In very rare cases, this version of EVo may result in a run failure if low pressure conditions are being used. If these low pressure runs are needed, an alternative EVo version which contains a multipoint precision module (GMPY2) is available upon emailing me at philippa.liggins@dtc.ox.ac.uk. However, for the vast majority of scenarios, this version is preferable, as it runs faster.**

### Prerequisites

This programme requires Python 3 to run.

If you use python virtual environments, or Anaconda, requirements files (requirements.txt and environment.yml for virtualenv and conda, respectively) can be found in the Data file.

Installation/Usage:
*******************

To install locally, EVo must be downloaded from GitHub using

`git clone git@github.com:pipliggins/EVo.git`

into the project directory where you wish to use EVo. EVo must then be locally pip-installed:
```
cd EVO
python -m pip install -e .
```

From this point, EVo can either be imported into your python scripts as a regular module using
`import evo`

and run using
`evo.main(<chem_file>, <env_file>, <output_options_file>, folder=<output folder name>)`

Or EVo can be run directly from the terminal from inside the `evo` directory:
```
cd EVO/evo
python dgs.py input/chem.yaml input/env.yaml --output-options input/output.yaml
```

The model should run and produce an output file `outputs/dgs_output_*.csv`, and a set of graphs in an 'outputs' folder, if a decompression run has been selected.

### Choosing run options in the env.yaml file

The different run types, model parameters and input values are all set in the env.yaml file. Run types are toggled on or off using True/False, as are the various requirements. You will be warned in the terminal window if the run type you have chosen, and the input parameters selected, do not match. A summary of which parameters are required for which run type will also be given below.

There are multiple run types which can be selected within EVo. At the highest level, a run can either be (1). single pressure, where equilibration between the gas and the melt only occurs at 1 pressure step set using P_START, or (2) a decompression run, where calculations start at P_START and run through pressure steps (the max and min size of which can be set using DP_MAX and DP_MIN respectively) until P_STOP is reached. These two high-level un options can be toggled between using SINGLE_STEP, where True sets up EVo to do a single pressure run, and False asks for decompression.

Within these two high-level run types, 3 options for selecting starting conditions are available:
* Standard: pick a pressure (P_START), a starting gas weight fraction (WgT) and some combination of either current melt volatile contents, or gas fugacities. EVo will calculate the missing variables and either stop at that point for a single pressure run, or continue on in decompression mode. This is the default option, and will be used if both FIND_SATURATION and ATOMIC_MASS_SET are False.
* Volatile saturation: Chosen by switching FIND_SATURATION to True, given only the melt volatile contents and magma fO2, EVo will calculate the volatile saturation pressure and start a run from there. Any values given in P_START and WgT will be ignored; WgT will be set to 1e-8 at the volatile saturation point.
* Atomic set: Chosen by switching ATOMIC_MASS_SET to True (and FIND_SATURATION to False). Given the melt fO2 and the atomic weight fractions of each of the other volatile species (H, +/- C, S & N, set using ATOMIC_H etc.), Evo will calculate both the appropriate distribution of each element across the different species considered, and the volatile saturation point of that composition. Particularly useful for studies where fO2 is varied but the amount of other elements should be held constant. Again, P_START and WgT will be ignored, WgT will be set to 1e-8 at the volatile saturation point.

Listed below is a table of all the options available in the env.yaml file, with a brief explanation of their function.

| Variable name                     | Input units                       | Description                                                                                                                   | Options/example values (default values are given; where multiple options are available, the value in bold is the default) | Required in:                                     |
|-----------------------------------|-----------------------------------|-------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------|
| COMPOSITION                       | N/A - string                      | Magma composition (links to the composition given in chem.yaml)                                                               | **basalt**, phonolite, rhyolite                                                                                             | All runs                                         |
| RUN_TYPE                          | N/A - string                      | Toggle between closed (gas & melt evolve together) and open system (gas is lost after each pressure step) degassing           | **closed**, open                                                                                                            | All runs                                         |
| SINGLE_STEP                       | True/False                        | See above                                                                                                                     | True/**False**                                                                                                              | All runs                                         |
| FIND_SATURATION                   | True/False                        | See above                                                                                                                     | True/**False**                                                                                                              | All runs                                         |
| ATOMIC_MASS_SET                   | True/False                        | See above                                                                                                                     | True/**False**                                                                                                              | All runs                                         |
| GAS_SYS                           | N/A - string                      | Dictate which elements are involved, out of C, O, H, S and N.                                                                 | **OH**, COH, SOH, COHS, COHSN                                                                                               | All runs                                         |
| FE_SYSTEM                         | True/False                        | If true, oxygen is equilibrated with iron in the melt to buffer fO2 change (Fe2O3_{(melt)} <=> 2FeO_{(melt)} + 0.5O2_{(gas)}) | True/**False**                                                                                                              | All runs                                         |
| S_SAT_WARN                        | True/False                        | If True, checks & warns of sulfide saturation, or excess (>10%) sulphate content                                              | **True**/False                                                                                                              | All runs                                         |
| T_START                           | Kelvin                            | System temperature                                                                                                            | 1473.15                                                                                                                   | All runs                                         |
| P_START                           | Bar                               | Starting pressure condition (for standard run)                                                                                | 3000                                                                                                                      | 'Standard' runs                                  |
| P_STOP                            | Bar                               | Final pressure                                                                                                                | 1 (> 0)                                                                                                                   | All runs                                         |
| DP_MIN                            | Bar                               | Minimum size of pressure step before non-convergence is assumed                                                               | 1                                                                                                                         | All runs                                         |
| DP_MAX                            | Bar                               | Maximum size of pressure step                                                                                                 | 100                                                                                                                       | All runs                                         |
| MASS                              | grams                             | Mass of the system                                                                                                            | 100                                                                                                                       | Optional for all                                 |
| WgT                               | Weight fraction                   | Weight fraction of exsolved gas at the starting pressure (for standard run)                                                   | 0.001, (1e-9 -> 0.99)                                                                                                     | 'Standard' runs                                  |
| LOSS_FRAC                         | Weight fraction                   | The fraction of the gas phase which is lost after each pressure step in open system degassing                                 | 0.99 (0<frac<1)                                                                                                           | Open degassing only                              |
| DENSITY_MODEL                     | N/A - string                      | Equation used to calculate melt density                                                                                       | spera2000                                                                                                                 | All runs                                         |
| FO2_MODEL                         | N/A - string                      | Equation to calculate fO2 from melt Fe2/Fe3 ratios                                                                            | **kc1991**, r2013                                                                                                           | All runs                                         |
| FMQ_MODEL                         | N/A - string                      | Rock buffer definitions                                                                                                       | frost1991                                                                                                                 | All runs                                         |
| H2O_MODEL                         | N/A - string                      | Solubility model for water                                                                                                    | burguisser2015                                                                                                            | All runs                                         |
| H2_MODEL                          | N/A - string                      | Solubility model used for H2                                                                                                  | burguisser2015                                                                                                            | All runs                                         |
| C_MODEL                           | N/A - string                      | Solubility model used for CO2/CO3/graphite system                                                                             | **burguisser2015**, eguchi2018                                                                                              | Runs including carbon species                    |
| CO_MODEL                          | N/A - string                      | Solubility model used for CO - recommend setting to None if the fO2 is above IW+1                                             | armstrong2015, **None**                                                                                                     | Runs including carbon species                    |
| CH4_MODEL                         | N/A - string                      | Solubility model used for CH4 - recommend setting to None if the fO2 is above IW+1                                            | ardia2013, **None**                                                                                                         | Runs including carbon species                    |
| SULFIDE_CAPACITY                  | N/A - string                      | Sulfide capacity law used for S2- solubility                                                                                  | oneill2002, **oniell2020**                                                                                                  | Runs including sulfur species                    |
| SULFATE_CAPACITY                  | N/A - string                      | Sulphate capacity used for SO4 solubility                                                                                     | nash2019                                                                                                                  | Runs including sulfur species                    |
| SCSS                              | N/A - string                      | Model used to calculate sulfide capacity at sulfide saturation                                                                | liu2007                                                                                                                   | Runs including sulfur species                    |
| N_MODEL                           | N/A - string                      | Solubility model used for N2                                                                                                  | libourel2003                                                                                                              | Runs including nitrogen                          |
| FO2_BUFFER_SET                    | True/False                        | Toggle option to set the initial fO2 using a value relative to a rock buffer (either IW, FMQ or NNO)                          | True/**False**                                                                                                              | Initial volatile content option, see Table 2     |
| FO2_BUFFER                        | N/A - string                      | Set which rock buffer is used if the above option is set to True                                                              | IW, **FMQ**, NNO                                                                                                            | If FO2_BUFFER_SET = True                           |
| FO2_BUFFER_START                  | log units relative to rock buffer | Starting fO2 value relative to rock buffer if above option(s) are True                                                        | 0                                                                                                                         | If FO2_BUFFER_SET = True                           |
| FO2_SET/FO2_START                 | True/False, bar                   | Toggle option to set starting fO2 condition as an absolute oxygen fugacity value in bar, and starting value.                  | True/**False**, 0                                                                                                           | Initial volatile content option, see Table 2     |
| ATOMIC_H                          | ppm by total system weight        | Total weight fraction of H in the system (melt + gas)                                                                         | 550                                                                                                                       | When SET_ATOMIC_MASS = True only                   |
| ATOMIC_C                          | ppm by total system weight        | Total weight fraction of C in the system (melt + gas)                                                                         | 200                                                                                                                       | When SET_ATOMIC_MASS = True only, if C is included |
| ATOMIC_S                          | ppm by total system weight        | Total weight fraction of S in the system (melt + gas)                                                                         | 4000                                                                                                                      | When SET_ATOMIC_MASS = True only, if S is included |
| ATOMIC_N                          | ppm by total system weight        | Total weight fraction of N in the system (melt + gas)                                                                         | 10                                                                                                                        | When SET_ATOMIC_MASS = True only, if N is included |
| FH2_SET/ FH2_START                 | True/False, bar                   | Toggle using hydrogen fugacity as a starting condition on/off, and set the starting value.                                    | True/**False**, 0.24                                                                                                        | Initial volatile content option, see Table 2     |
| FH2O_SET/ FH2O_START               | True/False, bar                   | Toggle using water fugacity as a starting condition on/off, and set the starting value.                                       | True/**False**, 1000                                                                                                        | Initial volatile content option, see Table 2     |
| FCO2_SET/ FCO2_START               | True/False, bar                   | Toggle using carbon dioxide fugacity as a starting condition on/off, and set the starting value.                              | True/**False**, 0.01                                                                                                        | Initial volatile content option, see Table 2     |
| WTH2O_SET/ WTH2O_START             | True/False, weight fraction       | Toggle using water content in the melt as a starting condition on/off, and set the starting value.                            | True/**False**, 0.03                                                                                                        | Initial volatile content option, see Table 2     |
| WTCO2_SET/ WTCO2_START             | True/False, weight fraction       | Toggle using carbon dioxide content in the melt as a starting condition on/off, and set the starting value.                   | True/**False**, 0.01                                                                                                        | Initial volatile content option, see Table 2     |
| SULFUR_SET/ SULFUR_START           | True/False, weight fraction       | Toggle using sulfur content in the melt as a starting condition on/off, and set the starting value.                           | True/**False**, 0.001                                                                                                       | Initial volatile content option, see Table 2     |
| NITROGEN_SET/ NITROGEN_START       | True/False, weight fraction       | Toggle using nitrogen content in the melt as a starting condition on/off, and set the starting value.                         | True/**False**, 0.01                                                                                                        | Initial volatile content option, see Table 2     |
| GRAPHITE_SATURATED/ GRAPHITE_START | True/False, weight fraction       | Toggle whether the melt is initially graphite saturated on/off, and set the starting value.                                   | True/**False**, 0                                                                                                           | Initial volatile content option, see Table 2     |

The various combinations of initial conditions available for each run type is tabulated below. In all cases unless stated otherwise, the starting fO2 must be set by selecting one from 'FO2_START' (use absolute fO2 value) or 'FO2_BUFFER_START' (use fO2 relative to buffer).

| Elements involved | Run type      | Non-optional values | (1) Initial volatile requirements - gas fugacities                     | (2) Initial volatile requirements - melt weight fractions                     | (3) Initial volatile requirements - total atomic masses                                                                          |
|-------------------|---------------|---------------------|------------------------------------------------------------------------|-------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------|
| OH                | 'Standard'    | P_START, WgT        | FH2, but only WITHOUT fO2 being set.                                   | Not an option.                                                                | Not an option - although it should be possible? set atomic masses but don't find saturation, just at current P (WgT is unknown). |
| OH                | Saturation    |                     | Not an option.                                                         | WTH2O                                                                         | Not an option                                                                                                                    |
| OH                | Atomic Values |                     | Not an option.                                                         | Not an option.                                                                | 'ATOMIC_H'.                                                                                                                      |
| COH               | 'Standard'    | P_START, WgT        | FH2, but only WITHOUT fO2 being set. Plus one from either FH2O or FCO2 | Pick one from WTH2O, WTCO2 or GRAPHITE_SATURATED                              | Not an option                                                                                                                    |
| COH               | Saturation    |                     | Not an option.                                                         | WTH2O and either WTCO2 or GRAPHITE_SATURATED                                  | Not an option                                                                                                                    |
| COH               | Atomic Values |                     | Not an option.                                                         | Not an option.                                                                | 'ATOMIC_H' and 'ATOMIC_C'.                                                                                                       |
| SOH               | 'Standard'    | P_START, WgT        | Not an option.                                                         | Pick one from either WTH2O or SULFUR_START                                    | Not an option                                                                                                                    |
| SOH               | Saturation    |                     | Not an option.                                                         | WTH2O and SULFUR_START                                                        | Not an option                                                                                                                    |
| SOH               | Atomic Values |                     | Not an option.                                                         | Not an option.                                                                | 'ATOMIC_H' and 'ATOMIC_S'.                                                                                                       |
| COHS              | 'Standard'    | P_START, WgT        | Not an option                                                          | Pick 2 from WTH2O, SULFUR_START, WTCO2 and GRAPHITE_SATURATED                 | Not an option                                                                                                                    |
| COHS              | Saturation    |                     | Not an option.                                                         | WTH2O, SULFUR_START and either WTCO2 or GRAPHITE_SATURATED                    | Not an option                                                                                                                    |
| COHS              | Atomic Values |                     | Not an option.                                                         | Not an option.                                                                | ATOMIC_H', 'ATOMIC_C' and 'ATOMIC_S'.                                                                                            |
| COHSN             | 'Standard'    | P_START, WgT        | Not an option.                                                         | Pick 3 from WTH2O, SULFUR_START, NITROGEN_START, WTCO2 and GRAPHITE_SATURATED | Not an option                                                                                                                    |
| COHSN             | Saturation    |                     | Not an option.                                                         | WTH2O, SULFUR_START, NITROGEN_START and either WTCO2 or GRAPHITE_SATURATED    | Not an option                                                                                                                    |
| COHSN             | Atomic Values |                     | Not an option.                                                         | Not an option.                                                                | 'ATOMIC_H', 'ATOMIC_C', 'ATOMIC_S' and 'ATOMIC_N'.                                                                               |


## Missing from EVo

This repository is still in development, and as such warnings about actions which are currently unsupported may not exist.

The following HAVE NOT been implemented yet and so will not work, even if no error is generated and the code appears to run - it's not doing what you think it's doing!
* Gas-only chemistry

## Reading the Output file

This table summarises the headers and meanings in the output CSV file.

| Name in output file  | Symbol | Units    | Description                                                                                      |
|----------------------|--------|----------|--------------------------------------------------------------------------------------------------|
| P              |    P    | bar      | Pressure                                                                                         |
| NNO/FMQ/IW           |        |          | Oxygen fugacity relative to the chosen buffer, Ni-NiO, Fayalite-Magnetite-Quartz or Iron-Wursite |
| fO2            |        | bar      | Raw oxygen fugacity                                                                               |
| F              |        | ratio      | mole fraction ratio of Fe2O3/FeO                                                                     |
| rho_bulk           |        | kg/m^3   | Density of the system (melt+gas)                                                                 |
| rho_melt           |        | kg/m^3   | Density of the silicate melt                                                                     |
| Exsol_vol%          |        | % | Gas volume fraction (of system)                                                                  |
| Gas_wt              |    WgT    | % | Gas weight fraction (of system)                                                                  |
| mol_mass      |        | kg/mol   | Mean molecular mass of gas phase                                                                 |
| mH2O, mH2, mO2...    |        | fraction | Species mole fraction in gas                                                                    |
| WH2O, WH2, WO2...    |        | fraction | Species weight fraction in gas                                                                   |
| fH2O, fH2...         |        | bar      | Species fugacities                                                                               |
| mCO2/CO, mCO2/H2O... |        |          | Molar ratios of gas species                                                                      |
| H2O_melt, H2_melt...   |        | % | Weight percent of species in the melt                                                            |
| S_melt                |        | % | Weight percent of S in the melt (sum of S dissolved as S2- and S6+)                              |
| tot_H, tot_C...      |        | ppm      | Total weight fraction of elements                                                                |
| tot_O_gas      |        | ppm      | Total weight fraction of elemental O in the volatile phases only, i.e. not counting O stored as FeO or Fe2O3                                                                |
