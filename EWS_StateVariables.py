# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 12:02:29 2022

@author: TEJan
"""

"""
EWS - Early Warning Signals
EWS State Variables

@authors: KoenvanLoon & TijmenJanssen
"""

import hourly_configuration as cfg


# Class StateVariable for Variable objects
"""
State variables used in all other EWS Python scripts are stored in the lists below. For each state variable, multiple
variables of these state variables, such as snapshot interval, window size and overlap, and datatype, can be defined
below as to guarantee the use of the same variables over the different EWS Python scripts.

    ! - Note that the 'full' sets of state variables are defined at the end of this file, and if state variables for EWS 
    are added, they also should be included here.
    
Args:
-----

name : string, short name for state variable (usually taken from model).

spatial : bool, selects whether data is spatial or not.

temporal : bool, selects whether data is temporal or not.

snapshot_interval : int, modelled time between saved spatial data.

window_size : int, amount of datapoints over which the early-warning signal is calculated (window).

window_overlap : int, amount of datapoints from previous window taken into the next window.

datatype : string, either 'map' for map-files or 'numpy' for numpy.txt-files.

full_name : string, full name of the state variable.

unit : string, unit of the state variable.
    

Example:
--------

Index Fund :

    INDF = StateVariable('INDF', temporal=True, datatype='numpy', window_size=30, window_overlap=0, 
        datatype='numpy.txt', full_name='Index fund closing value', unit="Dollars ($)")
    
"""

variables_weekly = []
variables_hourly = []


class StateVariable:
    def __init__(self, name, spatial=False, temporal=False, snapshot_interval=cfg.interval_map_snapshots,
                 window_size=3000, window_overlap=0, datatype='map', full_name='', unit='unit'):
        self.name = name
        self.spatial = spatial
        self.temporal = temporal
        self.snapshot_interval = snapshot_interval
        self.window_size = window_size
        self.window_overlap = window_overlap
        self.datatype = datatype
        self.full_name = full_name
        self.unit = unit


# State variables for EWS #

# Maximum interception store
micM = StateVariable('micM', spatial=True, full_name='Maximum interception storage spatial', unit="m")
micA = StateVariable('micA', temporal=True, datatype='numpy', full_name='Maximum interception storage temporal', unit="m")
micL = StateVariable('micL', temporal=True, datatype='numpy', full_name='Maximum interception storage at location', unit="m")

# LAI
laiM = StateVariable('laiM', spatial=True, full_name='LAI spatial', unit="-")
laiA = StateVariable('laiA', temporal=True, datatype='numpy', full_name='LAI temporal', unit="-")
laiL = StateVariable('laiL', temporal=True, datatype='numpy', full_name='LAI at location', unit="-")

# Soil moisture
moiM = StateVariable('moiM', spatial=True, full_name='Soil moisture spatial', unit="-")
moiA = StateVariable('moiA', temporal=True, datatype='numpy', full_name='Soil moisture temporal', unit="-")
moiL = StateVariable('moiL', temporal=True, datatype='numpy', full_name='Soil moisture at location', unit="-")

# Biomass
bioM = StateVariable('bioM', spatial=True, full_name='Biomass spatial', unit="kg m^-2")
bioA = StateVariable('bioA', temporal=True, datatype='numpy', full_name='Biomass temporal', unit="kg m^-2")
bioL = StateVariable('bioL', temporal=True, datatype='numpy', full_name='Biomass at location', unit="kg m^-2")

# Regolith thickness
regM = StateVariable('regM', spatial=True, full_name='Regolith thickness spatial', unit="m")
regA = StateVariable('regA', temporal=True, datatype='numpy', full_name='Regolith thickness temporal', unit="m")
regL = StateVariable('regL', temporal=True, datatype='numpy', full_name='Regolith thickness at location', unit="m")

# DEM
demM = StateVariable('demM', spatial=True, full_name='DEM spatial', unit="m")
demA = StateVariable('demA', temporal=True, datatype='numpy', full_name='DEM temporal', unit="m")
demL = StateVariable('demL', temporal=True, datatype='numpy', full_name='DEM at location', unit="m")

# Discharge
qA = StateVariable('qA', temporal=True, datatype='numpy', full_name='Discharge temporal', unit=1)

# Grazing rate
gA = StateVariable('gA', temporal=True, datatype='numpy', full_name='Grazing rate temporal', unit="kg m^-2 h^-1")

# Growth part
gpM = StateVariable('gpM', spatial=True, full_name='Growth part spatial', unit=1)
gpA = StateVariable('gpA', temporal=True, datatype='numpy', full_name='Growth part temporal', unit=1)

# Grazing part
grM = StateVariable('grM', spatial=True, full_name='Grazing part spatial', unit=1)
grA = StateVariable('grA', temporal=True, datatype='numpy', full_name='Grazing part temporal', unit=1)

# Net growth
grnM = StateVariable('grnM', spatial=True, full_name='Net growth spatial', unit=1)
grnA = StateVariable('grnA', temporal=True, datatype='numpy', full_name='Net growth temporal', unit=1)

# Net deposition
depM = StateVariable('depM', spatial=True, full_name='Net deposition spatial', unit=1)
depA = StateVariable('depA', temporal=True, datatype='numpy', full_name='Net deposition temporal', unit=1)
depL = StateVariable('depL', temporal=True, datatype='numpy', full_name='Net deposition at location', unit=1)

# Net weathering
weaM = StateVariable('weaM', spatial=True, full_name='Net weathering spatial', unit=1)
weaA = StateVariable('weaA', temporal=True, datatype='numpy', full_name='Net weathering temporal', unit=1)
weaL = StateVariable('weaL', temporal=True, datatype='numpy', full_name='Net weathering at location', unit=1)

# Net creep deposition
creM = StateVariable('creM', spatial=True, full_name='Net creep deposition spatial', unit=1)
creA = StateVariable('creA', temporal=True, datatype='numpy', full_name='Net creep deposition temporal', unit=1)
creL = StateVariable('creL', temporal=True, datatype='numpy', full_name='Net creep deposition at location', unit=1)


# Rq
Rq = StateVariable('Rq', temporal=True, datatype='numpy', full_name='Discharge', window_size=3000)
moiM = StateVariable('moiM', spatial=True, full_name='Soil moisture spatial', unit="-")
moiA = StateVariable('moiA', temporal=True, datatype='numpy', full_name='Soil moisture temporal', unit="-")



# Check which variables are present in the configuration and append these to the list of variables

full_set_of_variables_weekly = [micM, micA, micL, laiM, laiA, laiL, moiM, moiA, moiL, bioM, bioA, bioL, regM, regA,
                                regL, demM, demA, demL, qA, gA, gpM, gpA, grM, grA, grnM, grnA, depM, depA, depL, weaM,
                                weaA, weaL, creM, creA, creL]

full_set_of_variables_hourly = [Rq, moiM, moiA]

if cfg.state_variables_for_ews_weekly == 'full':
    variables_weekly = full_set_of_variables_weekly
else:
    for state_variable in cfg.state_variables_for_ews_weekly:
        for variable in full_set_of_variables_weekly:
            if variable.name == state_variable:
                variables_weekly.append(variable)

if cfg.state_variables_for_ews_hourly == 'full':
    variables_hourly = full_set_of_variables_hourly
else:
    for state_variable in cfg.state_variables_for_ews_hourly:
        for variable in full_set_of_variables_hourly:
            if variable.name == state_variable:
                variables_hourly.append(variable)