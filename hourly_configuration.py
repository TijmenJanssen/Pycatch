# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 15:33:38 2022

@author: TEJan
"""


import pathlib

## Time steps for hourly/weekly model runs
# - TODO replace numberOfTimeSteps for number_of_timesteps_hourly
# Define number of hourlyy time steps to run
number_of_timesteps_hourly = 3000  # in hours  AFRONDEN OP 100 TAL, ZODAT HET DEELBAAR IS DOOR WINDOWSIZE IN STATEVARIABLES 
# number_of_timesteps_hourly = 8760 # ~1y in hours

#Step size between h-model runs
#stepsizeSnapshots = 5200 mag er later uit

cutoff = False

#Loop over snapshots
stepsInShift = 100
stepsTotal = 100
# For 1 h-model run as test, use number_of_timesteps_weekly 


# Define number of weekly time steps to run
number_of_timesteps_weekly = 104000
#number_of_timesteps_weekly = 5200  # in weeks
# number_of_timesteps_weekly = 104000 # ~2000y in weeks

## Rate of forcing (grazing)
# Define the fraction of total time at which no grazing occurs at the beginning (baseline for initial state)
rel_start_grazing = 1/16
# rel_start_grazing = 0
# Define the total increase in grazing rate
tot_increase_grazing = 0.0003
# Select whether grazing rate returns to the initial value after the halfway point
# - note that this halfway point occurs on (1 - rel_start) * total time / 2
return_ini_grazing = True

## Number of Monte Carlo (MC) runs
# Define the number of MC samples or particles, results of realizations are written to the folder(s) 1, 2, ...
nrOfSamples = 1

## Particle filtering
# When True, a particle filtering run is done, usually False for first time users
filtering = False

## Create realizations
# Selects whether a single, given value is used for a number of parameters or whether a realization for that parameter
# is drawn randomly. Usually False for first time users.
createRealizations = False

## Calculate upstream totals
# Selects whether upstream totals are calculated (accuflux) in the subsurfacewateronelayer and interceptionuptomaxstore
# modules. May be needed for some reports and possibly budget checks (if one needs these). For normal use, this is set
# to False.
calculateUpstreamTotals = False





##############

# Parameters 

##############

## Reporting of variables
# Selects which set of variables are reported, either 'full' or 'filtering'. These are passed to the class of a
# component where the variables that can be reported can be defined.
# setOfVariablesToReport = 'full'
# setOfVariablesToReport = 'filtering'
setOfVariablesToReport = 'None'

# Option to call the methods that change the geomorphology in the weekly model, typically True
changeGeomorphology = True

# Option to fix both the regolith and the vegetation in the weekly model, typically False
fixedStates = False

## Timesteps to report variables for both hourly and weekly

# Saving map files which can be used for spatial EWS or hourly model
"Weekly only"  # TODO - Might be worthwhile to implement this into the hourly model
map_data = True
mean_timeseries_data = True
interval_map_snapshots = 100

# definition for components were all timesteps should be reported
timesteps_to_report_all_hourly = list(range(1, number_of_timesteps_hourly + 1, 1))
timesteps_to_report_all_weekly = list(range(interval_map_snapshots, number_of_timesteps_weekly + 1,
                                            interval_map_snapshots))
# timeStepsToReportAll = list(range(1, number_of_timesteps_hourly + 1, 1))

# Used for discharge (hourly model only)
timeStepsToReportRqs = list(range(20, number_of_timesteps_hourly + 1, 20))

# definition for components were a subset of timesteps should be reported
timesteps_to_report_some_hourly = list(range(100, number_of_timesteps_hourly + 1, 100))
timesteps_to_report_some_weekly = list(range(100, number_of_timesteps_weekly + 1, interval_map_snapshots))
# timeStepsToReportSome = list(range(100, number_of_timesteps_hourly + 1, 100))


## State variables for which to calculate Early Warning Signals (both hourly & weekly)
# - TODO change state_variables_for_ews to state_variables_for_ews_hourly
# - TODO check the 'full' list in EWS_StateVariables.py
#state_variables_for_ews_hourly = ['Gs']
state_variables_for_ews_hourly = 'full'
state_variables_for_ews_weekly = ['bioM', 'bioA', 'moiM', 'moiA', "regA"]
# state_variables_for_ews_weekly = 'full'

## Reporting for the model components (both hourly and weekly)
if setOfVariablesToReport == 'full':
    interception_report_rasters = ["Vo", "Vi", "Vgf", "Vms", "Vs"]
    #   reports of totals (Vot) only make sense if calculateUpstreamTotals is True
    infiltration_report_rasters_weekly = ["Ii", "Is", "Iks"]
    infiltration_report_rasters = ["Ii", "Ij", "Is", "Iks"]  # TODO - might want to rename this to ""_hourly, as above
    runoff_report_rasters = ["Rq", "Rqs"]
    subsurface_report_rasters = ["Gs", "Go"]    
    #   reports of totals (Gxt, Got) only make sense if calculateUpstreamTotals is True
    shading_report_rasters = ["Mfs", "Msc", "Msh"]
    surfacestore_report_rasters = ["Ss", "Sc"]
    rainfalleventsfromgammadistribution_report_rasters = ["Pf"]
    exchange_report_rasters = ["Xrc"]
    soilwashMMF_report_rasters = ["Wde", "Wdm", "Wfl"]
    regolith_report_rasters = ["Ast"]
    bedrockweathering_report_rasters = ["Cwe"]
    evapotrans_report_rasters = ["Ep", "Epc"]
    evapotranspirationsimple_report_rasters = ["Ep", "Ea"]
    biomassmodifiedmay_report_rasters = ["Xs"]
    baselevel_report_rasters = ["Ll"]
    creep_report_rasters = ["Ds"]
    randomparameters_report_rasters = ["RPic", "RPks", "RPrt", "RPsc", "RPmm"]
elif setOfVariablesToReport == 'filtering':
    interception_report_rasters = ["Vs", "Vgf"]
    #   reports of totals (Vot) only make sense if calculateUpstreamTotals is True
    infiltration_report_rasters_weekly = []
    infiltration_report_rasters = []
    runoff_report_rasters = ["Rq"]
    subsurface_report_rasters = []
    #   reports of totals (Gxt, Got) only make sense if calculateUpstreamTotals is True
    shading_report_rasters = []
    surfacestore_report_rasters = []
    rainfalleventsfromgammadistribution_report_rasters = []
    exchange_report_rasters = []
    soilwashMMF_report_rasters = []
    regolith_report_rasters = []
    bedrockweathering_report_rasters = []
    evapotrans_report_rasters = []
    evapotranspirationsimple_report_rasters = []
    biomassmodifiedmay_report_rasters = []
    baselevel_report_rasters = []
    creep_report_rasters = []
    randomparameters_report_rasters = []
elif setOfVariablesToReport == 'None':
    interception_report_rasters = []
    #   reports of totals (Vot) only make sense if calculateUpstreamTotals is True
    infiltration_report_rasters_weekly = []
    infiltration_report_rasters = []
    runoff_report_rasters = []
    subsurface_report_rasters = []
    #   reports of totals (Gxt, Got) only make sense if calculateUpstreamTotals is True
    shading_report_rasters = []
    surfacestore_report_rasters = []
    rainfalleventsfromgammadistribution_report_rasters = []
    exchange_report_rasters = []
    soilwashMMF_report_rasters = []
    regolith_report_rasters = []
    bedrockweathering_report_rasters = []
    evapotrans_report_rasters = []
    evapotranspirationsimple_report_rasters = []
    biomassmodifiedmay_report_rasters = []
    baselevel_report_rasters = []
    creep_report_rasters = []
    randomparameters_report_rasters = []

# TODO - Put the parts below above the reporting for the model components in their respective place, i.e. hour/week
#  model parts above the EWS stuff, after removal of unnecessary statements.

######################
# Hourly inputs only #
######################

# folder with input files (maps, timeseries)
inputFolder = "inputs_from_weekly"

# switch to report for locations as small numpy files
# mainly used for particle filtering
doReportComponentsDynamicAsNumpy = False

# switch to swap parameter values between two catchments
# first time users will need to set this to False
swapCatchments = False

# when True, one can read a set of parameters for all Monte Carlo realizations
# from disk (e.g. representing probability distributions from a calibration)
# first time users should have a False here
readDistributionOfParametersFromDisk = False


with_shading = True

if with_shading is False:
    fractionReceivedValue = 1.0
    fractionReceivedFlatSurfaceValue = 1.0


probabilityOfRainstorm=0.4
rainstormDuration=2
expectedRainfallIntensity=0.002
gammaShapeParameter=100




################
# model inputs #
################

# general ########

 
# set 
cloneString = str(pathlib.Path(inputFolder,"clone.map"))


# report locations, i.e. outflow points, for instance, at the outlet
locations = str(pathlib.Path(inputFolder, "clone.map"))


# meteorology #######

#Deze mogen er ook uit als we evaporatie met simple doen
# airTemperatureDetermString = str(pathlib.Path(inputFolder, "airTemperatureArnaJulAugSep0506.tss"))
# relativeHumidityDetermString = str(pathlib.Path(inputFolder, "relativeHumidityArnasJulAugSep0506.tss"))
# incomingShortwaveRadiationFlatSurfaceString = str(pathlib.Path(inputFolder, "incomingShortwaveRadiationArnasJulAugSep0506.tss"))
# windVelocityDetermString = str(pathlib.Path(inputFolder, "windVelocityArnasJulAugSep0506.tss"))
# elevationAboveSeaLevelOfMeteoStationValue = 900.0


# interception #######
# maximumInterceptionCapacityValue = mogelijk relevant voor w-model
# leafAreaIndexValue = mogelijk relevant voor w-model


# surface storage ######
maxSurfaceStoreValue = 0.0001 


# infiltration #######


# regolith geometry ########
# regolithThickness = mogelijk relevant voor w-model 


#TIJMEN reeds verandert in w-model parameters
# 'groundwater' (saturated flow) ##########
saturatedConductivityMetrePerDayValue = 12.5 
limitingPointFractionValue = 0.05
mergeWiltingPointFractionFSValue = 0.019
fieldCapacityFractionValue = 0.22

# green and ampt
# ksatValue = mogelijk relevant voor w-model
initialSoilMoistureFractionCFG = 0.22 # (= fieldCapacityFractionValue)
soilPorosityFractionValue = 0.43

# evapotranspiration ###########

# penman
# Physiclly rrelevant, a remnant of old model structure
multiplierMaxStomatalConductanceValue = 1.0


# # Irrelevant bij simple evaporation
# albedoSoil = 0.3
# albedoVeg = 0.2



# real time of first time step, duration of time step
# IMPORTANT NOTE: THIS IS NOW UTC TIME ALMOST CERTAINLY AT LEAST FOR SHADING
# print("# IMPORTANT NOTE: THIS IS NOW UTC TIME ALMOST CERTAINLY AT LEAST FOR SHADING")
startTimeYearValue = 2005
startTimeMonthValue = 7
startTimeDayValue = 1
timeStepDurationHoursFloatingPointValue = 1.0  # only tested for one hour!!!!


# lat long for shading (solar radiation)
latitudeOfCatchment = 52.12833333
longitudeOfCatchment = 5.19861111
timeZone = "Europe/Madrid"