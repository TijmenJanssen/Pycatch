# general
import datetime
import glob
import sys
import os
import pathlib
import shutil
import numpy as np
from collections import deque
import time
import threading
import queue
import math
import matplotlib.pyplot as plt

# add path to modules required
sys.path.append("./pcrasterModules/")

import EWS_main_configuration as cfg

# PCRaster itself
import pcraster as pcr
import pcraster.framework as pcrfw


# from pcrasterModules
import datetimePCRasterPython
import interceptionuptomaxstore
import surfacestore
import infiltrationgreenandampt
import subsurfacewateronelayer
import evapotranspirationpenman
import runoffaccuthreshold
import shading
import generalfunctions
import randomparameters
import rainfalleventsfromgammadistribution

# from this folder
import exchangevariables

# TIJMEN CONSTRUCTION 
# Toegevoegd
# - str(path.Path(inputFolder,"x" +"step"))
# - For loop bij einde
# -albedo equation
# -biomass plant height equation
# -probalistic precipitation (checken)

#Toevoegen
#-stomatal conductance
#meteo forcing (ERA-5)
#nakijken discharge totaal voor EWS
# Soil moisture wegschrijven zoals discharge
# Regen prob vergroten
# soilmoisture erbij
# cfg parameters erbij?


#OPMERKINGEN
# aguila --timesteps=[1,9,1] filename
#regolithThicknessHomogeneous heet nog zo, maar is wel gewoon uit het w-model


inputFolder = "inputs"


# only for advanced users
# uncomment following line and comment second line in case of particle filtering
# first time users should use class CatchmentModel(DynamicModel,MonteCarloModel):
# in case of particle filtering, CHANGE ALSO TIMESERIES FILE FOR SCENARIOS!!!!!!!!!
# class CatchmentModel(pcrfw.DynamicModel, pcrfw.MonteCarloModel, pcrfw.ParticleFilterModel):


class CatchmentModel(pcrfw.DynamicModel, pcrfw.MonteCarloModel):
  def __init__(self):
    pcrfw.DynamicModel.__init__(self)
    pcrfw.MonteCarloModel.__init__(self)
    pcr.setclone(cfg.cloneString)
    if cfg.filtering:
      pcrfw.ParticleFilterModel.__init__(self)

  def premcloop(self):
    self.clone = pcr.boolean(cfg.cloneString)
    variable_name  = generalfunctions.file_name_str("demM", '00', step)    
    self.dem = pcr.scalar(str(pathlib.Path(inputFolder,variable_name[0])))
    self.createInstancesPremcloop()
    
    self.durationHistory = cfg.number_of_timesteps_hourly

    # required for reporting as numpy
    self.locations = pcr.cover(pcr.nominal(cfg.locations), 0)
    pcr.report(self.locations, 'locats')
    
    # TIJMEN NU EVEN UIT
    # self.forestNoForest = pcr.boolean(cfg.forestNoForest)
    # idMap = pcr.uniqueid(self.clone)
    # oneLocationPerArea = pcr.areamaximum(idMap, self.forestNoForest) == idMap
    # self.locationsForParameters = pcr.cover(pcr.nominal(pcr.scalar(pcr.ifthen(oneLocationPerArea, self.forestNoForest)) + 1), 0)
    # end required for reporting as numpy

  def initial(self):
      self.timeStepDuration = cfg.timeStepDurationHoursFloatingPointValue
      self.initializeTime(cfg.startTimeYearValue, cfg.startTimeMonthValue, cfg.startTimeDayValue, self.timeStepDuration)
      self.createInstancesInitial()
      self.d_exchangevariables.upwardSeepageFlux = pcr.scalar(0)
      self.d_exchangevariables.evapFromSoilMultiplier = pcr.scalar(1)
      
      #functions ans settings for saving timeseries
      self.historyOfTotQ = deque([])
      self.historyOfSoilMoistureFraction = deque([])
      
      # budgets
      self.d_exchangevariables.cumulativePrecipitation = pcr.scalar(0)      
      
 
      
      # RAINFALL
      
      hours_per_week = 7 * 24
      nr_of_weeks = math.ceil(cfg.number_of_timesteps_hourly / hours_per_week)
              
    # 1 if this week contains rain, else 0 (uniform distribution)
      rairray = np.random.rand(nr_of_weeks)
      rairray[rairray >= (1-cfg.probabilityOfRainstorm)] = 1
      rairray[rairray < (1-cfg.probabilityOfRainstorm)] = 0
      #print(rairray)
    
    # 2D array of size nr_of_weeks x hours_per_week (= number of hours in the hourly model run)
      rairray_2D = np.zeros((nr_of_weeks, hours_per_week))
    
    # First day of the week equals 1 if it rains in this week
      rairray_2D[:, 0] = rairray
     #print(rairray_2D
    
    # Shuffles the hours in the week such that the rain does (most likely) not fall on the same day every week
      rairray_hours = np.copy(rairray_2D)
      idx = np.random.rand(*rairray_hours.shape).argsort(axis=1)
      rairray_hours = np.take_along_axis(rairray_hours, idx, axis=1)
    
    # Flatten the 2D array to a 1D array with the length of the hourly model run
      rairray_hours_flat = rairray_hours.flatten()
    #print(rairray_hours_flat)
    
    # For the set duration of the rainstorm, a 1 is added after an already given 1
      rairray_timeseries = rairray_hours_flat
      for i in range(1, cfg.rainstormDuration):
          rairray_timeseries += np.roll(rairray_hours_flat, i)
    
    #print(rairray_timeseries)
    
    # Boolean timeseries for rainstorm
      self.isRainingTest = rairray_timeseries.astype(bool)
    #print(rairray_timeseries_bool)

    


  def dynamic(self):
    
    # reports
    self.reportComponentsDynamic()
    save_maps = (self.currentTimeStep() % cfg.interval_map_snapshots) == 0 and cfg.map_data is True
    save_mean_timeseries = self.currentTimeStep() == cfg.number_of_timesteps_hourly and cfg.mean_timeseries_data is True
      
    
    # time
    self.d_dateTimePCRasterPython.update()
    timeDatetimeFormat = self.d_dateTimePCRasterPython.getTimeDatetimeFormat()

    # precipitation
    # for calibration
#TIJMEN

    # # OR: rainfallFluxDeterm = pcr.timeinputscalar(cfg.rainfallFluxDetermTimeSeries, pcr.nominal(cfg.rainfallFluxDetermTimeSeriesAreas))
    # rainfallFluxDeterm = pcr.timeinputscalar(cfg.rainfallFluxDetermTimeSeries, pcr.nominal(cfg.cloneString))
   
    # # for the experiments
    # rainfallFlux = rainfallFluxDeterm #generalfunctions.mapNormalRelativeError(rainfallFluxDeterm,0.25)
    # self.d_exchangevariables.cumulativePrecipitation = \
    #         self.d_exchangevariables.cumulativePrecipitation + rainfallFlux * self.timeStepDuration
    
    
    if self.isRainingTest[self.currentTimeStep()]:
        rainfallFlux = pcr.scalar(cfg.expectedRainfallIntensity * generalfunctions.mapgamma(cfg.gammaShapeParameter))
        
    else:
        rainfallFlux = pcr.scalar(0)
        
    rainfallAmount = rainfallFlux       
            
    # # precipitation
    # isRaining,rainfallFlux,rainfallAmount= self.d_rainfalleventsfromgammadistribution.getRainstorm()    
     

    # self.d_exchangevariables.cumulativePrecipitation= \
    #         self.d_exchangevariables.cumulativePrecipitation+rainfallFlux*self.timeStepDuration
            
    
    self.report(rainfallFlux, 'rF')
    # self.report(rainfallAmount, 'rA')
            
    # #####################################
    
    
    

    # interception store
    actualAdditionFluxToInterceptionStore = self.d_interceptionuptomaxstore.addWater(rainfallFlux)
    throughfallFlux = rainfallFlux - actualAdditionFluxToInterceptionStore
    
        
    # surface store
    totalToSurfaceFlux = throughfallFlux + self.d_exchangevariables.upwardSeepageFlux
    potentialToSurfaceStoreFlux = self.d_surfaceStore.potentialToFlux()

    # potential infiltration
    potentialHortonianInfiltrationFlux = self.d_infiltrationgreenandampt.potentialInfiltrationFluxFunction()
    maximumSaturatedOverlandFlowInfiltrationFlux = self.d_subsurfaceWaterOneLayer.getMaximumAdditionFlux()
    potentialInfiltrationFlux = pcr.min(potentialHortonianInfiltrationFlux, maximumSaturatedOverlandFlowInfiltrationFlux)

    # SOIL MOISTURE
    self.d_subsurfaceWaterOneLayer.calculateSoilMoistureFraction()
    variable = self.d_subsurfaceWaterOneLayer.soilMoistureFraction
    variable_mean = np.nanmean(pcr.pcr2numpy(variable, np.NaN))

    self.historyOfSoilMoistureFraction = generalfunctions.keepHistoryOfMaps(
        self.historyOfSoilMoistureFraction, variable_mean, self.durationHistory)

    if save_maps:
        generalfunctions.report_as_map(variable, 'moiM', self.currentSampleNumber(), self.currentTimeStep())
    if save_mean_timeseries:
        variable_mean_array = np.array(self.historyOfSoilMoistureFraction)
        generalfunctions.report_as_array(variable_mean_array, 'moiA', self.currentSampleNumber(),
                                                self.currentTimeStep())

    # abstraction from surface water
    potentialAbstractionFromSurfaceWaterFlux = potentialToSurfaceStoreFlux + potentialInfiltrationFlux
    actualAbstractionFromSurfaceWaterFlux, runoffCubicMetresPerHour = self.d_runoffAccuthreshold.update(
                                          totalToSurfaceFlux, potentialAbstractionFromSurfaceWaterFlux)
    potentialOutSurfaceStoreFlux = self.d_surfaceStore.potentialOutFlux()  
     
    
    self.report(runoffCubicMetresPerHour,'Ro')
    
    # TIJMEN added from w-model to save discharge as timeseries
    downstreamEdge = generalfunctions.edge(self.clone, 2, 0)
    pits = pcr.pcrne(pcr.pit(self.d_runoffAccuthreshold.ldd), 0)
    outflowPoints = pcr.pcrand(downstreamEdge, pits)
    totQ = pcr.ifthen(self.clone,
                 pcr.maptotal(
                      pcr.ifthenelse(outflowPoints, self.d_runoffAccuthreshold.RunoffCubicMetrePerHour, pcr.scalar(0))))

    variable = totQ
    variable_mean = np.nanmean(pcr.pcr2numpy(variable, np.NaN))

    self.historyOfTotQ = generalfunctions.keepHistoryOfMaps(self.historyOfTotQ, variable_mean,
                                                                   self.durationHistory)
    
    self.report(variable ,'test')
    
    if save_mean_timeseries:
        variable_mean_array = np.array(self.historyOfTotQ)
        generalfunctions.report_as_array(variable_mean_array, 'Rq', self.currentSampleNumber(),
                                                self.currentTimeStep())
    

    # infiltration
    availableForInfiltrationFlux = potentialOutSurfaceStoreFlux + actualAbstractionFromSurfaceWaterFlux
    availableForInfiltrationNotExceedingMaximumSaturatedOverlandFlowFlux = pcr.min(
                                          availableForInfiltrationFlux, maximumSaturatedOverlandFlowInfiltrationFlux)
    actualInfiltrationFlux = self.d_infiltrationgreenandampt.update(
                                          availableForInfiltrationNotExceedingMaximumSaturatedOverlandFlowFlux)

    # surface store
    surfaceStoreChange = actualAbstractionFromSurfaceWaterFlux - actualInfiltrationFlux
    self.d_surfaceStore.update(surfaceStoreChange)
    actualAdditionFlux = self.d_subsurfaceWaterOneLayer.addWater(actualInfiltrationFlux)

    if cfg.with_shading:
      # solar radiation (POTRAD, shading effect and inclination)
      fractionReceived, fractionReceivedFlatSurface, shaded = \
                                            self.d_shading.update(timeDatetimeFormat)

      # we assume all cells receive the same solar radiation as measured by the device
      # except for shading, if shading, there is nothing received
      fractionReceived = pcr.ifthenelse(shaded, pcr.scalar(0.0), pcr.scalar(1.0))
    else:
      fractionReceived = pcr.scalar(cfg.fractionReceivedValue)
      fractionReceivedFlatSurface = pcr.scalar(cfg.fractionReceivedFlatSurfaceValue)

    fWaterPotential = self.d_subsurfaceWaterOneLayer.getFWaterPotential()

    # potential evapotranspiration
    airTemperatureDeterm = pcr.timeinputscalar(cfg.airTemperatureDetermString, self.clone)
    airTemperature = airTemperatureDeterm #airTemperatureDeterm+mapnormal()

    relativeHumidityDeterm = pcr.timeinputscalar(cfg.relativeHumidityDetermString, self.clone)
    relativeHumidity = relativeHumidityDeterm #pcr.max(pcr.min(relativeHumidityDeterm+mapnormal()*0.1,pcr.scalar(1.0)),pcr.scalar(0))

    incomingShortwaveRadiationFlatSurface = pcr.timeinputscalar(cfg.incomingShortwaveRadiationFlatSurfaceString, self.clone)
    # incomingShortwaveRadiationFlatSurface = pcr.max(pcr.scalar(0),
    #                              generalfunctions.mapNormalRelativeError(incomingShortwaveRadiationFlatSurfaceDeterm,0.25))

    incomingShortwaveRadiationAtSurface = incomingShortwaveRadiationFlatSurface * fractionReceived

    windVelocityDeterm = pcr.timeinputscalar(cfg.windVelocityDetermString, self.clone)
    windVelocity = windVelocityDeterm #generalfunctions.mapNormalRelativeError(windVelocityDeterm,0.25)

    elevationAboveSeaLevelOfMeteoStation = cfg.elevationAboveSeaLevelOfMeteoStationValue

    potentialEvapotranspirationFlux, \
           potentialEvapotranspirationAmount, \
           potentialEvapotranspirationFromCanopyFlux, \
           potentialEvapotranspirationFromCanopyAmount = \
                            self.d_evapotranspirationPenman.potentialEvapotranspiration(
                            airTemperature,
                            relativeHumidity,
                            incomingShortwaveRadiationAtSurface,
                            incomingShortwaveRadiationFlatSurface,
                            fractionReceivedFlatSurface,
                            windVelocity,
                            elevationAboveSeaLevelOfMeteoStation,
                            fWaterPotential,
                            rainfallFlux < 0.000000000001)

    potentialEvapotranspirationFluxNoNegativeValues = pcr.max(0.0, potentialEvapotranspirationFlux)
    potentialEvapotranspirationFluxFromCanopyNoNegativeValues = pcr.max(0.0, potentialEvapotranspirationFromCanopyFlux)

    # evapotranspirate first from interception store
    actualAbstractionFluxFromInterceptionStore = self.d_interceptionuptomaxstore.abstractWater(
                            potentialEvapotranspirationFluxFromCanopyNoNegativeValues)

    # fraction of soil evapotranspiration depends on evapo from canopy
    evapFromSoilMultiplierMV = (potentialEvapotranspirationFluxFromCanopyNoNegativeValues -
                            actualAbstractionFluxFromInterceptionStore) / \
                            potentialEvapotranspirationFluxFromCanopyNoNegativeValues
    self.d_exchangevariables.evapFromSoilMultiplier = \
                           pcr.ifthenelse(potentialEvapotranspirationFluxNoNegativeValues < 0.0000000000001,
                           pcr.scalar(1), evapFromSoilMultiplierMV)

    # evapotranspirate from subsurface store
    # potentialEvapotranspirationFluxFromSubsurface= \
    #                       pcr.max(0.0,potentialEvapotranspirationFluxNoNegativeValues-actualAbstractionFluxFromInterceptionStore)
    potentialEvapotranspirationFluxFromSubsurface = self.d_exchangevariables.evapFromSoilMultiplier * \
                                                  potentialEvapotranspirationFluxNoNegativeValues
    actualAbstractionFluxFromSubsurface = self.d_subsurfaceWaterOneLayer.abstractWater(potentialEvapotranspirationFluxFromSubsurface)

    # upward seepage from subsurfacestore
    self.d_exchangevariables.upwardSeepageFlux = self.d_subsurfaceWaterOneLayer.lateralFlow()

    # reports
    self.reportComponentsDynamic()
    self.reportRandomParametersDynamic()
    self.printComponentsDynamic()
    if cfg.doReportComponentsDynamicAsNumpy:
      self.reportComponentsDynamicAsNumpy()

    #self.checkBudgets(self.currentSampleNumber(), self.currentTimeStep())

  def postmcloop(self):
    # required for reporting as numpy
    import generalfunctions
    self.timeStepDuration = cfg.timeStepDurationHoursFloatingPointValue # needed in case of forking, else the instances have been deleted
    self.initializeTime(cfg.startTimeYearValue, cfg.startTimeMonthValue, cfg.startTimeDayValue, self.timeStepDuration) # needed in case of forking, else the instances have been deleted
    self.createInstancesInitial() # needed in case of forking, else the instances have been deleted

    if cfg.doReportComponentsDynamicAsNumpy:
      self.reportAsNumpyComponentsPostmcloop()

  def createInstancesPremcloop(self):
    if cfg.with_shading:
       self.d_shading = shading.Shading(self.dem, cfg.latitudeOfCatchment, cfg.longitudeOfCatchment, cfg.timeZone, 1, cfg.timeStepsToReportRqs, cfg.shading_report_rasters)
      # print('no optimization of shading')


  def createInstancesInitial(self):
    import generalfunctions


    if cfg.readDistributionOfParametersFromDisk:
       path = '/home/derek/tmp/'
       maximumInterceptionCapacityPerLAI = pcr.scalar(path + pcrfw.generateNameS('RPic', self.currentSampleNumber()) + '.map')
       ksat = pcr.scalar(path + pcrfw.generateNameS('RPks', self.currentSampleNumber()) + '.map')
       regolithThicknessHomogeneous = pcr.scalar(path + pcrfw.generateNameS('RPrt', self.currentSampleNumber()) + '.map')
       saturatedConductivityMetrePerDay = pcr.scalar(path + pcrfw.generateNameS('RPsc', self.currentSampleNumber()) + '.map')
       multiplierMaxStomatalConductance = pcr.scalar(path + pcrfw.generateNameS('RPmm', self.currentSampleNumber()) + '.map')
    else:
      variable_name  = generalfunctions.file_name_str("micM", "00", step)
      maximumInterceptionCapacityPerLAI = generalfunctions.areauniformBounds(
                                  0.0001, 0.0005, pcr.nominal(1), pcr.scalar(str(pathlib.Path(inputFolder,variable_name[0]))), cfg.createRealizations)
      
      variable_name  = generalfunctions.file_name_str("Iks", '00', step)
      ksat = generalfunctions.areauniformBounds(
                                  0.025, 0.05, pcr.nominal(1), pcr.scalar(str(pathlib.Path(inputFolder,variable_name[0]))), cfg.createRealizations)
      
      variable_name  = generalfunctions.file_name_str("regM", '00', step)
      regolithThicknessHomogeneous = generalfunctions.areauniformBounds(
                                  1.0, 3.5, cfg.areas, pcr.scalar(str(pathlib.Path(inputFolder,variable_name[0]))), cfg.createRealizations)
      
      saturatedConductivityMetrePerDay = generalfunctions.mapuniformBounds(
                                  25.0, 40.0, pcr.scalar(cfg.saturatedConductivityMetrePerDayValue), cfg.createRealizations)
     
      multiplierMaxStomatalConductance = generalfunctions.mapuniformBounds(
                                  0.8, 1.1, pcr.scalar(cfg.multiplierMaxStomatalConductanceValue), cfg.createRealizations)

    if cfg.swapCatchments:
      regolithThicknessHomogeneous = generalfunctions.swapValuesOfTwoRegions(cfg.areas, regolithThicknessHomogeneous, True)

    self.d_randomparameters = randomparameters.RandomParameters(
                    cfg.timeStepsToReportRqs,
                    cfg.randomparameters_report_rasters,
                    maximumInterceptionCapacityPerLAI,
                    ksat,
                    regolithThicknessHomogeneous,
                    saturatedConductivityMetrePerDay,
                    multiplierMaxStomatalConductance)


    # class for exchange variables in initial and dynamic
    # introduced to make filtering possible
    self.d_exchangevariables = exchangevariables.ExchangeVariables(
                                    cfg.timesteps_to_report_some_hourly,
                                    cfg.exchange_report_rasters
                                    )


   



    ################
    # interception #
    ################

#TIJMEN 
    #o: self.ldd = cfg.lddMap
    
    # self.ldd = pcr.ldd(self.dem)
    #self.ldd = pcr.lddcreate(self.dem,10,1e31,1e31,1e31) 
    # of self.ldd = pcr.lddcreate(self.dem,0,0,0,0) 
    
    
    self.ldd = pcr.lddcreate(self.dem,10,1e31,1e31,1e31)       
    self.report(self.ldd, 'ldd')



    initialInterceptionStore = pcr.scalar(0.000001)
    variable_name  = generalfunctions.file_name_str("laiM","00", step)  
    leafAreaIndex = pcr.scalar(str(pathlib.Path(inputFolder,variable_name[0])))
    
    
    

    if cfg.swapCatchments:
      leafAreaIndex = generalfunctions.swapValuesOfTwoRegions(cfg.areas, leafAreaIndex, True)
    gapFraction = pcr.exp(-0.5 * leafAreaIndex)            # equation 40 in Brolsma et al 2010a
    maximumInterceptionStore = maximumInterceptionCapacityPerLAI * leafAreaIndex

    self.d_interceptionuptomaxstore = interceptionuptomaxstore.InterceptionUpToMaxStore(
                                    self.ldd,
                                    initialInterceptionStore,
                                    maximumInterceptionStore,
                                    gapFraction,
                                    cfg.calculateUpstreamTotals,
                                    self.timeStepDurationHours,
                                    cfg.timesteps_to_report_some_hourly,
                                    cfg.interception_report_rasters)

    #################
    # surface store #
    #################

    initialSurfaceStore = pcr.scalar(0.0)
    maxSurfaceStore = pcr.scalar(cfg.maxSurfaceStoreValue)
    self.d_surfaceStore = surfacestore.SurfaceStore(
                        initialSurfaceStore,
                        maxSurfaceStore,
                        self.timeStepDurationHours,
                        cfg.timesteps_to_report_some_hourly,
                        cfg.surfacestore_report_rasters)

    ################
    # infiltration #
    ################

    # N initialMoistureContentFraction taken from 1st July

    # DK
    # we do not use rts and Gs as input to calculate initial moisture fraction to avoid
    # problems when the initial regolith thickness is calibrated (it might be thinner than
    # initialMoistureThick -> problems!)
    # instead, we use initial moisture content fraction as input, read from disk, it is just calculated
    # by pcrcalc 'mergeInitialMoistureContentFraction=Gs000008.761/rts00008.761'
    # note that I also changed the name for the initial soil moisture as a fraction
    
    variable_name  = generalfunctions.file_name_str("moiM", '00', step)    
    initialSoilMoistureFractionFromDisk = pcr.scalar(str(pathlib.Path(inputFolder,variable_name[0])))
      
    if cfg.swapCatchments:
      initialSoilMoistureFractionFromDisk = generalfunctions.swapValuesOfTwoRegions(cfg.areas, initialSoilMoistureFractionFromDisk, True)


    # initial soil moisture as a fraction should not be above soil porosity as a fraction, just a check
    soilPorosityFraction = pcr.scalar(cfg.soilPorosityFractionValue)
    if cfg.swapCatchments:
      soilPorosityFraction = generalfunctions.swapValuesOfTwoRegions(cfg.areas, soilPorosityFraction, True)
    initialSoilMoistureFraction = pcr.min(soilPorosityFraction, initialSoilMoistureFractionFromDisk)
    hf = pcr.scalar(-0.0000001)
    self.d_infiltrationgreenandampt = infiltrationgreenandampt.InfiltrationGreenAndAmpt(
                                    soilPorosityFraction,
                                    initialSoilMoistureFraction,
                                    ksat,
                                    hf,
                                    self.timeStepDurationHours,
                                    cfg.timesteps_to_report_some_hourly,
                                    cfg.infiltration_report_rasters)

    ####################
    # subsurface water #
    ####################

    demOfBedrockTopography = self.dem

   # stream = pcr.boolean(cfg.streamValue)
    theSlope = pcr.slope(self.dem)
    regolithThickness = regolithThicknessHomogeneous

    self.multiplierWiltingPoint = pcr.scalar(1.0)
    limitingPointFraction = pcr.scalar(cfg.limitingPointFractionValue)

    if cfg.swapCatchments:
      limitingPointFraction = generalfunctions.swapValuesOfTwoRegions(cfg.areas, limitingPointFraction, True)
    mergeWiltingPointFractionFS = pcr.scalar(cfg.mergeWiltingPointFractionFSValue)
    if cfg.swapCatchments:
      mergeWiltingPointFractionFS = generalfunctions.swapValuesOfTwoRegions(cfg.areas, mergeWiltingPointFractionFS, True)
    wiltingPointFractionNotChecked = mergeWiltingPointFractionFS * self.multiplierWiltingPoint
    wiltingPointFraction = pcr.min(wiltingPointFractionNotChecked, limitingPointFraction)

    fieldCapacityFraction = pcr.scalar(cfg.fieldCapacityFractionValue)
    if cfg.swapCatchments:
      fieldCapacityFraction = generalfunctions.swapValuesOfTwoRegions(cfg.areas, fieldCapacityFraction, True)

    self.d_subsurfaceWaterOneLayer = subsurfacewateronelayer.SubsurfaceWaterOneLayer(
                                   self.ldd,
                                   demOfBedrockTopography,
                                   regolithThickness,
                                   initialSoilMoistureFraction,
                                   soilPorosityFraction,
                                   wiltingPointFraction,
                                   fieldCapacityFraction,
                                   limitingPointFraction,
                                   saturatedConductivityMetrePerDay,
                                   cfg.calculateUpstreamTotals,
                                   self.timeStepDurationHours,
                                   cfg.timesteps_to_report_some_hourly,
                                   cfg.subsurface_report_rasters)


    ##########
    # runoff #
    ##########

    self.d_runoffAccuthreshold = runoffaccuthreshold.RunoffAccuthreshold(
                               self.ldd,
                               self.timeStepDurationHours,
                               cfg.timeStepsToReportRqs,
                               cfg.runoff_report_rasters)

    ######################
    # evapotranspiration #
    ######################

    albedo = gapFraction * cfg.albedoSoil + (1-gapFraction) * cfg.albedoVeg    
    if cfg.swapCatchments:
      albedo = generalfunctions.swapValuesOfTwoRegions(cfg.areas, albedo, True)
      
#TIJMEN
    # O:     maxStomatalConductance = pcr.scalar(cfg.maxStomatalConductanceValue) * multiplierMaxStomatalConductance
    # Nu is het dummie getal, nakijken in literatuur
    
    maxStomatalConductance = pcr.scalar(0.004) * multiplierMaxStomatalConductance
    if cfg.swapCatchments:
      maxStomatalConductance = generalfunctions.swapValuesOfTwoRegions(cfg.areas, maxStomatalConductance, True)

    variable_name  = generalfunctions.file_name_str("bioM", "00", step)
    biomass = pcr.scalar(str(pathlib.Path(inputFolder,variable_name[0])))

    
    vegetationHeight = (biomass/0.679)**(1/0.997) #Proulx 2021         
    
    if cfg.swapCatchments:
      vegetationHeight = generalfunctions.swapValuesOfTwoRegions(cfg.areas, vegetationHeight, True)
    self.d_evapotranspirationPenman = evapotranspirationpenman.EvapotranspirationPenman(
                                         self.timeStepDurationHours,
                                         albedo,
                                         maxStomatalConductance,
                                         vegetationHeight,
                                         leafAreaIndex,
                                         cfg.timesteps_to_report_some_hourly,
                                         cfg.evapotrans_report_rasters)


  def reportComponentsDynamic(self):
    """report dynamic components as PCRaster maps
    components, the modules that are reported
    see also reportAsNumpyComponentsPostmcloop
    """
    components = [ \
                  self.d_exchangevariables, \
                 self.d_randomparameters, \
                 self.d_interceptionuptomaxstore, \
                 # self.d_surfaceStore, \
                 self.d_infiltrationgreenandampt, \
                 self.d_evapotranspirationPenman, \
                 self.d_runoffAccuthreshold, \
                 self.d_subsurfaceWaterOneLayer
                 ]

    if cfg.with_shading:
      components.append(self.d_shading)

    for component in components:
      component.reportAsMaps(self.currentSampleNumber(), self.currentTimeStep())


  def reportComponentsDynamicAsNumpy(self):
    """report dynamic components as PCRaster maps
    components, the modules that are reported
    see also reportAsNumpyComponentsPostmcloop
    """

    # report dynamic components as numpy, see also 'reportAsNupyComponentsPostmcloop'
    self.d_runoffAccuthreshold.reportAsNumpy(self.locations, self.currentSampleNumber(), self.currentTimeStep())
    self.d_subsurfaceWaterOneLayer.reportAsNumpy(self.locations, self.currentSampleNumber(), self.currentTimeStep())
    self.d_interceptionuptomaxstore.reportAsNumpy(self.locations, self.currentSampleNumber(), self.currentTimeStep())
    self.d_randomparameters.reportAsNumpy(self.locationsForParameters, self.currentSampleNumber(), self.currentTimeStep())


  def reportAsNumpyComponentsPostmcloop(self):
    """report dynamic components as PCRaster maps
    componentsToReportAsNumpy should correspond with the numpy one in reportComponentsDynamic
    """
    componentsToReportAsNumpy = [
                 self.d_runoffAccuthreshold,
                 self.d_subsurfaceWaterOneLayer,
                 self.d_interceptionuptomaxstore,
                 self.d_randomparameters
                                ]
    for component in componentsToReportAsNumpy:
      component.reportAsNumpyPostmcloop(range(1, cfg.nrOfSamples + 1), range(1, cfg.number_of_timesteps_hourly + 1))


  def reportRandomParametersDynamic(self):
    self.d_randomparameters.reportAtLastTimeStep(self.currentSampleNumber(), self.currentTimeStep(), self.nrTimeSteps())

  def printMemberVariables(self):
    import generalfunctions
    components = [
                 self.d_exchangevariables,
                 self.d_interceptionuptomaxstore,
                 self.d_surfaceStore,
                 self.d_infiltrationgreenandampt,
                 self.d_evapotranspirationPenman,
                 self.d_runoffAccuthreshold,
                 self.d_shading,
                 self.d_subsurfaceWaterOneLayer
                 ]

    for component in components:
      generalfunctions.printMemberVariables(component)

  def printComponentsDynamic(self):
    self.d_dateTimePCRasterPython.printit()

  def initializeTime(self, startTimeYear, startTimeMonth, startTimeDay, timeStepDurationHours):
    startTime = datetime.datetime(year=startTimeYear, month=startTimeMonth, day=startTimeDay)
    self.timeStepDurationHours = timeStepDurationHours
    self.timeStepDatetimeFormat = datetime.timedelta(hours=self.timeStepDurationHours)
    self.d_dateTimePCRasterPython = datetimePCRasterPython.DatetimePCRasterPython \
                                  (startTime, self.timeStepDatetimeFormat)

  def checkBudgets(self, currentSampleNumber, currentTimeStep):
    # DK not sure this is still correct

    increaseInPrecipitationStore = 0.0 - self.d_exchangevariables.cumulativePrecipitation
    pcr.report(increaseInPrecipitationStore, pcrfw.generateNameST('incP', currentSampleNumber, currentTimeStep))

    increaseInInterceptionStore = self.d_interceptionuptomaxstore.budgetCheck(currentSampleNumber, currentTimeStep)
    pcr.report(increaseInInterceptionStore, pcrfw.generateNameST('incI', currentSampleNumber, currentTimeStep))

    increaseInSurfaceStore = self.d_surfaceStore.budgetCheck(currentSampleNumber, currentTimeStep)
    pcr.report(increaseInSurfaceStore, pcrfw.generateNameST('incS', currentSampleNumber, currentTimeStep))
    increaseInSurfaceStoreQM = pcr.catchmenttotal(increaseInSurfaceStore, self.ldd) * pcr.cellarea()
    pcr.report(increaseInSurfaceStoreQM, pcrfw.generateNameST('testb', currentSampleNumber, currentTimeStep))

    # let op: infiltration store is directly passed to subsurface store, thus is not a real store
    increaseInInfiltrationStore = self.d_infiltrationgreenandampt.budgetCheck(currentSampleNumber, currentTimeStep)

    increaseInSubSurfaceWaterStore, lateralFlowInSubsurfaceStore, abstractionFromSubSurfaceWaterStore = \
                             self.d_subsurfaceWaterOneLayer.budgetCheck(currentSampleNumber, currentTimeStep)
    increaseInSubSurfaceStoreQM = pcr.catchmenttotal(increaseInSubSurfaceWaterStore, self.ldd) * pcr.cellarea()

    increaseInRunoffStoreCubicMetresInUpstreamArea = self.d_runoffAccuthreshold.budgetCheck()

    totalIncreaseInStoresCubicMetresInUpstreamArea = 0.0
    stores = [increaseInPrecipitationStore, increaseInInterceptionStore, increaseInSurfaceStore, increaseInSubSurfaceWaterStore]
    for store in stores:
      increaseInStoreCubicMetresInUpstreamArea = pcr.catchmenttotal(store, self.ldd) * pcr.cellarea()
      totalIncreaseInStoresCubicMetresInUpstreamArea = totalIncreaseInStoresCubicMetresInUpstreamArea + \
                                                     increaseInStoreCubicMetresInUpstreamArea

    pcr.report(totalIncreaseInStoresCubicMetresInUpstreamArea, pcrfw.generateNameST('inSt', currentSampleNumber, currentTimeStep))
    pcr.report(increaseInRunoffStoreCubicMetresInUpstreamArea, pcrfw.generateNameST('inRu', currentSampleNumber, currentTimeStep))
    pcr.report(pcr.catchmenttotal(self.d_exchangevariables.upwardSeepageFlux, self.ldd) * pcr.cellarea(), pcrfw.generateNameST('inSe', currentSampleNumber, currentTimeStep))
    # total budget is total increase in stores plus the upward seepage flux for each ts that is passed to the next
    # timestep and thus not taken into account in the current timestep budgets
    budget = totalIncreaseInStoresCubicMetresInUpstreamArea + increaseInRunoffStoreCubicMetresInUpstreamArea + \
           lateralFlowInSubsurfaceStore * pcr.cellarea() + pcr.catchmenttotal(abstractionFromSubSurfaceWaterStore, self.ldd) * pcr.cellarea() + \
           pcr.catchmenttotal(self.d_exchangevariables.upwardSeepageFlux, self.ldd) * pcr.cellarea()
    pcr.report(budget, pcrfw.generateNameST('B-tot', currentSampleNumber, currentTimeStep))
    budgetRel = budget / increaseInRunoffStoreCubicMetresInUpstreamArea
    pcr.report(budgetRel, pcrfw.generateNameST('B-rel', currentSampleNumber, currentTimeStep))

  def suspend(self):
    import generalfunctions
    if self.currentTimeStep() != cfg.number_of_timesteps_hourly:
      self.timeStepForResume = self.currentTimeStep()

      components = [ self.d_exchangevariables,
                   self.d_randomparameters,
                   self.d_interceptionuptomaxstore,
                   self.d_surfaceStore,
                   self.d_infiltrationgreenandampt,
                   self.d_evapotranspirationPenman,
                   self.d_runoffAccuthreshold, \
                   # self.d_shading, \
                   self.d_subsurfaceWaterOneLayer]

      for component in components:
        generalfunctions.reportMemberVariablesOfAClassForSuspend(component, self.currentTimeStep(), self.currentSampleNumber())

  def updateWeight(self):
    print('#### UPDATEWEIGHTING')
    print('filter period', self.filterPeriod())
    print('filter timestep ', self._d_filterTimesteps[self.filterPeriod() - 1])
    print('lijst ', self._d_filterTimesteps)
    print('filter sample ', self.currentSampleNumber())
    modelledData = self.readmap('Rqs')
    observations = self.readDeterministic('observations/Rqs')
    #observations=pcr.ifthen(pit(self.ldd) != 0,syntheticData)
    measurementErrorSD = 3.0 * observations + 1.0
    sum = pcr.maptotal(((modelledData - observations)**2) / (2.0 * (measurementErrorSD**2)))
    weight = pcr.exp(-sum)
    weightFloatingPoint, valid = pcr.cellvalue(weight, 1)
    return weightFloatingPoint

  def resume(self):
    print('#### RESUMING')
    print('filter timesteps', self._d_filterTimesteps)
    print('filter period', self.filterPeriod())
    print('filter timestep', self._d_filterTimesteps[self.filterPeriod() - 2])

    import generalfunctions

    # rerun initial
    self.timeStepDuration = cfg.timeStepDurationHoursFloatingPointValue
    self.initializeTime(cfg.startTimeYearValue, cfg.startTimeMonthValue, cfg.startTimeDayValue, self.timeStepDuration)
    self.createInstancesInitial()
    self.d_exchangevariables.upwardSeepageFlux = pcr.scalar(0)
    self.d_exchangevariables.evapFromSoilMultiplier = pcr.scalar(1)
    self.d_exchangevariables.cumulativePrecipitation = pcr.scalar(0)

    # resume time information
    self.d_dateTimePCRasterPython.resume(self._d_filterTimesteps[self.filterPeriod() - 2])

    components = [ self.d_exchangevariables,
                  self.d_randomparameters,
                  self.d_interceptionuptomaxstore,
                  self.d_surfaceStore,
                  self.d_infiltrationgreenandampt,
                  self.d_evapotranspirationPenman,
                  self.d_runoffAccuthreshold, \
                  # self.d_shading, \
                  self.d_subsurfaceWaterOneLayer]

    for component in components:
      generalfunctions.readMemberVariablesOfAClassForResume(
                       component, self._d_filterTimesteps[self.filterPeriod() - 2], self.currentSampleNumber())

    # remove files used to resume
    for filename in glob.glob(str(self.currentSampleNumber()) + '/stateVar/*/*'):
      os.remove(filename)
  

    





def timesteps(start_shift, end_shift, steps_in_shift= cfg.stepsInShift, steps_ba_shift=int((cfg.stepsTotal-cfg.stepsInShift)/2)):
    first_state_duration = (start_shift/cfg.interval_map_snapshots) - 1
    shift_duration = (end_shift/cfg.interval_map_snapshots) - (start_shift/cfg.interval_map_snapshots)
    second_state_duration = (cfg.number_of_timesteps_weekly/cfg.interval_map_snapshots) - (end_shift/cfg.interval_map_snapshots)

    num1, num2, num3 = steps_ba_shift, steps_in_shift, steps_ba_shift
    print(num1, num2, num3)
    if first_state_duration < num1:
        num1 = math.floor(first_state_duration)
        print("Timesteps before transition set to:", num1)
    if shift_duration < num2:
        num2 = math.ceil(shift_duration)
        print("Timesteps during transition set to:", num2)
    if second_state_duration < num3:
        num3 = math.floor(second_state_duration)
        print("Timesteps after transition set to:", num3)

    before_shift = np.linspace(1,
                               math.floor(start_shift/cfg.interval_map_snapshots) - 1,
                               num1, dtype='int') * cfg.interval_map_snapshots
    during_shift = np.linspace(math.floor(start_shift/cfg.interval_map_snapshots),
                               math.ceil(end_shift/cfg.interval_map_snapshots),
                               num2, dtype='int') * cfg.interval_map_snapshots
    after_shift = np.linspace(math.ceil(end_shift/cfg.interval_map_snapshots) + 1,
                               cfg.number_of_timesteps_weekly/cfg.interval_map_snapshots,
                               num3, dtype='int') * cfg.interval_map_snapshots
    timesteps = np.concatenate((before_shift, during_shift, after_shift))
    return timesteps
    
    
    
    # Plot & user input functions
## Threading to enter start & end point while plot is open
my_queue = queue.Queue()
def store_in_queue(f):
    def wrapper(*args):
        my_queue.put(f(*args))
    return wrapper

def plot():
    plt.plot(np.loadtxt('./inputs_from_weekly/1/bioA0005.200.numpy.txt'))
    plt.show()

@store_in_queue
def user_inputs():
    print("The total number of timesteps equals: ", cfg.number_of_timesteps_weekly, "steps.")
    start_shift = input("Enter starting point of shift:")
    end_shift = input("Enter end point of shift:")
    #plt.close()
    return start_shift, end_shift

mthread = threading.Thread(target=user_inputs)
mthread.start()
plot()
shift_data = my_queue.get()
mthread.join()

## end of threading

start_shift = shift_data[0]
end_shift = shift_data[1]
timesteps_ = timesteps(int(start_shift), int(end_shift))
print(timesteps_)

# For timer
start_time = time.time()


# TIJMEN: ORIGINEEL for k, step in enumerate(np.arange(100, cfg.number_of_timesteps_weekly + 100, 100)):
# In np.arragne(x,y,z) is x het eerste snapshot, y de laatste en z de stapgrootte ertussen. 
# Voor specifiek snapshot gebruik #for k, step in enumerate([x]): met voor x eg 100

# TO DO: LAATSTE SNAPSHOT STAAT NOG ONDER MAP 1, OPZICH MAAKT DAT NIET UIT

for k, step in enumerate(timesteps_):
#for k, step in enumerate(np.arange(100, cfg.number_of_timesteps_weekly + 100, cfg.stepsizeSnapshots)):
#for k, step in enumerate([100]):
     snapshot_number = str(k).zfill(2)
     

     if cfg.filtering:     

        myModel = CatchmentModel()
        dynamicModel = pcrfw.DynamicFramework(myModel, cfg.number_of_timesteps_hourly)
        mcModel = pcrfw.MonteCarloFramework(dynamicModel, cfg.nrOfSamples)
        mcModel.setForkSamples(True, 10)
        # pfModel = SequentialImportanceResamplingFramework(mcModel)
        pfModel = pcrfw.ResidualResamplingFramework(mcModel)
        filterTimestepsNoSelection = range(3750, cfg.number_of_timesteps_hourly + 1, 25)
        periodsToExclude = [
            [2617, 2976],
            [3649, 3689],
            [4173, 4416],
            [4046, 4366],
            [5281, 6075]
        ]
        filterTimesteps = generalfunctions.removePeriodsFromAListOfTimesteps(filterTimestepsNoSelection,
                                                                             periodsToExclude)
        pfModel.setFilterTimesteps(filterTimesteps)
        pfModel.run()
        
        dir_name = os.path.join('./h' + snapshot_number + '/')
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)
        os.rename('./1/', dir_name)

     else:                      
        
        myModel = CatchmentModel()
        dynamicModel = pcrfw.DynamicFramework(myModel, cfg.number_of_timesteps_hourly)
        mcModel = pcrfw.MonteCarloFramework(dynamicModel, cfg.nrOfSamples)
        mcModel.setForkSamples(True, 10)
        mcModel.run()
        dynamicModel.run()
        
        dir_name = os.path.join('./h' + snapshot_number + '/')
        if os.path.isdir(dir_name):
            shutil.rmtree(dir_name)
        os.rename('./1/', dir_name)
        
        
end_time = time.time()
print('Elapsed timee for Pychatch h-model loops is', ((end_time - start_time)/60), 'minutes' )