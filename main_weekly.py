# general
import sys
import datetime
from collections import deque
import glob
import math

sys.path.append("../pycatch/pcrasterModules/")

# from PCRaster modules
import generalfunctions
import datetimePCRasterPython
import interceptionuptomaxstore
import surfacestore
import infiltrationonlyksat
import subsurfacewateronelayer
import runoffaccuthreshold
import rainfalleventsfromgammadistribution
import exchangevariables_weekly
import soilwashMMF
import regolith
import bedrockweathering
import evapotranspirationsimple
import biomassmodifiedmay
import baselevel
import creep

import configuration_weekly as cfg

# PCRaster itself
from pcraster import *
from pcraster.framework import *

if cfg.fixedStates:
  cfg.numberOfTimeSteps=52*50
  fixedStatesReg=spatial(scalar(float(sys.argv[1])))
  fixedStatesBio=spatial(scalar(float(sys.argv[2])))

timeStepsWithStatsCalculated = range(cfg.intervalForStatsCalculated, \
                               cfg.numberOfTimeSteps,cfg.intervalForStatsCalculated)

def calculateGapFractionAndMaxIntStoreFromLAI(leafAreaIndex):
  maximumInterceptionCapacityPerLAI=scalar(0.001)
  gapFraction=exp(-0.5*leafAreaIndex)            # equation 40 in Brolsma et al 2010a
  maximumInterceptionStore=maximumInterceptionCapacityPerLAI*leafAreaIndex
  return gapFraction, maximumInterceptionStore

class CatchmentModel(DynamicModel,MonteCarloModel):

  def __init__(self):
    DynamicModel.__init__(self)
    MonteCarloModel.__init__(self)
    setclone('inputs_weekly/clone.map')

    # fix the seed for random functions
    setrandomseed(101)

    if cfg.filtering:
      ParticleFilterModel.__init__(self)


  def premcloop(self):

    self.clone=boolean("inputs_weekly/clone.map")
    self.numberOfCellsOnMap=maptotal(ifthenelse(self.clone,scalar(1),scalar(1)))

    # not essential this code
    # used for fixedStatesLoop.py
    # edgesMap=generalfunctions.edgeZone(self.clone,2.1)
    # self.mlocs=uniqueid(~ edgesMap)
    # self.report(self.mlocs,'testMap')

    # locations where values are reported as a numpy array to disk
    self.mlocs=nominal("inputs_weekly/mlocs")  # multiple report locations, read from 'mlocs.map'
    self.aLocation=self.mlocs

    # zone reported at each report location
    if cfg.fixedStates:
      # option 2 used for fixedStatesLoop.py
      # one big zone, resulting in the same value for each report
      # location
      edgesMap=generalfunctions.edgeZone(self.clone,2.1)
      self.zoneMap=ifthen(~ edgesMap,boolean(1))
      self.report(self.zoneMap,'testMap')
    else:
      # option 1 for normal runs, zones across the hillslope for each
      # location on mlocs
      self.zoneMap=nominal("inputs_weekly/zonsc.map")

    self.createInstancesPremcloop()

    self.durationHistory=207

    # time step duration in hours, typically (and only tested) one week, i.e. 7.0*24.0
    self.timeStepDuration = 7.0 * 24.0


  def initial(self):

    # TIME BEING DIVIDE BY 100 TO AVOID IT RUNS TOO LONG -  DK 210519 not clear what this comment is
    self.initializeTime(2001,2,26,self.timeStepDuration)

    self.createInstancesInitial()

    self.d_exchangevariables.upwardSeepageFlux=scalar(0)
    self.d_exchangevariables.evapFromSoilMultiplier=scalar(1)

    self.timeStepDurationYears=self.timeStepDuration/(365.0*24.0)

    self.actualAbstractionFluxFromSubsurface=0.0
   
    # functions and settings for calculating statistics
    self.historyOfSoilMoistureFraction=deque([])
    self.historyOfBiomass=deque([])
    self.historyOfRegolithThickness=deque([])
    self.historyOfDem=deque([])
    self.historyOfTotQ=deque([])
    nrSampleLocs=100
    fractionShortDistance=0.4
    separationDistance=3
    import generalfunctions  # not sure why this needs to be imported again
    self.samples=generalfunctions.samplingScheme(self.clone, nrSampleLocs, fractionShortDistance, separationDistance,0,0)
    self.report(self.samples,'samples')
    self.someLocs=pcrne(self.samples,0)

    # initial setting for saving grazing pressure 
    self.grazingPressureArray=numpy.empty([0])

    # budgets
    self.d_exchangevariables.cumulativePrecipitation=scalar(0)

    # initial values
    self.grazingRate=0.0
    self.runoffMetreWaterDepthPerHour=scalar(0.0)
    self.creepDeposition=spatial(scalar(0.0))



  def dynamic(self):

    import generalfunctions   # not sure why this needs to be imported again

    #option to print time info
    #print self.currentTimeStep()
    # time
    #self.d_dateTimePCRasterPython.update()
    #timeDatetimeFormat=self.d_dateTimePCRasterPython.getTimeDatetimeFormat()

    ## biomass
    if cfg.fixedStates:
      self.d_regolithdemandbedrock.setNewRegolith(spatial(scalar(fixedStatesReg)))
      self.d_biomassModifiedMay.setNewBiomass(spatial(scalar(fixedStatesBio)))

    
    # grazing pressure driver
    # this below increases grazing pressure and then reduces it again
    grazingRateIncreaseTotal=0.0003
    grazingRateIncrease=grazingRateIncreaseTotal/(cfg.numberOfTimeSteps/2.0)
    if self.currentTimeStep() < (cfg.numberOfTimeSteps/2.0):
      self.grazingRate=self.grazingRate+grazingRateIncrease
    else:
      self.grazingRate=self.grazingRate-grazingRateIncrease

    self.grazingPressureArray=numpy.append(self.grazingPressureArray,self.grazingRate)

    runoffMetreWaterDepthPerWeek=self.runoffMetreWaterDepthPerHour*cfg.theDurationOfRainstorm
    self.biomass,self.LAI=self.d_biomassModifiedMay.update(self.actualAbstractionFluxFromSubsurface, \
                          runoffMetreWaterDepthPerWeek,self.grazingRate)

    # update gap fraction and maximum interception store
    gapFraction,maxIntStore=calculateGapFractionAndMaxIntStoreFromLAI(self.LAI)
    self.d_interceptionuptomaxstore.setGapFraction(gapFraction)
    self.d_interceptionuptomaxstore.setMaximumStore(maxIntStore)

    # update stone cover for erosion
    fractionOfVegetationAboveSoil=0.7
    vegetationCoverForErosion=(1.0-gapFraction)*fractionOfVegetationAboveSoil
    self.d_soilwashMMF.updateStoneOrVegetationCover(vegetationCoverForErosion)
   
    # precipitation
    isRaining,rainfallFlux,rainfallAmount= self.d_rainfalleventsfromgammadistribution.getRainstorm() 

    self.d_exchangevariables.cumulativePrecipitation= \
            self.d_exchangevariables.cumulativePrecipitation+rainfallFlux*self.timeStepDuration
  

    if isRaining:
      # interception store
      actualAdditionFluxToInterceptionStore=self.d_interceptionuptomaxstore.addWater(rainfallFlux)
      throughfallFlux=rainfallFlux-actualAdditionFluxToInterceptionStore
    
      # surface store
      # time being no upward seepage
      #totalToSurfaceFlux=throughfallFlux+self.d_exchangevariables.upwardSeepageFlux
      totalToSurfaceFlux=throughfallFlux
      potentialToSurfaceStoreFlux=self.d_surfaceStore.potentialToFlux()

      # potential infiltration
      self.d_infiltrationonlyksat.setSaturatedConductivityFluxAsFunctionOfBiomass(self.biomass)
      potentialHortonianInfiltrationFlux=self.d_infiltrationonlyksat.potentialInfiltrationFluxFunction()
      maximumSaturatedOverlandFlowInfiltrationFlux=self.d_subsurfaceWaterOneLayer.getMaximumAdditionFlux()
      potentialInfiltrationFlux=min(potentialHortonianInfiltrationFlux,maximumSaturatedOverlandFlowInfiltrationFlux)

      # abstraction from surface water
      potentialAbstractionFromSurfaceWaterFlux=potentialToSurfaceStoreFlux + potentialInfiltrationFlux
      actualAbstractionFromSurfaceWaterFlux,runoffCubicMetrePerHour=self.d_runoffAccuthreshold.update( \
                                            totalToSurfaceFlux,potentialAbstractionFromSurfaceWaterFlux)
      potentialOutSurfaceStoreFlux=self.d_surfaceStore.potentialOutFlux()

      # infiltration
      availableForInfiltrationFlux=potentialOutSurfaceStoreFlux+actualAbstractionFromSurfaceWaterFlux
      availableForInfiltrationNotExceedingMaximumSaturatedOverlandFlowFlux=min( \
                                            availableForInfiltrationFlux,maximumSaturatedOverlandFlowInfiltrationFlux)
      actualInfiltrationFlux=self.d_infiltrationonlyksat.update( \
                                            availableForInfiltrationNotExceedingMaximumSaturatedOverlandFlowFlux)

      # surface store
      surfaceStoreChange=actualAbstractionFromSurfaceWaterFlux-actualInfiltrationFlux
      self.d_surfaceStore.update(surfaceStoreChange)
      actualAdditionFlux=self.d_subsurfaceWaterOneLayer.addWater(actualInfiltrationFlux)
      # empty it again
      self.d_surfaceStore.emptyIt()

      # surface wash
      self.runoffMetreWaterDepthPerHour=runoffCubicMetrePerHour/cellarea()
      netDeposition, netDepositionMetre, lateralFluxKg, totalDetachKgPerCell, transportCapacityKgPerCell= \
                                            self.d_soilwashMMF.calculateWash( \
                                            self.runoffMetreWaterDepthPerHour,rainfallFlux,throughfallFlux)
      # LET OP: dit is metre flux, maar dat zou het zijn als er slechts 1 regenbui is per jaar
      # het is dus de ene week uitgemiddeld over een jaar
      # om echt iets te krijgen met een eenheid m/jaar (dwz wat 'zou' de depositie zijn
      # als dit event elke week zou optreden), moet dit keer 52 weken
      # hetzelfde geldt voor actual deposition flux hieronder
      netDepositionMetreFlux=netDepositionMetre/self.timeStepDurationRegolithInYears
      ##LDD, surface##
      if cfg.changeGeomorphology:
        actualDepositionFlux=self.d_regolithdemandbedrock.updateWithDeposition(netDepositionMetreFlux)
        regolithThickness,demOfBedrock,dem,bedrockLdd,surfaceLdd=self.d_regolithdemandbedrock.getRegolithProperties()
        amountOfMoistureThickNetAdded=self.d_subsurfaceWaterOneLayer.updateRegolithThickness(regolithThickness)
        self.d_soilwashMMF.setSurfaceProperties(surfaceLdd,dem)
        self.d_runoffAccuthreshold.setSurfaceProperties(surfaceLdd)

    
    else:
      # surface wash
      netDeposition, netDepositionMetre, lateralFluxKg, totalDetachKgPerCell, transportCapacityKgPerCell= \
                                            self.d_soilwashMMF.noWash()
      actualDepositionFlux=spatial(scalar(0))
      self.runoffMetreWaterDepthPerHour=scalar(0)

    if cfg.changeGeomorphology:
      # random noise
      netDepositionMetreNoiseFlux=normal(1)/5000
      ##LDD, surface##
      actualDepositionNoiseFlux=self.d_regolithdemandbedrock.updateWithDeposition(netDepositionMetreNoiseFlux)
      regolithThickness,demOfBedrock,dem,bedrockLdd,surfaceLdd=self.d_regolithdemandbedrock.getRegolithProperties()

    # potential evapotranspiration, m/hour
    fWaterPotential=self.d_subsurfaceWaterOneLayer.getFWaterPotential()
    potentialEvapotranspirationFlux=self.d_evapotranspirationSimple.potentialEvapotranspiration(fWaterPotential,self.biomass)

    # evapotranspirate first from interception store
    # assume this does not depend on vegetation, and does not influence transpiration
    # assume it immediately empties (ie, within a week)
    potentialEvaporationFromInterceptionStore=99999.9
    actualAbstractionFluxFromInterceptionStore=self.d_interceptionuptomaxstore.abstractWater( \
                            potentialEvaporationFromInterceptionStore)

    # evapotranspirate from subsurface store
    potentialEvapotranspirationFluxFromSubsurface= \
                      max(0.0,potentialEvapotranspirationFlux)

    self.actualAbstractionFluxFromSubsurface= \
                      self.d_subsurfaceWaterOneLayer.abstractWater(potentialEvapotranspirationFluxFromSubsurface)

    # lateral flow in subsurface and upward seepage from subsurfacestore
    # typically switched off and never tested
    ##self.d_exchangevariables.upwardSeepageFlux=self.d_subsurfaceWaterOneLayer.lateralFlow()

    #self.checkBudgets(self.currentSampleNumber(), self.currentTimeStep())
    #self.printMemberVariables()

    ####
    # geomorpology
    ####
    if cfg.changeGeomorphology and (self.currentTimeStep() % 52 == 0):

      # bedrock weathering
      regolithThickness,demOfBedrock,dem,bedrockLdd,surfaceLdd=self.d_regolithdemandbedrock.getRegolithProperties()
      bedrockWeatheringFlux=self.d_bedrockweathering.weatheringRate(regolithThickness)
      ###LDD, bedrock###
      self.d_regolithdemandbedrock.updateWithBedrockWeathering(bedrockWeatheringFlux)

      # creep
      regolithThickness,demOfBedrock,dem,bedrockLdd,surfaceLdd=self.d_regolithdemandbedrock.getRegolithProperties()
      newRegolithThickness,outflow,flowOverBoundaries,correctedFactor,amountX,amountY,inflowX,inflowY= \
                               self.d_creep.diffuse(regolithThickness,dem,1)
      self.creepDeposition=newRegolithThickness-regolithThickness

      ###LDD, surface###
      self.d_regolithdemandbedrock.setNewRegolith(newRegolithThickness)

      # update bedrock with baselevel change
      baselevel=self.d_baselevel.getBaselevel(self.currentTimeStep())
      ###LDD, surface, bedrock###
      self.d_regolithdemandbedrock.setBaselevel(baselevel)

      # update subsurface store with new regolith thickness
      regolithThickness,demOfBedrock,dem,bedrockLdd,surfaceLdd=self.d_regolithdemandbedrock.getRegolithProperties()
      amountOfMoistureThickNetAdded=self.d_subsurfaceWaterOneLayer.updateRegolithThickness(regolithThickness)
      # no lateral flow, so bedrock does not need to be updated
      # self.d_subsurfaceWaterOneLayer=updateBedrock(self,bedRockLdd,demOfBedrock)

      # update soil wash and runoff with new surface properties
      self.d_soilwashMMF.setSurfaceProperties(surfaceLdd,dem)
      self.d_runoffAccuthreshold.setSurfaceProperties(surfaceLdd)

    # option to do a test run with a very thin regolith
    #print 'test with thin regolith'
    #if self.currentTimeStep() == 6000:
    #  self.d_regolithdemandbedrock.setNewRegolith(spatial(scalar(0.001)))

    # reports
    self.reportComponentsDynamic()
    #self.printComponentsDynamic()


    calculateStats = (self.currentTimeStep()% cfg.intervalForStatsCalculated ) == 0

    # some extra outputs
    if calculateStats:
      # growth part
      meanVariable=areaaverage(self.d_biomassModifiedMay.growthPart,self.zoneMap)
      generalfunctions.reportLocationsAsNumpyArray( \
                       self.aLocation,meanVariable,'gpA',self.currentSampleNumber(),self.currentTimeStep())
      # grazing part
      meanVariable=0.0-areaaverage(spatial(self.d_biomassModifiedMay.grazing),self.zoneMap)
      generalfunctions.reportLocationsAsNumpyArray( \
                       self.aLocation,meanVariable,'grA',self.currentSampleNumber(),self.currentTimeStep())
      # net growth 
      meanVariable=areaaverage(self.d_biomassModifiedMay.netGrowth,self.zoneMap)
      generalfunctions.reportLocationsAsNumpyArray( \
                       self.aLocation,meanVariable,'grNA',self.currentSampleNumber(),self.currentTimeStep())
      # net deposition
      #meanVariable=areaaverage(netDepositionMetre,self.zoneMap)
      meanVariable=areaaverage(actualDepositionFlux,self.zoneMap)
      generalfunctions.reportLocationsAsNumpyArray( \
                       self.aLocation,meanVariable,'depA',self.currentSampleNumber(),self.currentTimeStep())
      # net weathering
      meanVariable=areaaverage(self.d_bedrockweathering.weatheringMetrePerYear,self.zoneMap)
      generalfunctions.reportLocationsAsNumpyArray( \
                       self.aLocation,meanVariable,'weaA',self.currentSampleNumber(),self.currentTimeStep())
      # net creep deposition
      meanVariable=areaaverage(self.creepDeposition,self.zoneMap)
      generalfunctions.reportLocationsAsNumpyArray( \
                       self.aLocation,meanVariable,'creA',self.currentSampleNumber(),self.currentTimeStep())


    if self.currentTimeStep() == cfg.numberOfTimeSteps:
      name=generateNameS('grazing', self.currentSampleNumber()) 
      numpy.save(name,self.grazingPressureArray)
    
    
  def postmcloop(self):
    import generalfunctions
    names=['grA', 'gpA', 'grNA','depA','weaA','creA']
    for name in names:
      aVariable = generalfunctions.openSamplesAndTimestepsAsNumpyArraysAsNumpyArray( \
                  name,range(1,cfg.nrOfSamples+1),timeStepsWithStatsCalculated)
      numpy.save(name,aVariable)

  def createInstancesPremcloop(self):
    pass
    

  def createInstancesInitial(self):
    import generalfunctions

    timeStepsToReportAll = range(100,cfg.numberOfTimeSteps + 1,100)
    timeStepsToReportSome = range(3000,cfg.numberOfTimeSteps + 1,100)

    # class for exchange variables in initial and dynamic
    # introduced to make filtering possible
    
    self.d_exchangevariables=exchangevariables_weekly.ExchangeVariables( \
                                    timeStepsToReportSome, \
                                    cfg.exchange_report_rasters
                                    )

    # base level
    #deterministicDem=(ycoordinate(1)*0.4)
    deterministicDem=scalar('inputs_weekly/demini.map')
    #dem=deterministicDem+uniform(1)/2
    dem=deterministicDem
    baselevelRise=-0.0001
    self.d_baselevel=baselevel.Baselevel( \
                    generalfunctions.bottom(self.clone), \
                    deterministicDem, \
                    baselevelRise, \
                    self.timeStepDuration/(365.0*24.0), \
                    timeStepsToReportAll, \
                    cfg.baselevel_report_rasters)

    weatheringRateBareBedrock=0.0005
    weatheringExponentParameter=4.0
    self.d_bedrockweathering=bedrockweathering.BedrockWeathering( \
                             weatheringRateBareBedrock, \
                             weatheringExponentParameter, \
                             timeStepsToReportAll,\
                             cfg.bedrockweathering_report_rasters)
    steadyStateSoilDepth=self.d_bedrockweathering.steadyStateSoilDepth(0-baselevelRise)
    self.report(steadyStateSoilDepth,'sssd')

    # regolith
    regolithThickness=spatial(steadyStateSoilDepth)
    
    self.timeStepDurationRegolithInYears=1.0
    self.d_regolithdemandbedrock=regolith.RegolithDemAndBedrock(
                                 dem, \
                                 regolithThickness, \
                                 self.timeStepDurationRegolithInYears, \
                                 timeStepsToReportAll, \
                                 cfg.regolith_report_rasters)

    regolithThickness,demOfBedrock,dem,bedrockLdd,surfaceLdd=self.d_regolithdemandbedrock.getRegolithProperties()
    report(regolithThickness,'regIni.map')

    # biomass
    initialBiomass=2.0
    waterUseEfficiency=5.0
    maintenanceRate=0.5/(365.0*24.0)

    # original
    gamma=0.004 # runoff
    ## gamma zero
    #gamma=0.0001 # runoff
    ## gamma low
    #gamma=0.001 # runoff
    ## gamma medium
    #gamma=0.002 # runoff

    alpha = 0.4 # grazing
    dispersion = 0.01/(365.0*24)
    runoff = 0.0  # niet gebruikt ! (?) kan dus weg
    sdOfNoise = 0.000000000001
    LAIPerBiomass = 2.5
    self.d_biomassModifiedMay = biomassmodifiedmay.BiomassModifiedMay( \
                 initialBiomass, \
                 waterUseEfficiency, \
                 maintenanceRate, \
                 gamma, \
                 alpha, \
                 dispersion, \
                 runoff, \
                 sdOfNoise, \
                 LAIPerBiomass, \
                 self.timeStepDuration, \
                 timeStepsToReportAll, \
                 cfg.biomassmodifiedmay_report_rasters)

    # precipitation
    # scenario: original
    probabilityOfARainstorm=0.4
    durationOfRainstorm=cfg.theDurationOfRainstorm
    expectedRainfallIntensity=0.002
    gammaShapeParameter=100

    # scenario: higher intensity
    #probabilityOfARainstorm=0.4*0.75
    #durationOfRainstorm=cfg.theDurationOfRainstorm
    #expectedRainfallIntensity=0.002/0.75
    #gammaShapeParameter=100

    # scenario: much higher intensity
    #probabilityOfARainstorm=0.4*0.25
    #durationOfRainstorm=cfg.theDurationOfRainstorm
    #expectedRainfallIntensity=0.002/0.25
    #gammaShapeParameter=100

    # scenario: less rainstorms, longer duration
    #probabilityOfARainstorm=0.4*0.50
    #durationOfRainstorm=cfg.theDurationOfRainstorm/0.50
    #expectedRainfallIntensity=0.002
    #gammaShapeParameter=100

    # scenario: shorter rainstorm
    #probabilityOfARainstorm=0.4
    #durationOfRainstorm=cfg.theDurationOfRainstorm/2.0
    #expectedRainfallIntensity=0.002*2.0
    #gammaShapeParameter=100

    # scenario: more rainstorms (and also more rain in total)
    #probabilityOfARainstorm=0.999
    #durationOfRainstorm=cfg.theDurationOfRainstorm
    #expectedRainfallIntensity=0.002
    #gammaShapeParameter=100

    # scenario: all more 
    #probabilityOfARainstorm=0.4
    #durationOfRainstorm=cfg.theDurationOfRainstorm*2.0
    #expectedRainfallIntensity=0.004
    #gammaShapeParameter=100

    self.d_rainfalleventsfromgammadistribution= \
                                               rainfalleventsfromgammadistribution.RainfallEventsFromGammaDistribution( \
                                               probabilityOfARainstorm, \
                                               durationOfRainstorm, \
                                               expectedRainfallIntensity, \
                                               gammaShapeParameter, \
                                               timeStepsToReportAll,
                                               cfg.rainfalleventsfromgammadistribution_report_rasters)


    # interception
    initialLeafAreaIndex=initialBiomass*LAIPerBiomass
    initialInterceptionStore=scalar(0.000001)
    gapFraction,maximumInterceptionStore=calculateGapFractionAndMaxIntStoreFromLAI(initialLeafAreaIndex)

    self.d_interceptionuptomaxstore=interceptionuptomaxstore.InterceptionUpToMaxStore( \
                                    spatial(ldd(5)), \
                                    initialInterceptionStore, \
                                    maximumInterceptionStore, \
                                    gapFraction, \
                                    cfg.calculateUpstreamTotals, \
                                    durationOfRainstorm,
                                    timeStepsToReportAll,
                                    cfg.interception_report_rasters)

    # surface store
    initialSurfaceStore=scalar(0.0)
    maxSurfaceStore=scalar(0.0001)
    self.d_surfaceStore=surfacestore.SurfaceStore( \
                        initialSurfaceStore, \
                        maxSurfaceStore, \
                        durationOfRainstorm,
                        timeStepsToReportAll,
                        cfg.surfacestore_report_rasters)

    # infiltration
    bareSoilSaturatedConductivityFlux=scalar(0.0001)
    #maxSaturatedConductivityFluxFromVegetation=scalar(0.01)
    maxSaturatedConductivityFluxFromVegetation=scalar(0.1)
    biomassHalfSaturation=scalar(1.0)
    ksat=bareSoilSaturatedConductivityFlux
    self.d_infiltrationonlyksat=infiltrationonlyksat.InfiltrationOnlyKsat( \
                                    ksat, \
                                    bareSoilSaturatedConductivityFlux, \
                                    maxSaturatedConductivityFluxFromVegetation, \
                                    biomassHalfSaturation, \
                                    durationOfRainstorm,  \
                                    timeStepsToReportAll, \
                                    cfg.infiltration_report_rasters)


    # subsurface water
    # loam values from Niko Wanders, see mac disk articles/crittransGeom table
    #initialSoilMoistureFraction=scalar(0.03)
    initialSoilMoistureFraction=scalar(0.43)

    soilPorosityFraction= scalar(0.43)
    fieldCapacityFraction=scalar(0.22)
    limitingPointFraction=scalar(0.05)
    wiltingPointFraction= scalar(0.019)

    saturatedConductivityMetrePerDay=generalfunctions.mapuniformBounds( \
                                     2,8,scalar(12.5),cfg.createRealizations)

    self.d_subsurfaceWaterOneLayer=subsurfacewateronelayer.SubsurfaceWaterOneLayer(
                                   bedrockLdd,
                                   demOfBedrock,
                                   regolithThickness,
                                   initialSoilMoistureFraction,
                                   soilPorosityFraction,
                                   wiltingPointFraction,
                                   fieldCapacityFraction,
                                   limitingPointFraction,
                                   saturatedConductivityMetrePerDay,
                                   cfg.calculateUpstreamTotals, 
                                   self.timeStepDurationHours,
                                   timeStepsToReportAll,
                                   cfg.subsurface_report_rasters)

    # evapotranspiration
    beta=1.0
    maximumEvapotranspirationFlux=0.8/(365.0*24.0)
    self.d_evapotranspirationSimple=evapotranspirationsimple.EvapotranspirationSimple( \
                                    self.timeStepDuration, \
                                    beta, \
                                    maximumEvapotranspirationFlux, \
                                    timeStepsToReportAll, \
                                    cfg.evapotranspirationsimple_report_rasters) \

    # runoff
    self.d_runoffAccuthreshold=runoffaccuthreshold.RunoffAccuthreshold(
                               surfaceLdd,
                               durationOfRainstorm,
                               timeStepsToReportAll,
                               cfg.runoff_report_rasters)

    # soilwash
    plantHeightMetres=5.0
    stoneCoverFraction=0.1
    vegetationCoverOfSoilFraction=0.1
    manningsN=0.03 # 'original'

    # standard erosion scenario
    detachabilityOfSoilRaindrops=1.6 # 'original'  (used for all scenarios)
    detachabilityOfSoilRunoff=6.4 #'original'

    ## more erosion scenario
    #detachabilityOfSoilRaindrops=16
    #detachabilityOfSoilRunoff=64

    self.d_soilwashMMF=soilwashMMF.SoilWashMMF( \
                                              surfaceLdd,
                                              dem,
                                              durationOfRainstorm,
                                              plantHeightMetres,
                                              detachabilityOfSoilRaindrops,
                                              stoneCoverFraction, 
                                              detachabilityOfSoilRunoff,
                                              vegetationCoverOfSoilFraction,
                                              manningsN,
                                              soilPorosityFraction,
                                              timeStepsToReportAll,
                                              cfg.soilwashMMF_report_rasters)


    # creep
    diffusion=0.01
    self.d_creep=creep.Creep( \
                 dem, \
                 self.timeStepDurationRegolithInYears, \
                 diffusion, \
                 timeStepsToReportAll, \
                 cfg.creep_report_rasters)

  def reportComponentsDynamic(self):
    components =[ \
                 self.d_exchangevariables, \
                 self.d_evapotranspirationSimple, \
                 self.d_regolithdemandbedrock, \
                 self.d_bedrockweathering,
                 self.d_baselevel, \
                 self.d_rainfalleventsfromgammadistribution , \
                 self.d_interceptionuptomaxstore, \
                 self.d_surfaceStore, \
                 self.d_infiltrationonlyksat, \
                 self.d_runoffAccuthreshold, \
                 self.d_subsurfaceWaterOneLayer, \
                 self.d_soilwashMMF, \
                 self.d_creep,
                 self.d_biomassModifiedMay
                 ] 

    for component in components:
      component.reportAsMaps(self.currentSampleNumber(), self.currentTimeStep())


  def printMemberVariables(self):
    import generalfunctions
    components =[ \
                 self.d_exchangevariables, \
                 self.d_interceptionuptomaxstore, \
                 self.d_surfaceStore, \
                 self.d_infiltrationonlyksat, \
                 self.d_runoffAccuthreshold, \
                 self.d_subsurfaceWaterOneLayer
                 ] 

    for component in components:
      generalfunctions.printMemberVariables(component)
 
    
  def printComponentsDynamic(self):
    self.d_dateTimePCRasterPython.printit()

 
  def initializeTime(self,startTimeYear, startTimeMonth, startTimeDay, timeStepDurationHours):
    startTime=datetime.datetime(year=startTimeYear, month=startTimeMonth, day=startTimeDay)
    self.timeStepDurationHours = timeStepDurationHours
    self.timeStepDatetimeFormat=datetime.timedelta(hours=self.timeStepDurationHours)
    self.d_dateTimePCRasterPython=datetimePCRasterPython.DatetimePCRasterPython \
                                  (startTime, self.timeStepDatetimeFormat)


  def checkBudgets(self,currentSampleNumber,currentTimeStep):

    increaseInPrecipitationStore=0.0-self.d_exchangevariables.cumulativePrecipitation
    report(increaseInPrecipitationStore,generateNameST('incP', currentSampleNumber, currentTimeStep))

    increaseInInterceptionStore=self.d_interceptionuptomaxstore.budgetCheck(currentSampleNumber, currentTimeStep)
    report(increaseInInterceptionStore,generateNameST('incI', currentSampleNumber, currentTimeStep))

    increaseInSurfaceStore=self.d_surfaceStore.budgetCheck(currentSampleNumber, currentTimeStep)
    report(increaseInSurfaceStore,generateNameST('incS', currentSampleNumber, currentTimeStep))
    increaseInSurfaceStoreQM=catchmenttotal(increaseInSurfaceStore,self.ldd)*cellarea()
    report(increaseInSurfaceStoreQM,generateNameST('testb', currentSampleNumber, currentTimeStep))

    # let op: infiltration store is directly passed to subsurface store, thus is not a real store
    increaseInInfiltrationStore=self.d_infiltrationonlyksat.budgetCheck(currentSampleNumber, currentTimeStep)

    increaseInSubSurfaceWaterStore,lateralFlowInSubsurfaceStore,abstractionFromSubSurfaceWaterStore= \
                             self.d_subsurfaceWaterOneLayer.budgetCheck(currentSampleNumber , currentTimeStep)
    increaseInSubSurfaceStoreQM=catchmenttotal(increaseInSubSurfaceWaterStore,self.ldd)*cellarea()

    increaseInRunoffStoreCubicMetresInUpstreamArea=self.d_runoffAccuthreshold.budgetCheck()

    totalIncreaseInStoresCubicMetresInUpstreamArea=0.0
    stores=[increaseInPrecipitationStore, increaseInInterceptionStore, increaseInSurfaceStore, increaseInSubSurfaceWaterStore]
    for store in stores:
      increaseInStoreCubicMetresInUpstreamArea=catchmenttotal(store,self.ldd)*cellarea()
      totalIncreaseInStoresCubicMetresInUpstreamArea=totalIncreaseInStoresCubicMetresInUpstreamArea+ \
                                                     increaseInStoreCubicMetresInUpstreamArea
    
    report(totalIncreaseInStoresCubicMetresInUpstreamArea,generateNameST('inSt',currentSampleNumber,currentTimeStep))
    report(increaseInRunoffStoreCubicMetresInUpstreamArea ,generateNameST('inRu',currentSampleNumber,currentTimeStep))
    report(catchmenttotal(self.d_exchangevariables.upwardSeepageFlux,self.ldd)*cellarea(),generateNameST('inSe',currentSampleNumber,currentTimeStep))
    # total budget is total increase in stores plus the upward seepage flux for each ts that is passed to the next
    # timestep and thus not taken into account in the current timestep budgets
    budget=totalIncreaseInStoresCubicMetresInUpstreamArea+increaseInRunoffStoreCubicMetresInUpstreamArea + \
           lateralFlowInSubsurfaceStore*cellarea() + catchmenttotal(abstractionFromSubSurfaceWaterStore,self.ldd)*cellarea() + \
           catchmenttotal(self.d_exchangevariables.upwardSeepageFlux,self.ldd)*cellarea()
    report(budget,generateNameST('B-tot', currentSampleNumber, currentTimeStep))
    budgetRel=budget/increaseInRunoffStoreCubicMetresInUpstreamArea
    report(budgetRel,generateNameST('B-rel', currentSampleNumber, currentTimeStep))

myModel = CatchmentModel()
dynamicModel = DynamicFramework(myModel, cfg.numberOfTimeSteps)
mcModel = MonteCarloFramework(dynamicModel, cfg.nrOfSamples)
mcModel.setForkSamples(True,10)
mcModel.run()
