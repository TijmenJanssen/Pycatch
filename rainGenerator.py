# # -*- coding: utf-8 -*-
# """
# Created on Tue Mar  1 13:06:15 2022

# @author: TEJan
# """

import math
import numpy as np
import hourly_configuration as cfg
import random
import pcraster as pcr
import matplotlib.pyplot as plt


def mapgamma(shapeParameter):
  '''Returns a realization from the gamma distribution with a mean of one
  shapeParameter is a Python floating point
  return value is a Python floating point
  '''

  scaleParameter = 1.0 / shapeParameter
  realization = random.gammavariate(shapeParameter, scaleParameter)
  return realization


hours_per_week = 7 * 24
nr_of_weeks = math.ceil((cfg.number_of_timesteps_hourly + 1)/ hours_per_week)
            
# 1 if this week contains rain, else 0 (uniform distribution)
rairray = np.random.rand(nr_of_weeks) 
rairray[rairray >= (1-cfg.probabilityOfRainstorm)] =   1
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
rairray_timeseries = [0.0] * len(rairray_hours_flat)
for k, event in enumerate(rairray_hours_flat):
    if event == 1.0:
        rairray_timeseries[k] = cfg.expectedRainfallIntensity * mapgamma(cfg.gammaShapeParameter)
    else:
        rairray_timeseries[k] = 0.0

#rairray_timeseries = rairray_hours_flat * cfg.expectedRainfallIntensity * generalfunctions.mapgamma(cfg.gammaShapeParameter)
for i in range(1, cfg.rainstormDuration):
    rairray_timeseries += np.roll(rairray_timeseries, i)
    
np.savetxt('./rain_array.numpy.txt', rairray_timeseries)

rain = np.loadtxt('./rain_array.numpy.txt') 




plt.plot(np.arange(0, len(rain)), rain)
plt.xlabel("time (h)")
plt.ylabel("rain")

plt.show()