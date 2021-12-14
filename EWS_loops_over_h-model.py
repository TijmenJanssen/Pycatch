# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 15:22:12 2021

@author: TEJan
"""     

import EWSPy as ews
from pcraster import *
import numpy as np
import os
import time
from scipy import ndimage

import configuration_weekly as cfg
import NULL_models_timeseries as temp_NULL
import NULL_models_spatial as spat_NULL
import EWS_StateVariables as ews_sv



def ews_calculations(variable, run='/h01', path='./1/', timer_on=False):
    ## Temporal EWS calculations ##
    
    if variable.temporal:
        if temporal_ews_present:
            print(f"Started temporal EWS calculations for {variable.name}")

            ## Start timer if set to True ##
            if timer_on:
                start_time = time.time()

            ## Timeseries file loading ##
            state_variable_timeseries = []
            if variable.datatype == 'numpy':
                file_name = ews.file_name_str(variable.name, temporal_ews_interval)
                state_variable_timeseries = np.loadtxt(run + file_name + ".numpy.txt")
            else:
                print(f"Datatype for {variable.name} currently not supported.")

            ## Splitting timeseries into (overlapping) windows ##
            stack_of_windows = time_series2time_windows(state_variable_timeseries, variable.window_size,
                                                        variable.window_overlap)

            ## EWS calculations ###
            # Temporal mean #
            fpath = os.path.join(path + variable.name + '.t.mn')
            np.savetxt(fpath + '.numpy.txt', ews.temporal_mean(stack_of_windows))

            # Temporal std #
            fpath = os.path.join(path + variable.name + '.t.std')
            np.savetxt(fpath + '.numpy.txt', ews.temporal_std(stack_of_windows))

            # Temporal var #
            fpath = os.path.join(path + variable.name + '.t.var')
            np.savetxt(fpath + '.numpy.txt', ews.temporal_var(stack_of_windows))

            # Temporal cv #
            fpath = os.path.join(path + variable.name + '.t.cv')
            np.savetxt(fpath + '.numpy.txt', ews.temporal_cv(stack_of_windows))

            # Temporal skw #
            fpath = os.path.join(path + variable.name + '.t.skw')
            np.savetxt(fpath + '.numpy.txt', ews.temporal_skw(stack_of_windows))

            # Temporal krt #
            fpath = os.path.join(path + variable.name + '.t.krt')
            np.savetxt(fpath + '.numpy.txt', ews.temporal_krt(stack_of_windows))

            # # Temporal dfa # TODO - returns 3 values (tuple?)
            # fpath = os.path.join(path + variable.name + '.t.dfa')
            # np.savetxt(fpath + '.numpy.txt', ews.temporal_dfa(stack_of_windows))

            # Temporal autocorr. #
            fpath = os.path.join(path + variable.name + '.t.acr')
            np.savetxt(fpath + '.numpy.txt', ews.temporal_autocorrelation(stack_of_windows))

            # Temporal AR1 #
            fpath = os.path.join(path + variable.name + '.t.AR1')
            np.savetxt(fpath + '.numpy.txt', ews.temporal_AR1(stack_of_windows))

            # Temporal return rate #
            fpath = os.path.join(path + variable.name + '.t.rr')
            np.savetxt(fpath + '.numpy.txt', ews.temporal_returnrate(stack_of_windows))

            # Temporal cond. het. #
            fpath = os.path.join(path + variable.name + '.t.coh')
            np.savetxt(fpath + '.numpy.txt', ews.temporal_cond_het(stack_of_windows))

            ## End timer if set to True##
            if timer_on:
                end_time = time.time()
                print("Elapsed time for temporal EWS calculations equals:", end_time - start_time, '\n')

        elif temporal_ews_present == False:
            print(
                f"Mean timeseries data == False in configuration_weekly.py, could not calculate EWS for {variable.name}.")

    ## Temporal EWS calculations ##
    if variable.spatial:
        if spatial_ews_present:
            print(f"Started spatial EWS calculations for {variable.name}")

            ## Start timer if set to True ##
            if timer_on:
                start_time = time.time()

            ## Spatial maps file loading ##
            state_variable_snapshots = [0.0] * len(spatial_ews_interval)
            if variable.datatype == 'numpy':
                for k, timestep in enumerate(spatial_ews_interval):
                    file_name = ews.file_name_str(variable.name, timestep)
                    state_variable_snapshots[k] = np.loadtxt(run + file_name + 'numpy.txt')
            if variable.datatype == 'map':
                for k, timestep in enumerate(spatial_ews_interval):
                    file_name = ews.file_name_str(variable.name, timestep)
                    state_variable_snapshots[k] = pcr2numpy(readmap(run + file_name), np.NaN)
            else:
                print(f"Datatype for {variable.name} currently not supported.")
            state_variable_snapshots = np.asarray(state_variable_snapshots)

            ## EWS calculations ##
            # Spatial mean #
            fpath = os.path.join(path + variable.name + '.s.mn')
            np.savetxt(fpath + '.numpy.txt', ews.spatial_mean(state_variable_snapshots))

            # Spatial std #
            fpath = os.path.join(path + variable.name + '.s.std')
            np.savetxt(fpath + '.numpy.txt', ews.spatial_std(state_variable_snapshots))

            # Spatial var #
            fpath = os.path.join(path + variable.name + '.s.var')
            np.savetxt(fpath + '.numpy.txt', ews.spatial_var(state_variable_snapshots))

            # Spatial skw #
            fpath = os.path.join(path + variable.name + '.s.skw')
            np.savetxt(fpath + '.numpy.txt', ews.spatial_skw(state_variable_snapshots))

            # Spatial krt #
            # fpath = os.path.join(path + variable.name + '.s.krt')
            # np.savetxt(fpath + '.numpy.txt', ews.spatial_krt(state_variable_snapshots))

            # Spatial correlation (Moran's I) #
            fpath = os.path.join(path + variable.name + '.s.mI')
            np.savetxt(fpath + '.numpy.txt', ews.spatial_corr(state_variable_snapshots))

            # # spatial DFT #
            # fpath = os.path.join(path + variable.name + '.s.dft')
            # np.savetxt(fpath + '.numpy.txt', ews.spatial_DFT(state_variable_snapshots))

            ## End timer if set to True##
            if timer_on:
                end_time = time.time()
                print("Elapsed time for temporal EWS calculations equals:", end_time - start_time, '\n')

        elif spatial_ews_present == False:
            print(f"Map data == False in configuration_weekly.py, could not calculate EWS for {variable.name}.")
 



# Temporal mean 
fpath = os.path.join(path + variable.name + '.h.t.mn')
if run == 'h01':
    np.savetxt(fpath + '.numpy.txt', ews.temporal_mean(stack_of_windows))
else:
    a = np.loadtxt(fpath + '.numpy.txt')
    a = np.append(a, ews.temporal_mean(stack_of_windows))
    np.savetxt(fpath + '.numpy.txt', ews.temporal_mean((stack_of_windows)))
    
           
# Temporal std    
fpath = os.path.join(path + variable.name + '.h.t.std')
if run == 'h01':
    np.savetxt(fpath + '.numpy.txt', ews.temporal_std(stack_of_windows))
else:
    a = np.loadtxt(fpath + '.numpy.txt')
    a = np.append(a, ews.temporal_std(stack_of_windows)) 
    np.savetxt(fpath + '.numpy.txt', a)
    
    
# Temporal var       
fpath = os.path.join(path + variable.name + '.h.t.var')
if run == 'h01':
    np.savetxt(fpath + '.numpy.txt', ews.temporal_var(stack_of_windows))
else:
    a = np.loadtxt(fpath + '.numpy.txt')
    a = np.append(a, ews.temporal_var(stack_of_windows)) 
    np.savetxt(fpath + '.numpy.txt', a)  
    
    
 # Temporal cv                
fpath = os.path.join(path + variable.name + '.h.t.cv')
if run == 'h01':
    np.savetxt(fpath + '.numpy.txt', ews.temporal_cv(stack_of_windows))
else:
    a = np.loadtxt(fpath + '.numpy.txt')
    a = np.append(a, ews.temporal_cv(stack_of_windows)) 
    np.savetxt(fpath + '.numpy.txt', a) 

 # Temporal skw          
fpath = os.path.join(path + variable.name + '.h.t.skw')
if run == 'h01':
    np.savetxt(fpath + '.numpy.txt', ews.temporal_skw(stack_of_windows))
else:
    a = np.loadtxt(fpath + '.numpy.txt')
    a = np.append(a, ews.temporal_skw(stack_of_windows)) 
    np.savetxt(fpath + '.numpy.txt', a) 

  # Temporal krt   
fpath = os.path.join(path + variable.name + '.h.t.krt')
if run == 'h01':
    np.savetxt(fpath + '.numpy.txt', ews.temporal_krt(stack_of_windows))
else:
    a = np.loadtxt(fpath + '.numpy.txt')
    a = np.append(a, ews.temporal_krt(stack_of_windows)) 
    np.savetxt(fpath + '.numpy.txt', a) 

        # # # Temporal dfa # TODO - returns 3 values (tuple?)
        # # fpath = os.path.join(path + variable.name + '.t.dfa')
        # # np.savetxt(fpath + '.numpy.txt', ews.temporal_dfa(stack_of_windows))


# Temporal autocorr.         
fpath = os.path.join(path + variable.name + '.h.t.acr')
if run == 'h01':
    np.savetxt(fpath + '.numpy.txt', ews.temporal_acr(stack_of_windows))
else:
    a = np.loadtxt(fpath + '.numpy.txt')
    a = np.append(a, ews.temporal_acr(stack_of_windows)) 
    np.savetxt(fpath + '.numpy.txt', a) 

 # Temporal AR1         
fpath = os.path.join(path + variable.name + '.h.t.AR1')
if run == 'h01':
    np.savetxt(fpath + '.numpy.txt', ews.temporal_AR1(stack_of_windows))
else:
    a = np.loadtxt(fpath + '.numpy.txt')
    a = np.append(a, ews.temporal_AR1(stack_of_windows)) 
    np.savetxt(fpath + '.numpy.txt', a) 

        # # Temporal return rate #    
fpath = os.path.join(path + variable.name + '.h.t.rr')
if run == 'h01':
    np.savetxt(fpath + '.numpy.txt', ews.temporal_returnrate(stack_of_windows))
else:
    a = np.loadtxt(fpath + '.numpy.txt')
    a = np.append(a, ews.temporal_returnrate(stack_of_windows)) 
    np.savetxt(fpath + '.numpy.txt', a) 

# Temporal cond. het. 
        

fpath = os.path.join(path + variable.name + '.h.t.coh')
if run == 'h01':
    np.savetxt(fpath + '.numpy.txt', ews.temporal_cond_het(stack_of_windows))
else:
    a = np.loadtxt(fpath + '.numpy.txt')
    a = np.append(a, ews.temporal_cond_het(stack_of_windows)) 
    np.savetxt(fpath + '.numpy.txt', a)     
    
    
    
    