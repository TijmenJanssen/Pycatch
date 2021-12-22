import EWSPy as ews
from pcraster import *
import numpy as np
import os
import time
from scipy import ndimage

import EWS_main_configuration as cfg
# import NULL_models_timeseries_weekly as temp_NULL
# import NULL_models_spatial_weekly as spat_NULL
import EWS_StateVariables as ews_sv

#TIJMEN: IN DIT DOCUMENT MAP AANPASSEN VOOR RUN

### User input ###

## State variables for EWS ##
variables = ews_sv.variables_hourly  # State variables present in EWS_StateVariables can be added through configuration_weekly

## Generate dummy datasets for Kendall tau? ##
generate_dummy_datasets = False
save_detrended_data = True # Temporal only
method_1 = True
method_2 = True
method_3 = True
nr_generated_datasets = 1

### End user input ###

### Information from configuration weekly ###
## Realizations/MC runs ##
realizations = cfg.nrOfSamples
# realizations = 1 # for test cases

## Timesteps, intervals ##
spatial_ews_present = cfg.map_data
spatial_ews_interval = np.arange(cfg.interval_map_snapshots, cfg.number_of_timesteps_hourly + cfg.interval_map_snapshots,
                                 cfg.interval_map_snapshots)

temporal_ews_present = cfg.mean_timeseries_data
temporal_ews_interval = cfg.number_of_timesteps_hourly  # the interval is defined as t=0 to the last timestep.


### Functions ###

def time_series2time_windows(time_series, window_size=100, window_overlap=0):
    return np.array([time_series[i:i + window_size] for i in range(0, len(time_series), window_size - window_overlap)])


def generate_datasets(variable, path='./1/', nr_realizations=1, detrending_temp='None', sigma=50,
                      method1=False, method2=False, method3=False):
    if variable.temporal:
        ## Load data ##
        state_variable_timeseries = []
        if variable.datatype == 'numpy':
            file_name = ews.file_name_str(variable.name, temporal_ews_interval)
            state_variable_timeseries = np.loadtxt(path + file_name + ".numpy.txt")
        else:
            print(f"Datatype for {variable.name} currently not supported.")

        ## Detrending: 'None', 'Gaussian' ##
        if detrending_temp == 'None':
            state_variable_timeseries = state_variable_timeseries
        if detrending_temp == 'Gaussian':  # TODO - Multiple sigmas?
            state_variable_timeseries = state_variable_timeseries - ndimage.gaussian_filter1d(state_variable_timeseries,
                                                                                              sigma)
            if save_detrended_data:
                generated_number_length = 4
                if len(str(realizations)) > 4:
                    generated_number_length = len(str(realizations))
                generated_number_string = 'dtr' + str(realization).zfill(generated_number_length)
                dir_name = os.path.join(path + generated_number_string)
                fname = ews.file_name_str(variable.name, cfg.number_of_timesteps_hourly)
                fpath = os.path.join(dir_name, fname)
                np.savetxt(fpath + '.numpy.txt', state_variable_timeseries)

        ## Generate dummy datasets ##
        if method1:
            temp_NULL.method1_(state_variable_timeseries, realizations=nr_realizations, path=path,
                               file_name=variable.name)
        if method2:
            temp_NULL.method2_(state_variable_timeseries, realizations=nr_realizations, path=path,
                               file_name=variable.name)
        if method3:
            temp_NULL.method3_(state_variable_timeseries, realizations=nr_realizations, path=path,
                               file_name=variable.name)

    if variable.spatial:
        ## Load data ##
        state_variable_snapshots = [0.0] * len(spatial_ews_interval)
        if variable.datatype == 'numpy':
            for k, timestep in enumerate(spatial_ews_interval):
                file_name = ews.file_name_str(variable.name, timestep)
                state_variable_snapshots[k] = np.loadtxt(path + file_name + 'numpy.txt')
        if variable.datatype == 'map':
            for k, timestep in enumerate(spatial_ews_interval):
                file_name = ews.file_name_str(variable.name, timestep)
                state_variable_snapshots[k] = pcr2numpy(readmap(path + file_name), np.NaN)
        else:
            print(f"Datatype for {variable.name} currently not supported.")
        state_variable_snapshots = np.asarray(state_variable_snapshots)

        ## Generate dummy datasets ##
        if method1:
            spat_NULL.method1_(state_variable_snapshots, realizations=nr_realizations, path=path,
                               file_name=variable.name)
        if method2:
            spat_NULL.method2_(state_variable_snapshots, realizations=nr_realizations, path=path,
                               file_name=variable.name)
        if method3:
            spat_NULL.method3_(state_variable_snapshots, realizations=nr_realizations, path=path,
                               file_name=variable.name)


def ews_calculations_generated_datasets(variable, path='./1/', nr_realizations=1, timer_on=False, method1=False,
                                        method2=False, method3=False):
    generated_number_length = 4
    if len(str(realizations)) > 4:
        generated_number_length = len(str(realizations))

    if method1:
        for realization in range(nr_realizations):
            generated_number_string = 'm1g' + str(realization).zfill(generated_number_length) + '/'
            dir_name = os.path.join(path + generated_number_string)
            ews_calculations(variable, path=dir_name, timer_on=timer_on)
    if method2:
        for realization in range(nr_realizations):
            generated_number_string = 'm2g' + str(realization).zfill(generated_number_length) + '/'
            dir_name = os.path.join(path + generated_number_string)
            ews_calculations(variable, path=dir_name, timer_on=timer_on)
    if method3:
        for realization in range(nr_realizations):
            generated_number_string = 'm3g' + str(realization).zfill(generated_number_length) + '/'
            dir_name = os.path.join(path + generated_number_string)
            ews_calculations(variable, path=dir_name, timer_on=timer_on)


def ews_calculations(variable, path='./1/', timer_on=False):
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
                state_variable_timeseries = np.loadtxt(path + file_name + ".numpy.txt")
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
                    state_variable_snapshots[k] = np.loadtxt(path + file_name + 'numpy.txt')
            if variable.datatype == 'map':
                for k, timestep in enumerate(spatial_ews_interval):
                    file_name = ews.file_name_str(variable.name, timestep)
                    state_variable_snapshots[k] = pcr2numpy(readmap(path + file_name), np.NaN)
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
            fpath = os.path.join(path + variable.name + '.s.krt')
            np.savetxt(fpath + '.numpy.txt', ews.spatial_krt(state_variable_snapshots))

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


### Running the functions for given variables ###

for i in range(cfg.stepsTotal):
    fpath = str(i).zfill(2)
    for variable in variables:
        ews_calculations(variable, path=f'./h{fpath}/', timer_on=True)
        if generate_dummy_datasets:
            generate_datasets(variable, path=f'./h{fpath}/', nr_realizations=nr_generated_datasets,
                              detrending_temp='Gaussian', method1=method_1, method2=method_2, method3=method_3)
            ews_calculations_generated_datasets(variable, path=f'./h{path}/',
                                                nr_realizations=nr_generated_datasets,
                                                timer_on=True, method1=method_1, method2=method_2, method3=method_3)
            
            
# end_time = time.time()
# print('Elapsed timee for Pychatch h-model loops is', ((end_time - start_time)/60), 'minutes' )