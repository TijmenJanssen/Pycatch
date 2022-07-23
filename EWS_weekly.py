# -*- coding: utf-8 -*-
"""
EWS - Early Warning Signals
EWS_weekly

@authors: KoenvanLoon & TijmenJanssen
"""

from pcraster import *
import numpy as np
import os
import time
from scipy import ndimage

import EWSPy as ews
import hourly_configuration as cfg
# import NULL_models_timeseries_weekly as temp_NULL
# import NULL_models_spatial_weekly as spat_NULL
import EWS_StateVariables as ews_sv

### User input ###

## State variables for EWS ##
variables = ews_sv.variables_hourly  # State variables present in EWS_StateVariables can be added through configuration


### End user input ###

### Information from configuration weekly ###

## Timesteps, intervals ##

if not cfg.cutoff:
    spatial_ews_interval = np.arange(cfg.interval_map_snapshots, cfg.number_of_timesteps_hourly +
                                     cfg.interval_map_snapshots, cfg.interval_map_snapshots)
elif cfg.cutoff:
    spatial_ews_interval = np.arange(cfg.interval_map_snapshots, cfg.cutoff_point + cfg.interval_map_snapshots,
                                     cfg.interval_map_snapshots)

temporal_ews_interval = cfg.number_of_timesteps_hourly  # the interval is defined as t=0 : the last timestep.


### Functions ###

def time_series2time_windows(timeseries, window_size=100, window_overlap=0):
    actual_window_overlap = window_size - window_overlap
    sh = (timeseries.size - window_size + 1, window_size)
    st = timeseries.strides * 2
    if window_overlap != 0:
        view = np.lib.stride_tricks.as_strided(timeseries, strides=st, shape=sh)[0::actual_window_overlap]
    elif window_overlap == 0:
        view = np.lib.stride_tricks.as_strided(timeseries, strides=st, shape=sh)[0::window_size]
    return view


def generate_datasets_main(variable, state_variable, method, nr_realizations=1, path='./1/'):
    print(f"Started generating dataset(s) for {variable.name} using {method.__name__}")
    detrended_data = method(state_variable, realizations=nr_realizations, path=path, file_name=variable.name)
    print(f"Finished generating dataset(s) for {variable.name} using {method.__name__} \n")
    if method.__name__ == temp_NULL.detrend_.__name__:
        return detrended_data


def generate_datasets(variable, path='./1/', nr_realizations=1, detrending_temp='None', sigma=50,
                      method1=False, method2=False, method3=False):
    if variable.temporal:
        state_variable_timeseries, files_present = temporal_data_file_loading(variable, path=path)

        if files_present:
            if state_variable_timeseries.ndim == 1:
                # Detrending: 'None', 'Gaussian'
                state_variable_timeseries = generate_datasets_main(variable, state_variable_timeseries, temp_NULL.detrend_, nr_realizations=nr_realizations, path=path)
                # Generate dummy datasets
                if method1:
                    generate_datasets_main(variable, state_variable_timeseries, temp_NULL.method1_,
                                           nr_realizations=nr_realizations, path=path)
                if method2:
                    generate_datasets_main(variable, state_variable_timeseries, temp_NULL.method2_,
                                           nr_realizations=nr_realizations, path=path)
                if method3:
                    generate_datasets_main(variable, state_variable_timeseries, temp_NULL.method3_,
                                           nr_realizations=nr_realizations, path=path)
            else:
                print(f"Multiple dimensions are currently not supported for generated datasets, so no datasets are being "
                      f"generated for {variable.name}.")

    if variable.spatial:
        state_variable_snapshots, files_present = spatial_data_file_loading(variable, path=path)
        if files_present:
            state_variable_snapshots = np.asarray(state_variable_snapshots)

            # Generate dummy datasets
            if method1:
                generate_datasets_main(variable, state_variable_snapshots, spat_NULL.method1_,
                                       nr_realizations=nr_realizations, path=path)
            if method2:
                generate_datasets_main(variable, state_variable_snapshots, spat_NULL.method2_,
                                       nr_realizations=nr_realizations, path=path)
            if method3:
                generate_datasets_main(variable, state_variable_snapshots, spat_NULL.method3_,
                                       nr_realizations=nr_realizations, path=path)


def ews_calculations_generated_datasets_init(variable, path='./1/', nr_realizations=1, timer_on=False, method1=False,
                                             method2=False, method3=False):
    generated_number_length = ews.generated_number_length(nr_realizations)

    if cfg.save_detrended_data and variable.temporal:
        ews_calculations_generated_datasets_main(variable, 'dtr', gen_nr_len=generated_number_length, path=path,
                                                 nr_realizations=1, timer_on=timer_on)
    if method1:
        ews_calculations_generated_datasets_main(variable, 'm1g', gen_nr_len=generated_number_length, path=path,
                                                 nr_realizations=nr_realizations, timer_on=timer_on)
    if method2:
        ews_calculations_generated_datasets_main(variable, 'm2g', gen_nr_len=generated_number_length, path=path,
                                                 nr_realizations=nr_realizations, timer_on=timer_on)
    if method3:
        ews_calculations_generated_datasets_main(variable, 'm3g', gen_nr_len=generated_number_length, path=path,
                                                 nr_realizations=nr_realizations, timer_on=timer_on)


def ews_calculations_generated_datasets_main(variable, method, gen_nr_len=4, path='./1/', nr_realizations=1, timer_on=False):
    for realization in range(nr_realizations):
        generated_number_string = method + str(realization).zfill(gen_nr_len) + '/'
        dir_name = os.path.join(path + generated_number_string)
        ews_calculations_init(variable, path=dir_name, timer_on=timer_on)


def temporal_data_file_loading(variable, path='./1/'):
    state_variable_timeseries = []
    EWS_calculations = True
    if variable.datatype == 'numpy':
        file_name = ews.file_name_str(variable.name, temporal_ews_interval)
        if os.path.exists(path + file_name + ".numpy.txt"):
            state_variable_timeseries = np.loadtxt(path + file_name + ".numpy.txt")
        else:
            print(f"{file_name}.numpy.txt not found in {path}")
            EWS_calculations = False
    else:
        print(f"Datatype for {variable.name} currently not supported.")
        EWS_calculations = False

    return state_variable_timeseries, EWS_calculations


def window_stacker(variable, state_variable_timeseries):
    nr_dim = state_variable_timeseries.ndim
    if nr_dim == 1:
        if cfg.cutoff:
            state_variable_timeseries = state_variable_timeseries[:cfg.cutoff_point]
        stack_of_windows = time_series2time_windows(state_variable_timeseries, variable.window_size,
                                                    variable.window_overlap)
        stack_x, stack_y = np.asarray(stack_of_windows).shape
    else:
        stack_of_windows = [0.0] * np.asarray(state_variable_timeseries).shape[1]
        for k, timeseries in enumerate(state_variable_timeseries.T):
            if cfg.cutoff:
                stack_of_windows[k] = time_series2time_windows(timeseries[:cfg.cutoff_point],
                                                               variable.window_size, variable.window_overlap)
            elif not cfg.cutoff:
                stack_of_windows[k] = time_series2time_windows(timeseries, variable.window_size,
                                                               variable.window_overlap)
        stack_x, stack_y, stack_z = np.asarray(stack_of_windows).shape
        stack_of_windows = np.asarray(stack_of_windows).reshape(-1, stack_z)
    return stack_of_windows, nr_dim, stack_x, stack_y


def spatial_data_file_loading(variable, path='./1/'):
    state_variable_snapshots = [0.0] * len(spatial_ews_interval)
    EWS_calculations = True
    if variable.datatype == 'numpy':
        for k, timestep in enumerate(spatial_ews_interval):
            file_name = ews.file_name_str(variable.name, timestep)
            if os.path.exists(path + file_name + ".numpy.txt"):
                state_variable_snapshots[k] = np.loadtxt(path + file_name + 'numpy.txt')
            else:
                print(f"{file_name}.numpy.txt not found in {path}.")
                EWS_calculations = False

    if variable.datatype == 'map':
        for k, timestep in enumerate(spatial_ews_interval):
            file_name = ews.file_name_str(variable.name, timestep)
            if os.path.exists(path + file_name):
                state_variable_snapshots[k] = pcr2numpy(readmap(path + file_name), np.NaN)
            else:
                print(f"{file_name} not found in {path}.")
                EWS_calculations = False
    else:
        print(f"Datatype for {variable.name} currently not supported.")
        EWS_calculations = False

    return state_variable_snapshots, EWS_calculations


def ews_calc_and_save(variable, data, method_name, method_function, path='./1/', nr_dim=1, stack_x=1, stack_y=1):
    fpath = os.path.join(path + variable.name + method_name)
    signal = method_function(data)
    if nr_dim > 1:
        signal = signal.reshape(stack_x, stack_y)
    np.savetxt(fpath + '.numpy.txt', signal)


def temporal_ews_calculations(variable, state_variable_timeseries, path='./1/'):
    stack_of_windows, nr_dim, stack_x, stack_y = window_stacker(variable, state_variable_timeseries)

    ews_calc_and_save(variable, stack_of_windows, '.t.mn', ews.temporal_mean, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.std', ews.temporal_std, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.var', ews.temporal_var, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.cv', ews.temporal_cv, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.skw', ews.temporal_skw, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.krt', ews.temporal_krt, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.acr', ews.temporal_autocorrelation, path=path, nr_dim=nr_dim,
                      stack_x=stack_x, stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.AR1', ews.temporal_AR1, path=path, nr_dim=nr_dim, stack_x=stack_x,
                      stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, '.t.rr', ews.temporal_returnrate, path=path, nr_dim=nr_dim,
                      stack_x=stack_x, stack_y=stack_y)
    ews_calc_and_save(variable, stack_of_windows, ".t.sum", ews.temporal_sum, path = path,  nr_dim=nr_dim,
                      stack_x=stack_x, stack_y=stack_y)

    # Temporal dfa TODO - returns 3-4 values, save only 1?
    fpath = os.path.join(path + variable.name + '.t.dfa')
    _, _, _, temporal_statistic = ews.temporal_dfa(stack_of_windows, window_size=variable.window_size)
    # scales, fluct, coeff, propagator (== temporal statistic)
    if nr_dim > 1:
        temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
    np.savetxt(fpath + '.numpy.txt', temporal_statistic)

    # Temporal cond. het. TODO - returns 2 values, save only 1?
    fpath = os.path.join(path + variable.name + '.t.coh')
    save_p = True
    if save_p and nr_dim == 1:
        temporal_statistic = [[0.0], [0.0]]
        statistic, p_val = ews.temporal_cond_het(stack_of_windows)
        temporal_statistic[0] = statistic
        temporal_statistic[1] = p_val
    else:
        temporal_statistic, _ = ews.temporal_cond_het(stack_of_windows)  # _ is the p-value of the test, not saved
        if nr_dim > 1:
            temporal_statistic = temporal_statistic.reshape(stack_x, stack_y)
    np.savetxt(fpath + '.numpy.txt', temporal_statistic)


def spatial_ews_calculations(variable, state_variable_maps, path='./1/'):
    ews_calc_and_save(variable, state_variable_maps, '.s.mn', ews.spatial_mean, path=path)
    ews_calc_and_save(variable, state_variable_maps, '.s.std', ews.spatial_std, path=path)
    ews_calc_and_save(variable, state_variable_maps, '.s.var', ews.spatial_var, path=path)
    ews_calc_and_save(variable, state_variable_maps, '.s.skw', ews.spatial_skw, path=path)
    ews_calc_and_save(variable, state_variable_maps, '.s.krt', ews.spatial_krt, path=path)
    ews_calc_and_save(variable, state_variable_maps, '.s.mI', ews.spatial_corr, path=path)
    # ews_calc_and_save(variable, state_variable_maps, '.s.dft', ews.spatial_DFT, path=path)


def ews_calculations_init(variable, path='./1/', timer_on=False):
    if variable.temporal:
        if cfg.mean_timeseries_data:
            ews_calculations_main(variable, temporal_data_file_loading, temporal_ews_calculations, path=path,
                                  timer_on=timer_on)
        elif not cfg.mean_timeseries_data:
            print(f"Mean timeseries data == False in configuration_weekly.py, could not calculate EWS for "
                  f"{variable.name}.")
    elif variable.spatial:
        if cfg.map_data:
            ews_calculations_main(variable, spatial_data_file_loading, spatial_ews_calculations, path=path,
                                  timer_on=timer_on)
        elif not cfg.map_data:
            print(f"Map data == False in configuration, could not calculate EWS for {variable.name}.")


def ews_calculations_main(variable, loading_function, calculation_function, path='./1/', timer_on=False):
    state_variable, files_present = loading_function(variable, path=path)

    if files_present:
        print(f"Started EWS calculations for {variable.name} in {path}")
        if timer_on:
            start_time = time.time()

        calculation_function(variable, state_variable, path=path)

        if timer_on:
            end_time = time.time()
            print(f"Elapsed time for EWS calculations for {variable.name} in {path} equals:", end_time - start_time,
                  '\n')


def EWS_weekly_calculations():
    start_time = time.time()
    for realization in range(1, cfg.nrOfSamples + 1):
        for variable in variables:
            ews_calculations_init(variable, path=f'./{realization}/', timer_on=True)
            if cfg.generate_dummy_datasets:
                generate_datasets(variable, path=f'./{realization}/', nr_realizations=cfg.nr_generated_datasets,
                                  detrending_temp='None', sigma=100, method1=cfg.method_1, method2=cfg.method_2,
                                  method3=cfg.method_3)
                ews_calculations_generated_datasets_init(variable, path=f'./{realization}/',
                                                         nr_realizations=cfg.nr_generated_datasets,
                                                         timer_on=True, method1=cfg.method_1, method2=cfg.method_2,
                                                         method3=cfg.method_3)
    end_time = time.time() - start_time
    print(f"Total elapsed time equals: {end_time} seconds")
    
    
def EWS_hourly_calculations():    
    for i in range(cfg.stepsTotal):
        fpath = str(i).zfill(2)
        for variable in variables:
            ews_calculations_init(variable, path=f'./h{fpath}/', timer_on=True)          


EWS_hourly_calculations()
