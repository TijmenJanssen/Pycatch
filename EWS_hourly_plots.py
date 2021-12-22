import numpy as np
import os

import EWS_main_configuration as cfg
import EWSPy as ews
import EWS_StateVariables as ews_sv

import matplotlib.pyplot as plt

## State variables for EWS ##
# State variables present in EWS_StateVariables can be added through EWS_main_configuration.py
# variables = ews_sv.variables_weekly
variables = ews_sv.variables_hourly

names = []
for variable in variables:
    names.append([f'{variable.full_name} as {variable.name}'])

## Number of timesteps over which EWS can be calculated ##
# This number can be different for the weekly and hourly model
# number_of_timesteps = cfg.number_of_timesteps_weekly
number_of_timesteps = cfg.number_of_timesteps_hourly


## Statistical EWS ##
ews_temporal_signals = {'t.mn': "mean", 't.std': "standard deviation", 't.var': "variance",
                        't.cv': "coefficient of variation", 't.skw': "skewness", 't.krt': "kurtosis",
                        't.dfa': "detrended fluctuation analysis", 't.acr': "autocorrelation", 't.rr': "return rate",
                        't.coh': "conditional heteroskedasticity", 'timeseries': "timeseries"}
ews_spatial_signals = {'s.mn': "mean", 's.std': "standard deviation", 's.var': "variance", 's.skw': "skewness",
                       's.krt': "kurtosis", 's.mI': "Moran's I"}


def plot2(variable1, signal1='None', variable2='None', signal2='None', path='./1/', save=False, show=False):
    if variable1.spatial:
        x_axis1 = np.arange(cfg.interval_map_snapshots, number_of_timesteps + cfg.interval_map_snapshots,
                            cfg.interval_map_snapshots)
    if variable1.temporal:
        x_axis1 = np.arange(0, number_of_timesteps, variable1.window_size - variable1.window_overlap)

    if signal1 == 'timeseries':
        fname = ews.file_name_str(variable1.name, number_of_timesteps)
        fpath = os.path.join(path + fname)
        timeseries_y_axis = np.loadtxt(fpath + '.numpy.txt')
        timeseries_x_axis = np.arange(0, number_of_timesteps, 1)
        plt.plot(timeseries_x_axis, timeseries_y_axis, label=f'Continues measurement of {variable1.full_name}')
    elif signal1 != 'None':
        fpath = os.path.join(path + variable1.name + '.' + signal1)
        signal1_array = np.loadtxt(fpath + '.numpy.txt')
        if variable1.spatial:
            plt.plot(x_axis1, signal1_array, label=f'{variable1.full_name} {ews_spatial_signals[signal1]}')
        if variable1.temporal:
            plt.plot(x_axis1, signal1_array, label=f'{variable1.full_name} {ews_temporal_signals[signal1]}')

    if variable2 != 'None':
        if variable2.spatial:
            x_axis2 = np.arange(cfg.interval_map_snapshots, number_of_timesteps + cfg.interval_map_snapshots,
                                cfg.interval_map_snapshots)
        if variable2.temporal:
            x_axis2 = np.arange(0, number_of_timesteps, variable2.window_size - variable2.window_overlap)

        if signal2 == 'timeseries':
            fname = ews.file_name_str(variable2.name, number_of_timesteps)
            fpath = os.path.join(path + fname)
            timeseries_y_axis = np.loadtxt(fpath + '.numpy.txt')
            timeseries_x_axis = np.arange(0, number_of_timesteps, 1)
            plt.plot(timeseries_x_axis, timeseries_y_axis, label=f'Continues measurement of {variable2.full_name}')
        elif signal2 != 'None':
            fpath = os.path.join(path + variable2.name + '.' + signal2)
            signal2_array = np.loadtxt(fpath + '.numpy.txt')
            if variable2.spatial:
                plt.plot(x_axis2, signal2_array, label=f'{variable2.full_name} {ews_spatial_signals[signal2]}')
            if variable2.temporal:
                plt.plot(x_axis2, signal2_array, label=f'{variable2.full_name} {ews_temporal_signals[signal2]}')

    if variable2 != 'None':
        plt.xlabel('x - axis')
        plt.ylabel('y - axis')
        if variable1.temporal:
            if variable2.temporal:
                plt.title(f"{variable1.full_name} {ews_temporal_signals[signal1]} and {variable2.full_name} "
                          f"{ews_temporal_signals[signal2]}")
            if variable2.spatial:
                plt.title(f"{variable1.full_name} {ews_temporal_signals[signal1]} and {variable2.full_name} "
                          f"{ews_spatial_signals[signal2]}")
        if variable1.spatial:
            if variable2.temporal:
                plt.title(f"{variable1.full_name} {ews_spatial_signals[signal1]} and {variable2.full_name} "
                          f"{ews_temporal_signals[signal2]}")
            if variable2.spatial:
                plt.title(f"{variable1.full_name} {ews_spatial_signals[signal1]} and {variable2.full_name} "
                          f"{ews_spatial_signals[signal2]}")

    else:
        plt.xlabel('x - axis')
        plt.ylabel('y - axis')
        if variable1.spatial:
            plt.title(f"{variable1.full_name} {ews_spatial_signals[signal1]}")
        if variable1.temporal:
            plt.title(f"{variable1.full_name} {ews_temporal_signals[signal1]}")
    plt.legend()

    if save:
        if variable2 != 'None':
            plt.savefig(path + f"{variable1.name}_{signal1}_and_{variable2.name}_{signal2}.pdf", format="pdf")
        else:
            plt.savefig(path + f"{variable1.name}_{signal1}.pdf", format="pdf")
    if show:
        plt.show()


def user_plotmaker(path='./1/'):
    print("Variables present in the current run are:", names)
    print("Enter the short name for variable 1:")
    variable1_input = input()
    variable1 = [variable for variable in variables if variable.name == variable1_input][0]
    if variable1.temporal:
        print("EW signals present are:", ews_temporal_signals)
    if variable1.spatial:
        print("EW signals present are:", ews_spatial_signals)
    print("Enter the short name for the signal for variable 1:")
    signal1_input = input()

    print("Include a second state variable? [Y/n]")
    second_variable_input = input()
    if second_variable_input == 'Y' or second_variable_input == 'y':
        print("Enter the short name for variable 2:")
        variable2_input = input()

        variable2 = [variable for variable in variables if variable.name == variable2_input][0]
        if variable2.temporal:
            print("EW signals present are:", ews_temporal_signals)
        if variable2.spatial:
            print("EW signals present are:", ews_spatial_signals)
        print("Enter the short name for the signal for variable 1:")
        signal2_input = input()
    else:
        variable2 = 'None'
        signal2_input = 'None'

    print("Save the plot as a .pdf? [Y/n]")
    save_plot = input()
    if save_plot == 'Y' or save_plot == 'y':
        save = True
    else:
        save = False

    print("Show the plot when finished? [Y/n]")
    print("Note that the program is still running if the plot stays open.")
    show_plot = input()
    if show_plot == 'Y' or show_plot == 'y':
        show = True
    else:
        show = False

    plot2(variable1=variable1, signal1=signal1_input, variable2=variable2, signal2=signal2_input, path=path,
          save=save, show=show)


# TIJMEN: Deze path maakt niet uit

def user_plotmaker_looper(path='./1/'):  
    user_plotmaker(path=path)
    print("Would you like to make another plot in this directory? [Y/n]")
    answer = input()
    if answer == 'Y' or answer == 'y':
        user_plotmaker_looper(path=path)
    if answer == 'N' or answer == 'n':
        print("next file.")


#Hier map invoeren voor die je geplot wil hebben.
user_plotmaker_looper(path='./h02/')

# for i in range(cfg.stepsTotal):
#     fpath = str(i).zfill(2)
#     user_plotmaker_looper(path=f'./h{fpath}/')
