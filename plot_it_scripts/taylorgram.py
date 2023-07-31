# Import functions from other scripts
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from lmfit import Model
import peakutils

from plot_it_scripts.plotting_script import *
from plot_it_scripts.scatter_module_v2 import *

def taylor_main(df, master_dict, sample_idx, user_input_dict, plot_dict, xs, ys, param_dict):

    # Appending dummy values to the master dict to avoid mispairing of cells in the output table
    master_dict['ID_list_new'].append(user_input_dict['ID_list'][sample_idx])
    master_dict['notes_list'].append(user_input_dict['sample_notes'][sample_idx])
    master_dict['model_list'].append(' ')
    master_dict['fit_parameters'].append(' ')
    master_dict['R_square'].append(' ')
    master_dict['KD_fit'].append(' ')
    master_dict['RMSE'].append(' ')
    master_dict['chi_square'].append(' ')

    # Loading local variables for plotting
    figure = plot_dict['figure']
    plot_figure = plot_dict['plot_figure']
    graph_name = plot_dict['graph_names'][sample_idx]
    plot_marker = user_input_dict['plot_markers'][sample_idx]
    x_title = user_input_dict['x_titles'][sample_idx]
    y_title = user_input_dict['y_titles'][sample_idx]
    subplot_row = user_input_dict['subplot_row'][sample_idx]
    subplot_col = user_input_dict['subplot_col'][sample_idx]
    color_list = plot_dict['color_list']
    color_count = plot_dict['color_count']
    error_id = plot_dict['error_id']

    # Make sure data is numpy arrays
    xs = np.array(xs)
    ys = np.array(ys)

    # Find the tR value
    tR = find_tR(xs, ys)
    master_dict['tR'].append(round(tR, 2))

    # Normalize the x values
    xs_norm = xs - tR

    plot_func(
        figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col, 'None',
        sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
    plot_func(
        plot_figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col, 'None',
        sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

    f = interp1d(xs, ys)

    xs_new = np.linspace(list(xs)[0], list(xs)[-1], len(xs)*10)
    ys_new = f(list(xs_new))

    if param_dict['peak_analysis'] == 'yes':
        xs_tda, ys_tda = fit_Rh(xs, ys, param_dict, master_dict)

        plot_func(
            figure, graph_name + '_tda', xs_tda, ys_tda, 'None', plot_marker, x_title, y_title,
            subplot_row, subplot_col, 'AKTA_baseline', sample_idx, param_dict, color_list, color_count,
            user_input_dict, master_dict)
    else:
        master_dict['Rh'].append('N/A')

    if param_dict['peak_fitting'] == 'yes':
        xs_fit, ys_fit = fit_gauss(xs_norm, ys)

        plot_func(
            figure, graph_name + '_fit', xs_fit, ys_fit, 'None', plot_marker, x_title, y_title,
            subplot_row, subplot_col, 'AKTA_baseline', sample_idx, param_dict, color_list, color_count,
            user_input_dict, master_dict)
        plot_func(
            plot_figure, graph_name + '_fit', xs_fit, ys_fit, 'None', plot_marker, x_title, y_title,
            subplot_row, subplot_col, 'AKTA_baseline', sample_idx, param_dict, color_list, color_count,
            user_input_dict, master_dict)

    # Add the appropriate color to the table coloring list
    table_color_list_manager(user_input_dict['plot_color'], plot_dict['table_color_list'])

def find_tR(xs, ys):

    # Smooth the data to ensure spikes is not seen as max values
    ys_smooth = savgol_filter(ys, 31, 2)

    y_max = max(ys_smooth)
    tR = xs[ys_smooth == y_max][0]

    return tR

def fit_gauss(xs, ys):

    gmodel = Model(gaussian)

    result = gmodel.fit(ys-3.35, x=xs, cen=0, amp=5, wid=0.5)

    #print(result.fit_report())

    return xs, result.best_fit+3.35

def fit_Rh(xs, ys, param_dict, master_dict):

    if param_dict['time_unit'] == 'min':
        conv = 60
        xs_out = xs
    else:
        conv = 1

    # Pressure in Pa
    pres = float(param_dict['pressure'])*100

    # capillary radius in m
    rc = float(param_dict['cap_diameter'])/10**6/2

    # Viscocity
    visc = float(param_dict['visc'])

    # Capillary length in m
    cap_len = float(param_dict['cap_length'])

    # Convert x values to appropriate time unit
    xs = xs*conv

    # Reference tine
    tR = find_tR(xs, ys)

    model1 = Model(mod_gaussian)
    result1 = model1.fit(ys, x=xs, base=0, cen=tR, amp=5, wid=0.5)

    # Set sigma from the gaussian fit
    sigma = result1.params['wid'].value

    # Mu is the volumetric flow rate
    mu = pres*np.pi*(rc**4)/8/visc/cap_len

    k = ((mu**2)*(sigma**2))/(2*tR)

    D = ((rc**2)*(mu**2))/(48*k)

    kB = 1.380649*(10**(-23))

    T = 298.15

    # Hydrodynamic radius in m
    Rh = (kB*T)/(6*np.pi*visc*D)

    # Convert from m to nm
    Rh = Rh*10**9
    Rh = str(round(Rh, 2))

    master_dict['Rh'].append(Rh)

    return xs_out, result1.best_fit

def polydis_eval(xs, ys, tR):

    tR = tR - 0.02

    xs_new = xs[(xs <= tR) & (ys > 3.43)]
    ys_new = ys[(xs <= tR) & (ys > 3.43)]

    xs_mod = (xs_new - tR)**2
    ys_mod = np.log(ys_new)

    return xs_mod, ys_mod

def gaussian(x, amp, cen, wid):
    return amp * np.exp(-((x - cen)**2) / (2*(wid**2)))

def mod_gaussian(x, base, amp, cen, wid):
    return base+(amp * np.exp(-((x - cen)**2) / (2*(wid**2))))
