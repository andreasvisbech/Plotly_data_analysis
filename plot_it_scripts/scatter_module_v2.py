import lmfit
from sklearn.metrics import r2_score
#from scipy import stats
#import scipy

# Import functions from other scripts
from plot_it_scripts.plotting_script import *
from plot_it_scripts.main_functions import *

def scatter_plot_without_fit(master_dict, sample_idx, user_input_dict, plot_dict, xs, ys, param_dict):

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

    if error_id in master_dict['columns']:

        std_dev = plot_dict['errors']

        plot_func(
            figure, graph_name, xs, ys, std_dev, plot_marker, x_title, y_title, subplot_row,
            subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
            user_input_dict, master_dict)
        plot_func(
            plot_figure, graph_name, xs, ys, std_dev, plot_marker, x_title, y_title, subplot_row,
            subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
            user_input_dict, master_dict)


    else:
        plot_func(
            figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col, 'None',
            sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
        plot_func(
            plot_figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col, 'None',
            sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

    # Add the appropriate color to the table coloring list
    table_color_list_manager(color_list[color_count], plot_dict['table_color_list'])

def scatter_plot_with_fit(sample_idx, user_input_dict, plot_dict, param_dict, master_dict, xs, ys, unique_x):

    fit_intervals = user_input_dict['fit_intervals']
    fit_models = user_input_dict['fit_models']
    fit_modes = user_input_dict['fit_modes']

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
    sample_notes = user_input_dict['sample_notes']

    # Define the mode of fitting (local or global)
    fit_mode = fit_modes[sample_idx].lower()

    # Define the model
    if fit_models[sample_idx].lower() in ['fida_1to1', '1to1', 'one2one']:
        fitting_model = 'fida_1to1'
    elif fit_models[sample_idx].lower() in ['fida_excess', 'excess']:
        fitting_model = 'fida_excess'
    elif fit_models[sample_idx].lower() in ['4pl']:
        fitting_model = '4pl'
    elif fit_models[sample_idx].lower() in ['3pl', 'hill', 'hill_simple']:
        fitting_model = '3pl'
    else:
        fitting_model = None
        print('Fitting model is not reconngnized')

    y_mean, std_dev = replicate_mean_error(unique_x, xs, ys)

    plot_func(
        figure, graph_name, unique_x, y_mean, std_dev, plot_marker, x_title, y_title, subplot_row,
        subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
        user_input_dict, master_dict)
    plot_func(
        plot_figure, graph_name, unique_x, y_mean, std_dev, plot_marker, x_title, y_title, subplot_row,
        subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
        user_input_dict, master_dict)

    # Fit parameters based on experimental data
    fit_result = scatter_fit(xs, ys, fit_mode, fitting_model, param_dict, master_dict)

    # Calculate the fit by using the estimated parameters
    # Extract boundaries for the fit from the excel sheet.
    # # If no interval is supplied code will automatically assign 0;0.
    fitting_interval = fit_intervals[sample_idx].split(';')
    x_fit_plot, y_fit_plot = calculate_the_fit(xs, ys, fitting_interval, fitting_model, fit_result, fit_mode)

    # Plot the calculated fit
    plot_func(
        figure, graph_name+'_fit', x_fit_plot, y_fit_plot, 'None', 'line', x_title, y_title, subplot_row,
        subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
        user_input_dict, master_dict)
    plot_func(
        plot_figure, graph_name + '_fit', x_fit_plot, y_fit_plot, 'None', 'line', x_title, y_title, subplot_row,
        subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
        user_input_dict, master_dict)

    fit_statistics(xs, ys, fitting_model, fit_result, fit_mode, master_dict)

    fit_params_to_master_dict(sample_idx, master_dict, user_input_dict, plot_dict, fit_result, fitting_model, fit_mode)

    # Calculate residuals. If specified by user we also plot the residuals.
    residuals = calc_residuals(xs, ys, fitting_model, fit_result)

    if param_dict['scatter_residuals'] == 'yes':
        plot_func(figure, graph_name + '_residuals', xs, residuals, 'None', 'dots', x_title, y_title,
                  subplot_row, subplot_col, 'scatter_residuals', sample_idx, param_dict, color_list,
                  color_count,
                  user_input_dict, master_dict)

    # Add the appropriate color to the table coloring list
    table_color_list_manager(color_list[color_count], plot_dict['table_color_list'])

def scatter_fit(xs, ys, fit_mode, model, param_dict, master_dict):

    # Make sure values are arrays
    xs = np.array(xs)
    ys = np.array(ys)

    if fit_mode == 'local':

        # In local fitting mode the list of unique x values is the only x list needed.
        xs_unique = np.unique((xs))
        ys_mean, ys_std = replicate_mean_error(xs_unique, xs, ys)

        fit_result = local_fitting(model, xs_unique, ys_mean, param_dict)
        # print(fit_result.params.pretty_print())

    elif fit_mode == 'global':

        # If the data contains replicate values these should be split into separate data arrays
        xs, ys = split_replicates(xs, ys)

        # Save the number of replicates in the data
        master_dict['replicates'] = len(ys)

        fit_result = global_fitting(model, xs, ys, param_dict)
        # print(fit_result.params.pretty_print())

    else:
        fit_result = None

    return fit_result

def local_fitting(model_name, xs, ys, param_dict):

    # Defining the fitting parameters.
    fit_params = lmfit.Parameters()

    if model_name == 'fida_1to1':

        fit_params.add('RI_1', value=0.01, min=param_dict['RI min'], max=param_dict['RI max'])
        fit_params.add('RIA_1', value=0.01, min=param_dict['RIA min'], max=param_dict['RIA max'])
        fit_params.add('KD_1', value=0.01, min=param_dict['KD min'], max=param_dict['KD max'])
        fit_params.add('CI_1', value=0.01, min=param_dict['CI min'], max=param_dict['CI max'])

        out = lmfit.minimize(local_minimizer_fida1to1, fit_params, args=(xs, ys), method='least_squares')

    elif model_name == 'fida_excess':

        fit_params.add('RI_1', value=0.01, min=param_dict['RI min'], max=param_dict['RI max'])
        fit_params.add('RIA_1', value=0.01, min=param_dict['RIA min'], max=param_dict['RIA max'])
        fit_params.add('KD_1', value=0.01, min=param_dict['KD min'], max=param_dict['KD max'])
        fit_params.add('CI_1', value=0.01, min=param_dict['CI min'], max=param_dict['CI max'])

        out = lmfit.minimize(local_minimizer_fida_excess, fit_params, args=(xs, ys), method='least_squares')

    elif model_name == '4pl':

        fit_params.add('Bmin_1', value=0.01, min=param_dict['Bmin_min'], max=param_dict['Bmin_max'])
        fit_params.add('Bmax_1', value=0.01, min=param_dict['Bmax_min'], max=param_dict['Bmax_max'])
        fit_params.add('KD_1', value=0.01, min=param_dict['KD min'], max=param_dict['KD max'])
        fit_params.add('kcoop_1', value=0.01, min=param_dict['kcoop_min'], max=param_dict['kcoop_max'])

        out = lmfit.minimize(local_minimizer_4pl, fit_params, args=(xs, ys), method='least_squares')

    elif model_name == '3pl':

        fit_params.add('Bmin_1', value=0.01, min=param_dict['Bmin_min'], max=param_dict['Bmin_max'])
        fit_params.add('Bmax_1', value=0.01, min=param_dict['Bmax_min'], max=param_dict['Bmax_max'])
        fit_params.add('KD_1', value=0.01, min=param_dict['KD min'], max=param_dict['KD max'])

        out = lmfit.minimize(local_minimizer_3pl, fit_params, args=(xs, ys), method='least_squares')

    else:
        out = None

    return out

def global_fitting(model_name, xs, ys, param_dict):

    # Defining the fitting parameters.
    fit_params = lmfit.Parameters()

    if model_name == 'fida_1to1':

        for iy, y in enumerate(ys):
            # Loading in parameters for the fit
            fit_params.add(f'RI_{iy + 1}', value=0.01, min=param_dict['RI min'], max=param_dict['RI max'])
            fit_params.add(f'RIA_{iy + 1}', value=0.01, min=param_dict['RIA min'], max=param_dict['RIA max'])
            fit_params.add(f'KD_{iy + 1}', value=0.01, min=param_dict['KD min'], max=param_dict['KD max'])
            fit_params.add(f'CI_{iy + 1}', value=0.01, min=param_dict['CI min'], max=param_dict['CI max'])

        for iy in range(1, len(ys)):
            fit_params[f'KD_{iy + 1}'].expr = 'KD_1'
            fit_params[f'RI_{iy + 1}'].expr = 'RI_1'
            fit_params[f'RIA_{iy + 1}'].expr = 'RIA_1'
        out = lmfit.minimize(global_minimizer_fida1to1, fit_params, args=(xs, ys), method='least_squares')

    elif model_name == 'fida_excess':

        for iy, y in enumerate(ys):
            # Loading in parameters for the fit
            fit_params.add(f'RI_{iy + 1}', value=0.01, min=param_dict['RI min'], max=param_dict['RI max'])
            fit_params.add(f'RIA_{iy + 1}', value=0.01, min=param_dict['RIA min'], max=param_dict['RIA max'])
            fit_params.add(f'KD_{iy + 1}', value=0.01, min=param_dict['KD min'], max=param_dict['KD max'])
            fit_params.add(f'CI_{iy + 1}', value=0.01, min=param_dict['CI min'], max=param_dict['CI max'])

        for iy in range(1, len(ys)):
            fit_params[f'KD_{iy + 1}'].expr = 'KD_1'
            fit_params[f'RI_{iy + 1}'].expr = 'RI_1'
            fit_params[f'RIA_{iy + 1}'].expr = 'RIA_1'
            fit_params[f'CI_{iy + 1}'].expr = 'CI_1'

        out = lmfit.minimize(global_minimizer_fida_excess, fit_params, args=(xs, ys), method='least_squares')

    elif model_name == '4pl':

        for iy, y in enumerate(ys):
            # Loading in parameters for the fit
            fit_params.add(f'Bmin_{iy + 1}', value=0.01, min=param_dict['Bmin_min'], max=param_dict['Bmin_max'])
            fit_params.add(f'Bmax_{iy + 1}', value=0.01, min=param_dict['Bmax_min'], max=param_dict['Bmax_max'])
            fit_params.add(f'KD_{iy + 1}', value=0.01, min=param_dict['KD min'], max=param_dict['KD max'])
            fit_params.add(f'kcoop_{iy + 1}', value=0.01, min=param_dict['kcoop_min'], max=param_dict['kcoop_max'])

        for iy in range(1, len(ys)):
            fit_params[f'KD_{iy + 1}'].expr = 'KD_1'
            fit_params[f'Bmin_{iy + 1}'].expr = 'Bmin_1'
            fit_params[f'Bmax_{iy + 1}'].expr = 'Bmax_1'
            fit_params[f'kcoop_{iy + 1}'].expr = 'kcoop_1'

        out = lmfit.minimize(global_minimizer_4pl, fit_params, args=(xs, ys), method='least_squares')

    elif model_name == '3pl':

        for iy, y in enumerate(ys):
            # Loading in parameters for the fit
            fit_params.add(f'Bmin_{iy + 1}', value=0.01, min=param_dict['Bmin_min'], max=param_dict['Bmin_max'])
            fit_params.add(f'Bmax_{iy + 1}', value=0.01, min=param_dict['Bmax_min'], max=param_dict['Bmax_max'])
            fit_params.add(f'KD_{iy + 1}', value=0.01, min=param_dict['KD min'], max=param_dict['KD max'])
            fit_params.add(f'kcoop_{iy + 1}', value=0.01, min=param_dict['kcoop_min'], max=param_dict['kcoop_max'])

        for iy in range(1, len(ys)):
            fit_params[f'KD_{iy + 1}'].expr = 'KD_1'
            fit_params[f'Bmin_{iy + 1}'].expr = 'Bmin_1'
            fit_params[f'Bmax_{iy + 1}'].expr = 'Bmax_1'

        out = lmfit.minimize(global_minimizer_3pl, fit_params, args=(xs, ys), method='least_squares')

    else:
        out = None

    return out

def calculate_the_fit(xs, ys, fitting_interval, model, fit_result, fit_mode):

    for j in range(0, len(fitting_interval)):
        fitting_interval[j] = float(
            fitting_interval[j].replace(",", "."))  # Replace commas with proper dots to get proper numbers.
    fit_start = fitting_interval[0]
    fit_end = fitting_interval[1]
    x_fit_plot = interval_generator(fit_start, fit_end)

    if model == 'fida_1to1':
        RI = fit_result.params['RI_1'].value
        RIA = fit_result.params['RIA_1'].value
        KD = fit_result.params['KD_1'].value
        y_fit_plot = model_fida_1to1(x_fit_plot, RI, RIA, KD)

    elif model == 'fida_excess':
        RI = fit_result.params['RI_1'].value
        RIA = fit_result.params['RIA_1'].value
        KD = fit_result.params['KD_1'].value
        CI = fit_result.params['CI_1'].value
        y_fit_plot = model_fida_excess(x_fit_plot, RI, RIA, KD, CI)

    elif model == '4pl':
        Bmin = fit_result.params['Bmin_1'].value
        Bmax = fit_result.params['Bmax_1'].value
        KD = fit_result.params['KD_1'].value
        kcoop = fit_result.params['kcoop_1'].value
        y_fit_plot = model_4pl(x_fit_plot , Bmin, Bmax, KD, kcoop)

    elif model == '3pl':
        Bmin = fit_result.params['Bmin_1'].value
        Bmax = fit_result.params['Bmax_1'].value
        KD = fit_result.params['KD_1'].value
        y_fit_plot = model_3pl(x_fit_plot , Bmin, Bmax, KD)

    else:
        y_fit_plot = None

    return x_fit_plot, y_fit_plot

def fit_statistics(xs, ys, model, fit_result, fit_mode, master_dict):

    # Make sure data is array
    xs = np.array(xs)
    ys = np.array(ys)

    # Define the values to be used for fitting statistics calculations
    if fit_mode == 'local':
        # In local fitting mode the list of unique x values is the only x list needed.
        xs_unique = np.unique((xs))
        x_val = xs_unique
        y_val, y_val_std = replicate_mean_error(xs_unique, xs, ys)

    elif fit_mode == 'global':
        x_val = xs
        y_val = ys

    else:
        x_val = None
        y_val = None

    # Length of data
    n = len(y_val)

    if model == 'fida_1to1':
        RI = fit_result.params['RI_1'].value
        RIA = fit_result.params['RIA_1'].value
        KD = fit_result.params['KD_1'].value
        p = 3
        y_fit = model_fida_1to1(x_val, RI, RIA, KD)
        y_fit_full = model_fida_1to1(xs, RI, RIA, KD)

    elif model == 'fida_excess':
        RI = fit_result.params['RI_1'].value
        RIA = fit_result.params['RIA_1'].value
        KD = fit_result.params['KD_1'].value
        CI = fit_result.params['CI_1'].value
        p = 4
        y_fit = model_fida_excess(x_val, RI, RIA, KD, CI)
        y_fit_full = model_fida_excess(xs, RI, RIA, KD, CI)

    elif model == '4pl':
        Bmin = fit_result.params['Bmin_1'].value
        Bmax = fit_result.params['Bmax_1'].value
        KD = fit_result.params['KD_1'].value
        kcoop = fit_result.params['kcoop_1'].value
        p = 4
        y_fit = model_4pl(x_val, Bmin, Bmax, KD, kcoop)
        y_fit_full = model_4pl(xs, Bmin, Bmax, KD, kcoop)

    elif model == '3pl':
        Bmin = fit_result.params['Bmin_1'].value
        Bmax = fit_result.params['Bmax_1'].value
        KD = fit_result.params['KD_1'].value
        p = 3
        y_fit = model_3pl(x_val, Bmin, Bmax, KD)
        y_fit_full = model_3pl(xs, Bmin, Bmax, KD)

    else:
        p = None
        y_fit = None
        y_fit_full = None

    # Calculate R square based on the fitted values and the observed values
    r_square = r2_score(y_val, y_fit)

    # Calculate the adjusted R square
    term1 = 1 - r_square
    term2 = n - 1
    term3 = n - p - 1
    r_square_adj = 1 - (term1 * term2) / (term3)

    # Append R square values to master dict
    r_square = str(round(r_square, 3))
    r_square_adj = str(round(r_square_adj, 3))
    master_dict['R_square'].append('R<sup>2</sup>=' + r_square + '<br>' + 'R<sup>2</sup><sub>adj</sub>=' + r_square_adj)

    # Calculate RMSE.
    # The RMSE value is calculated on the full data i.e. not the means of replicates
    MSE = np.square(np.subtract(ys, y_fit_full)).mean()
    RMSE = np.sqrt(MSE)
    RMSE = str(round(RMSE, 3))
    #MSE_full = np.square(np.subtract(ys, y_fit_full)).mean()
    #RMSE_full = np.sqrt(MSE_full)
    #RMSE_full = str(round(RMSE_full, 3))

    # Append RMSE to master dict
    master_dict['RMSE'].append('RMSE=' + RMSE)

def calc_residuals(xs, ys, model, fit_result):

    if model == 'fida_1to1':
        RI = fit_result.params['RI_1'].value
        RIA = fit_result.params['RIA_1'].value
        KD = fit_result.params['KD_1'].value
        y_fit = model_fida_1to1(xs, RI, RIA, KD)
        residuals = ys - y_fit

    elif model == 'fida_excess':
        RI = fit_result.params['RI_1'].value
        RIA = fit_result.params['RIA_1'].value
        KD = fit_result.params['KD_1'].value
        CI = fit_result.params['CI_1'].value
        y_fit = model_fida_excess(xs, RI, RIA, KD, CI)
        residuals = ys - y_fit

    elif model == '4pl':
        Bmin = fit_result.params['Bmin_1'].value
        Bmax = fit_result.params['Bmax_1'].value
        KD = fit_result.params['KD_1'].value
        kcoop = fit_result.params['kcoop_1'].value
        y_fit = model_4pl(xs, Bmin, Bmax, KD, kcoop)
        residuals = ys - y_fit

    elif model == '3pl':
        Bmin = fit_result.params['Bmin_1'].value
        Bmax = fit_result.params['Bmax_1'].value
        KD = fit_result.params['KD_1'].value
        y_fit = model_3pl(xs, Bmin, Bmax, KD)
        residuals = ys - y_fit

    else:
        residuals = None

    return residuals

def fit_params_to_master_dict(sample_idx, master_dict, user_input_dict, plot_dict, fit_result, fitting_model, fit_mode):

    graph_name = plot_dict['graph_names'][sample_idx]
    sample_notes = user_input_dict['sample_notes'][sample_idx]

    master_dict['ID_list_new'].append(graph_name)
    master_dict['notes_list'].append(sample_notes)
    master_dict['model_list'].append(fitting_model + ', ' + fit_mode)

    KD = fit_result.params['KD_1'].value
    KD = str(round(KD, 3))
    master_dict['KD_fit'].append(KD)

    master_dict['chi_square'].append('test')

    if fitting_model == 'fida_1to1':
        RI = fit_result.params['RI_1'].value
        RI = str(round(RI, 3))
        RIA = fit_result.params['RIA_1'].value
        RIA = str(round(RIA, 3))
        CI = fit_result.params['CI_1'].value
        CI = str(round(CI, 3))
        master_dict['fit_parameters'].append('RI=' + RI + ', RIA=' + RIA)

    elif fitting_model == 'fida_excess':
        RI = fit_result.params['RI_1'].value
        RI = str(round(RI, 3))
        RIA = fit_result.params['RIA_1'].value
        RIA = str(round(RIA, 3))
        CI = fit_result.params['CI_1'].value
        CI = str(round(CI, 3))
        master_dict['fit_parameters'].append('RI=' + RI + ', RIA=' + RIA + ', CI=' + CI)

    elif fitting_model == '4pl':
        Bmin = fit_result.params['Bmin_1'].value
        Bmin = str(round(Bmin, 3))
        Bmax = fit_result.params['Bmax_1'].value
        Bmax = str(round(Bmax, 3))
        kcoop = fit_result.params['kcoop_1'].value
        kcoop = str(round(kcoop, 3))
        master_dict['fit_parameters'].append('Bmin=' + Bmin + ', Bmax=' + Bmax + ', kcoop=' + kcoop)

    elif fitting_model == '3pl':
        Bmin = fit_result.params['Bmin_1'].value
        Bmin = str(round(Bmin, 3))
        Bmax = fit_result.params['Bmax_1'].value
        Bmax = str(round(Bmax, 3))
        master_dict['fit_parameters'].append('Bmin=' + Bmin + ', Bmax=' + Bmax)

def local_minimizer_fida1to1(pars, xs, ys):

    params = pars.valuesdict()

    RI = params['RI_1']
    RIA = params['RIA_1']
    KD = params['KD_1']

    resid = ys - model_fida_1to1(xs, RI, RIA, KD)

    return resid

def local_minimizer_fida_excess(pars, xs, ys):

    params = pars.valuesdict()

    RI = params['RI_1']
    RIA = params['RIA_1']
    KD = params['KD_1']
    CI = params['CI_1']

    resid = ys - model_fida_excess(xs, RI, RIA, KD, CI)

    return resid

def local_minimizer_4pl(pars, xs, ys):

    params = pars.valuesdict()

    Bmin = params['Bmin_1']
    Bmax = params['Bmax_1']
    KD = params['KD_1']
    kcoop = params['kcoop_1']

    resid = ys - model_4pl(xs, Bmin, Bmax, KD, kcoop)

    return resid

def local_minimizer_3pl(pars, xs, ys):

    params = pars.valuesdict()

    Bmin = params['Bmin_1']
    Bmax = params['Bmax_1']
    KD = params['KD_1']

    resid = ys - model_3pl(xs, Bmin, Bmax, KD)

    return resid

def global_minimizer_fida1to1(params, xs, ys):

    resid_list = []

    for i, y in enumerate(ys):

        RI = params['RI_%i' % (i + 1)].value
        RIA = params['RIA_%i' % (i + 1)].value
        KD = params['KD_%i' % (i + 1)].value

        resid = y - model_fida_1to1(xs[i], RI, RIA, KD)
        resid = resid.tolist()

        resid_list = resid_list + resid

    return np.array(resid_list)

def global_minimizer_fida_excess(params, xs, ys):

    resid_list = []

    for i, y in enumerate(ys):
        RI = params['RI_%i' % (i + 1)].value
        RIA = params['RIA_%i' % (i + 1)].value
        KD = params['KD_%i' % (i + 1)].value
        CI = params['CI_%i' % (i + 1)].value

        resid = y - model_fida_excess(xs[i], RI, RIA, KD, CI)
        resid = resid.tolist()

        resid_list = resid_list + resid

    return np.array(resid_list)

def global_minimizer_4pl(params, xs, ys):
    resid_list = []

    for i, y in enumerate(ys):
        Bmin = params['Bmin_%i' % (i + 1)].value
        Bmax = params['Bmax_%i' % (i + 1)].value
        KD = params['KD_%i' % (i + 1)].value
        kcoop = params['kcoop_%i' % (i + 1)].value

        resid = y - model_4pl(xs[i], Bmin, Bmax, KD, kcoop)
        resid = resid.tolist()

        resid_list = resid_list + resid

    return np.array(resid_list)

def global_minimizer_3pl(params, xs, ys):
    resid_list = []

    for i, y in enumerate(ys):
        Bmin = params['Bmin_%i' % (i + 1)].value
        Bmax = params['Bmax_%i' % (i + 1)].value
        KD = params['KD_%i' % (i + 1)].value

        resid = y - model_3pl(xs[i], Bmin, Bmax, KD)
        resid = resid.tolist()

        resid_list = resid_list + resid

    return np.array(resid_list)

def fida_text2floats(x_id, y_id, df):
    func_list = []

    # Slice the dataframe to only include the rows that contain a concentration
    data_slice = df[df[x_id].str.contains('nM]', na=False)]
    xdata_list = data_slice[x_id].astype(str).str.replace(",", ".").tolist()
    ydata_list = pd.to_numeric(data_slice[y_id].astype(str).str.replace(",", ".")).tolist()
    for m in range(len(xdata_list)):
        conc = xdata_list[m][xdata_list[m].index('[') + 1: xdata_list[m].index(' nM]')]
        func_list.append(float(conc))

    return func_list, ydata_list

def scatter_data_slice(df, sample_idx, x_id, y_id, user_input_dict):

    data_interval = user_input_dict['data_interval']

    # Extracting the raw x and y values from the excel sheet.
    # Replace commas with dots in case user has e.g. Danish Excel.
    xs = df[x_id][pd.to_numeric(df[x_id], errors='coerce').notnull()]
    ys = df[y_id][pd.to_numeric(df[y_id], errors='coerce').notnull()]
    xs = pd.to_numeric(xs.astype(str).str.replace(",", "."))
    ys = pd.to_numeric(ys.astype(str).str.replace(",", "."))

    # Convert to arrays
    xs = xs.to_numpy()
    ys = ys.to_numpy()

    # Slicing the data so only data within the specified data interval is included.
    if data_interval[sample_idx] != 0:

        interval_var = list(data_interval[sample_idx].split(';'))
        interval_min, interval_max = float(interval_var[0]) , float(interval_var[1])

        xs_slice = xs[(xs >= interval_min) & (xs <= interval_max)]
        ys_slice = ys[(xs >= interval_min) & (xs <= interval_max)]

    else:

        xs_slice = xs
        ys_slice = ys

    # Create a list with unique concentration values. These should include ALL relevant x values.
    unique_x = np.unique(xs)
    #unique_x = list(dict.fromkeys(xs))

    return xs_slice, ys_slice, unique_x

def replicate_mean_error(unique_x, redundant_x, y_val):

    # Make sure data values are arrays
    unique_x = np.array(unique_x)
    redundant_x = np.array(redundant_x)
    y_val = np.array(y_val)

    # Allocate memory for mean values and standard deviations
    mean_list = np.zeros(len(unique_x))
    std_list = np.zeros(len(unique_x))

    for a, x_val in enumerate(unique_x):

        # Isolate the y values on all the indexes where x values are identical.
        ys = y_val[redundant_x==x_val]

        # Calculuate mean and standard deviation
        ys_mean = np.mean(ys)
        ys_std = np.std(ys)

        # Add values to
        mean_list[a] = ys_mean
        std_list[a] = ys_std

    return mean_list, std_list

def interval_generator(start, end):
    # The function is used to generate values of an interval with varying step size.
    # This is to have a nice distribution of data points throughout the interval without having
    # excessive number of points which make the html files large and laggy.

    counter = start
    interval = [counter]

    if counter == 0:
        counter = counter + 0.00001

    while counter < end:
        counter = counter + counter * 0.005  # Setting the step size increment
        interval.append(counter)

    return np.array(interval)

def split_replicates(xs, ys):

    x_out = []
    y_out = []

    xs_unique = np.unique((xs))

    count_list = []
    dummy_list = []
    for i, x in enumerate(xs):

        count = dummy_list.count(x)
        count_list.append(count)

        dummy_list.append(x)

    count_list = np.array(count_list)

    for j in range(max(count_list) + 1):
        x_new = xs[count_list == j]
        x_out.append(x_new)

        y_new = ys[count_list == j]
        y_out.append(y_new)

    return x_out , y_out

def model_fida_1to1(x, ri, ria, kd):
    y = (1 + (1 / kd) * x) / (((1 / ri) - (1 / ria)) + (1 + ((1 / kd) * x)) * (1 / ria))
    return y

def model_fida_excess(x, ri, ria, kd, ci):
    y = 1 / ((ci + x + kd - np.sqrt((ci + x + kd) ** 2 - 4 * ci * x)) / (2 * ria * ci) + (
            ci - x - kd + np.sqrt((ci + x + kd) ** 2 - 4 * ci * x)) / (2 * ri * ci))
    return y

def model_4pl(x, bmin, bmax, kd, k_coop):
    y = bmax + (bmin - bmax) / (1 + (x / kd) ** k_coop)
    return y

def model_3pl(x, bmin, bmax, kd):
    y = bmin + (x * (bmax - bmin)) / (kd + x)
    return y