# Import relevant modules
import statistics

# import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import scipy.stats

# Import functions from other scripts
from plot_it_scripts.plotting_script import *
from plot_it_scripts.main_functions import *


def linear_model(x, a, b):
	y = x * a + b
	return y


def hill_equation(x, bmin, bmax, kd, k_coop):
	y = bmin + ((x ** k_coop) * (bmax - bmin)) / ((kd ** k_coop) + (x ** k_coop))
	return y


def hill_simple(x, bmin, bmax, kd):
	y = bmin + (x * (bmax - bmin)) / (kd + x)
	return y


def fit_4pl(x, bmin, bmax, kd, k_coop):
	y = bmax + (bmin - bmax) / (1 + (x / kd) ** k_coop)
	return y


def model_fida_1to1(x, ri, ria, kd):
	y = (1 + (1 / kd) * x) / (((1 / ri) - (1 / ria)) + (1 + ((1 / kd) * x)) * (1 / ria))
	return y


def model_fida_excess(x, ri, ria, kd, ci):
	y = 1 / ((ci + x + kd - np.sqrt((ci + x + kd) ** 2 - 4 * ci * x)) / (2 * ria * ci) + (
			ci - x - kd + np.sqrt((ci + x + kd) ** 2 - 4 * ci * x)) / (2 * ri * ci))
	return y


def scatter_data_slice(df, sample_idx, x_id, y_id, user_input_dict):
	data_interval = user_input_dict['data_interval']

	# Extracting the raw x and y values from the excel sheet.
	# Replace commas with dots in case user has e.g. Danish Excel.
	xs = df[x_id][pd.to_numeric(df[x_id], errors='coerce').notnull()]
	ys = df[y_id][pd.to_numeric(df[y_id], errors='coerce').notnull()]
	xs = pd.to_numeric(xs.astype(str).str.replace(",", "."))
	ys = pd.to_numeric(ys.astype(str).str.replace(",", "."))

	# Create a list with unique concentration values. These should include ALL relevant x values.
	unique_x = list(dict.fromkeys(xs))

	# Slicing the data so only data within the specified data interval is included.
	if data_interval[sample_idx] != 0:
		interval_var = list(data_interval[sample_idx].split(';'))
		xsys_interval = pd.concat([xs, ys], axis=1)
		xsys_interval_slice = xsys_interval[
			(xsys_interval[x_id] >= float(interval_var[0])) & (xsys_interval[x_id] <= float(interval_var[1]))]
		xs = xsys_interval_slice[x_id]
		ys = xsys_interval_slice[y_id]

	return xs, ys, unique_x


def scatter_plot_with_fit(sample_idx, user_input_dict, plot_dict, param_dict, master_dict, xs, ys, unique_x):
	fit_intervals = user_input_dict['fit_intervals']
	fit_models = user_input_dict['fit_models']
	fit_modes = user_input_dict['fit_modes']

	# Loading local variables for analysis
	# ID_list = user_input_dict['ID_list']
	# sample_notes = user_input_dict['sample_notes']

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

	# If local fitting approach the data will be averaged and plotted with error bars.
	# Else the full data set will be plotted.
	if fit_modes[sample_idx] == 'Local':
		y_mean, std_dev = replicate_mean_error(unique_x, xs, ys)

		plot_func(
			figure, graph_name, unique_x, y_mean, std_dev, plot_marker, x_title, y_title, subplot_row,
			subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
			user_input_dict, master_dict)
		plot_func(
			plot_figure, graph_name, unique_x, y_mean, std_dev, plot_marker, x_title, y_title, subplot_row,
			subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
			user_input_dict, master_dict)

	else:
		plot_func(
			figure, graph_name, xs, ys, np.zeros(len(xs)), plot_marker, x_title, y_title, subplot_row,
			subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
			user_input_dict, master_dict)
		plot_func(
			plot_figure, graph_name, xs, ys, np.zeros(len(xs)), plot_marker, x_title, y_title, subplot_row,
			subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
			user_input_dict, master_dict)

	# Extract boundaries for the fit from the excel sheet and turn to list of floats.
	# If no interval is supplied code will automaticallt assign 0;0.
	fitting_interval = fit_intervals[sample_idx].split(';')
	for j in range(0, len(fitting_interval)):
		fitting_interval[j] = float(
			fitting_interval[j].replace(",", "."))  # Replace commas with proper dots to get proper numbers.

	# Running fitting function and extracting parameters.
	if fit_models[sample_idx] == 'Hill':
		scatter_fit(
			hill_equation, 'Hill', xs, ys, fit_modes[sample_idx], fitting_interval[0], fitting_interval[1],
			sample_idx, param_dict, master_dict, user_input_dict, plot_dict)

	elif fit_models[sample_idx].lower() == 'hill_simple':
		scatter_fit(
			hill_simple, 'Hill_simple', xs, ys, fit_modes[sample_idx], fitting_interval[0], fitting_interval[1],
			sample_idx, param_dict, master_dict, user_input_dict, plot_dict)

	elif fit_models[sample_idx].lower() == '4pl':
		scatter_fit(
			fit_4pl, '4PL', xs, ys, fit_modes[sample_idx], fitting_interval[0], fitting_interval[1], sample_idx,
			param_dict, master_dict, user_input_dict, plot_dict)

	elif fit_models[sample_idx].lower() in ['fida_1to1', '1to1', 'one2one']:
		scatter_fit(
			model_fida_1to1, 'FIDA_1to1', xs, ys, fit_modes[sample_idx], fitting_interval[0],
			fitting_interval[1], sample_idx, param_dict, master_dict, user_input_dict, plot_dict)

	elif fit_models[sample_idx].lower() in ['fida_excess', 'excess']:
		scatter_fit(
			model_fida_excess, 'FIDA_excess', xs, ys, fit_modes[sample_idx], fitting_interval[0],
			fitting_interval[1], sample_idx, param_dict, master_dict, user_input_dict, plot_dict)


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

	plot_func(
		figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col, 'None',
		sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
	plot_func(
		plot_figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col, 'None',
		sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)


def replicate_mean_error(unique_x, redundant_x, y_val):
	"""
	The function takes a list of unique concentrations, the redundant concentration list and a list of y values
	corresponding to the redundant concentration list.
	The function then calculates appropriate mean values and the standard deviation on the values used for the
	mean calculation.

	:param unique_x:
	:param redundant_x:
	:param y_val:
	:return:
	"""

	func_mean = []
	func_std = []

	y_val = list(y_val)

	for a in range(len(unique_x)):
		# Getting indices of the specific concentration
		conc_indices = [k for k, x in enumerate(redundant_x) if x == unique_x[a]]

		func_var = []
		for b in range(len(conc_indices)):
			func_var.append(y_val[conc_indices[b]])
		func_mean.append(sum(func_var) / len(func_var))

		# If more than one replicate is present in the dataset calculate standard deviation.
		if len(conc_indices) == 1:
			func_std.append(0)
		elif len(conc_indices) > 1:
			func_std.append(statistics.stdev(func_var))

	return func_mean, func_std


def scatter_fit(
		model_function, model_name, x_val, y_val, fitting_mode, fitting_min, fitting_max, sample_idx,
		param_dict, master_dict, user_input_dict, plot_dict):
	# Loading local variables for plotting
	figure = plot_dict['figure']
	graph_name = plot_dict['graph_names'][sample_idx]
	x_title = user_input_dict['x_titles'][sample_idx]
	y_title = user_input_dict['y_titles'][sample_idx]
	subplot_row = user_input_dict['subplot_row'][sample_idx]
	subplot_col = user_input_dict['subplot_col'][sample_idx]
	color_list = plot_dict['color_list']
	color_count = plot_dict['color_count']

	# Loading local variables for analysis
	ID_list = user_input_dict['ID_list']
	sample_notes = user_input_dict['sample_notes']
	fit_models = user_input_dict['fit_models']

	x_val_list = []
	y_val_list = []

	unique_x = list(dict.fromkeys(x_val))

	alpha_level = param_dict['alpha_level']

	# Setting fitting boundaries depending on which model is used
	if model_name == 'Hill':
		bound_param = [
			(param_dict['Bmin_min'], param_dict['Bmax_min'], param_dict['KD_fit_min'], param_dict['k_coop_min']),
			(param_dict['Bmin_max'], param_dict['Bmax_max'], param_dict['KD_fit_max'], param_dict['k_coop_max'])]

	elif model_name == 'Hill_simple':
		bound_param = [
			(param_dict['Bmin_min'], param_dict['Bmax_min'], param_dict['KD_fit_min']),
			(param_dict['Bmin_max'], param_dict['Bmax_max'], param_dict['KD_fit_max'])
		]

	elif model_name == '4PL':
		bound_param = [
			(param_dict['Bmin_min'], param_dict['Bmax_min'], param_dict['KD_fit_min'], param_dict['k_coop_min']),
			(param_dict['Bmin_max'], param_dict['Bmax_max'], param_dict['KD_fit_max'], param_dict['k_coop_max'])
		]

	elif model_name == 'FIDA_1to1':
		bound_param = [
			(param_dict['RI min'], param_dict['RIA min'], param_dict['KD min']),
			(param_dict['RI max'], param_dict['RIA max'], param_dict['KD max'])
		]

	elif model_name == 'FIDA_excess':
		bound_param = (
			(param_dict['RI min'], param_dict['RIA min'], param_dict['KD min'], param_dict['CI min']),
			(param_dict['RI max'], param_dict['RIA max'], param_dict['KD max'], param_dict['CI max'])
		)

	else:
		bound_param = ()
		ValueError('Model is not known!?')

	if fitting_mode == 'Local':

		# In local fitting mode the list of unique x values is the only x list needed.
		x_val_list.append(unique_x)

		# Calculate mean and std error on replicate values in the data.
		# Replace zero values in standard errors for use in curve fitting
		y_mean = replicate_mean_error(unique_x, x_val, y_val)[0]
		y_val_list.append(y_mean)

	elif fitting_mode == 'Global':

		# Get the max number of times a concentration value occurs in raw data. This is the number of graphs to calculate.
		func_conc_count = []
		for m in unique_x:
			func_conc_count.append(list(x_val).count(m))
		graph_number = max(func_conc_count)

		# Create list of empty lists for filling with the individual graphs.
		x_val_list = [[] for _ in range(graph_number)]
		y_val_list = [[] for _ in range(graph_number)]

		# Sort the data values into the appropriate lists
		for j in unique_x:
			# Getting indices of the specific concentration
			conc_indices = [k for k, x in enumerate(list(x_val)) if x == j]
			for k in range(len(conc_indices)):
				x_val_list[k].append(list(x_val)[conc_indices[k]])
				y_val_list[k].append(list(y_val)[conc_indices[k]])

	# Define a list with sigma values in case the user wants to weight the curve fit based on errors.
	sigma_list = sigma_list_for_scatter_fit(unique_x, x_val, y_val, x_val_list, fitting_mode)

	# Emptying the lists before running the fitting
	master_dict['func_KD_list'] = []
	master_dict['func_R2_list'] = []

	func_KD_list = []
	func_R2_list = []

	for c in range(len(x_val_list)):

		# The fitting process is based on whether or not the user has requested weighting of the fit
		# using sigma values.
		if param_dict['scatter_fit_error_weighing'] == 'yes, as relative sigma':
			parameters, covariance = curve_fit(
				model_function, list(x_val_list[c]), list(y_val_list[c]),
				sigma=sigma_list[c], bounds=bound_param, method='trf',
				maxfev=10000)
		elif param_dict['scatter_fit_error_weighing'] == 'yes, as absolute sigma':
			parameters, covariance = curve_fit(
				model_function, list(x_val_list[c]), list(y_val_list[c]),
				sigma=sigma_list[c], absolute_sigma=True, bounds=bound_param, method='trf',
				maxfev=10000)
		else:
			parameters, covariance = curve_fit(
				model_function, list(x_val_list[c]), list(y_val_list[c]),
				bounds=bound_param, method='trf',
				maxfev=10000)

		if len(x_val_list) == 1:
			master_dict['ID_list_new'].append(str(ID_list[sample_idx]))
			master_dict['notes_list'].append(sample_notes[sample_idx])
			master_dict['model_list'].append(fit_models[sample_idx] + ', ' + fitting_mode)

		elif len(x_val_list) > 1:
			master_dict['ID_list_new'].append(str(ID_list[sample_idx]) + '_fit' + str(c + 1))
			master_dict['notes_list'].append(' ')
			master_dict['model_list'].append(' ')

		if model_name == 'Hill':

			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = hill_equation(x_fit, parameters[0], parameters[1], parameters[2], parameters[3])

			plot_func(
				figure, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'None', 'line', x_title, y_title,
				subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
				user_input_dict, master_dict)

			# Calculate statistics and apply to master dict
			scatter_statistics(x_val_list[c], y_val_list[c], parameters, model_name, alpha_level, master_dict)

			# Write parameters from fit to the master dict
			write_fit_params_to_master_dict(model_name, parameters, master_dict, func_KD_list)

		elif model_name == 'Hill_simple':
			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = hill_simple(x_fit, parameters[0], parameters[1], parameters[2])
			plot_func(figure, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'None', 'line', x_title, y_title,
					  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
					  user_input_dict, master_dict)

			# Calculate statistics and apply to master dict
			scatter_statistics(x_val_list[c], y_val_list[c], parameters, model_name, alpha_level, master_dict)

			# Write parameters from fit to the master dict
			write_fit_params_to_master_dict(model_name, parameters, master_dict, func_KD_list)

		elif model_name == '4PL':
			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = fit_4pl(x_fit, parameters[0], parameters[1], parameters[2], parameters[3])
			plot_func(figure, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'None', 'line', x_title, y_title,
					  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
					  user_input_dict, master_dict)

			# Calculate statistics and apply to master dict
			scatter_statistics(x_val_list[c], y_val_list[c], parameters, model_name, alpha_level, master_dict)

			# Write parameters from fit to the master dict
			write_fit_params_to_master_dict(model_name, parameters, master_dict, func_KD_list)

		elif model_name == 'FIDA_1to1':
			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = model_fida_1to1(x_fit, parameters[0], parameters[1], parameters[2])
			plot_func(figure, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'None', 'line', x_title, y_title,
					  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
					  user_input_dict, master_dict)

			# Calculate statistics and apply to master dict
			scatter_statistics(x_val_list[c], y_val_list[c], parameters, model_name, alpha_level, master_dict)

			# Write parameters from fit to the master dict
			write_fit_params_to_master_dict(model_name, parameters, master_dict, func_KD_list)

		elif model_name == 'FIDA_excess':

			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = model_fida_excess(x_fit, parameters[0], parameters[1], parameters[2], parameters[3])
			plot_func(figure, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'None', 'line', x_title, y_title,
					  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
					  user_input_dict, master_dict)

			# Calculate statistics and apply to master dict
			scatter_statistics(x_val_list[c], y_val_list[c], parameters, model_name, alpha_level, master_dict)

			# Write parameters from fit to the master dict
			write_fit_params_to_master_dict(model_name, parameters, master_dict, func_KD_list)

		if fitting_mode == 'Local':

			# Calculate RMSE for the full data
			rmse = rmse_calculation(x_val, y_val, parameters, model_name)
			master_dict['RMSE'].append(rmse)

		elif fitting_mode == 'Global':

			# Calculate RMSE for the "local" data i.e. meaning the data subset used in the specific fit
			rmse = rmse_calculation(x_val_list[c], y_val_list[c], parameters, model_name)
			master_dict['RMSE'].append(rmse)

		# Add the appropriate color to the table coloring list
		table_color_list_manager(color_list[color_count], plot_dict['table_color_list'])

	# If global fitting (i.e. more than one KD and R^2 has been calculated for the data) add final result to table
	if fitting_mode == 'Global':
		master_dict['ID_list_new'].append(str(ID_list[sample_idx]))
		master_dict['notes_list'].append(sample_notes[sample_idx])
		master_dict['model_list'].append(fit_models[sample_idx] + ', ' + fitting_mode)
		master_dict['fit_parameters'].append(' ')

		global_KD = statistics.mean(master_dict['func_KD_list'])
		global_R2 = statistics.stdev(master_dict['func_KD_list'])
		master_dict['KD_fit'].append(
			"{:.3f}".format(global_KD) + ' ' + u"\u00B1" + ' ' + str("{:.3f}".format(global_R2)))

		master_dict['R_square'].append(' ')
		master_dict['chi_square'].append(' ')

		master_dict['RMSE'].append(' ')

		# Add the appropriate color to the table coloring list
		table_color_list_manager(color_list[color_count], plot_dict['table_color_list'])

	# If specified by user calculate the residuals on the full data set and include in the plot.
	if param_dict['scatter_residuals'] == 'yes':
		residuals = calculate_and_plot_residuals(x_val, y_val, parameters, model_name, )

		plot_func(figure, graph_name + '_residuals', x_val, residuals, 'None', 'dots', x_title, y_title,
				  subplot_row, subplot_col, 'scatter_residuals', sample_idx, param_dict, color_list, color_count,
				  user_input_dict, master_dict)


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


def fida_text2floats(x_id, y_id, df):
	func_list = []
	func_list2 = []

	# Slice the dataframe to only include the rows that contain a concentration
	data_slice = df[df[x_id].str.contains('nM]', na=False)]
	xdata_list = data_slice[x_id].astype(str).str.replace(",", ".").tolist()
	ydata_list = pd.to_numeric(data_slice[y_id].astype(str).str.replace(",", ".")).tolist()
	for m in range(len(xdata_list)):
		conc = xdata_list[m][xdata_list[m].index('[') + 1: xdata_list[m].index(' nM]')]
		func_list.append(float(conc))

	# func_list2.append(func_list)
	# func_list2.append(ydata_list)
	return func_list, ydata_list


def scatter_statistics(x_val, y_val, parameters, model_name, alpha_level, master_dict):

	# Calculate r square and adjusted r square and put to master dict
	r_square, r_square_adj = r_square_function(x_val, y_val, parameters, model_name)

	master_dict['R_square'].append('R<sup>2</sup>=' + "{:.3f}".format(r_square) + '<br>' +
								   'R<sup>2</sup><sub>adj</sub>=' + "{:.3f}".format(r_square_adj))

	# Calculate chi square stats and put to master dict
	chi2, DoF, tail_left, tail_right, p_val = chi_square_function(x_val, y_val, parameters, model_name, alpha_level)

	master_dict['chi_square'].append('\u03A7<sup>2</sup>=' + str(chi2) + ', DoF=' + str(DoF) + \
									 '<br>' + 'Tails=' + str(tail_left) + ';' + str(tail_right) + \
									 '<br>' + 'p-value' + str(p_val) + '<br>alpha=' + str(alpha_level))


def write_fit_params_to_master_dict(model_name, parameters, master_dict, func_KD_list):

	if model_name == 'FIDA_1to1':
		func_RI = "{:.2f}".format(parameters[0])
		func_RIA = "{:.2f}".format(parameters[1])
		func_KD = "{:.2f}".format(parameters[2])
		master_dict['func_KD_list'].append(parameters[2])

		master_dict['KD_fit'].append(func_KD)
		master_dict['fit_parameters'].append('RI=' + str(func_RI) + ', RIA=' + str(func_RIA))

	elif model_name == 'FIDA_excess':
		func_RI = "{:.2f}".format(parameters[0])
		func_RIA = "{:.2f}".format(parameters[1])
		func_KD = "{:.2f}".format(parameters[2])
		master_dict['func_KD_list'].append(parameters[2])
		func_CI = "{:.3f}".format(parameters[3])

		master_dict['KD_fit'].append(func_KD)
		master_dict['fit_parameters'].append('RI=' + str(func_RI) + ', RIA=' + str(func_RIA) + ', CI=' + str(func_CI))

	elif model_name == 'Hill_simple':
		Bmin = "{:.2f}".format(parameters[0])
		Bmax = "{:.2f}".format(parameters[1])
		func_KD = "{:.2f}".format(parameters[2])
		master_dict['func_KD_list'].append(parameters[2])

		master_dict['KD_fit'].append(func_KD)
		master_dict['fit_parameters'].append('Bmin=' + str(Bmin) + ', Bmax=' + str(Bmax))

	if model_name == 'Hill':
		Bmin = "{:.2f}".format(parameters[0])
		Bmax = "{:.2f}".format(parameters[1])
		func_KD = "{:.2f}".format(parameters[2])
		master_dict['func_KD_list'].append(parameters[2])
		k_coop = "{:.3f}".format(parameters[3])

		master_dict['KD_fit'].append(func_KD)
		master_dict['fit_parameters'].append('Bmin=' + str(Bmin) + ', Bmax=' + str(Bmax) + ', k_coop=' + str(k_coop))

	elif model_name == '4PL':
		Bmin = "{:.2f}".format(parameters[0])
		Bmax = "{:.2f}".format(parameters[1])
		func_KD = "{:.2f}".format(parameters[2])
		master_dict['func_KD_list'].append(parameters[2])
		k_coop = "{:.3f}".format(parameters[3])

		master_dict['KD_fit'].append(func_KD)
		master_dict['fit_parameters'].append('Bmin=' + str(Bmin) + ', Bmax=' + str(Bmax) + ', k_coop=' + str(k_coop))


def r_square_function(x_val, y_val, parameters, model_name):
	# Converting values to arrays
	x_val = np.asarray(x_val)
	y_val = np.asarray(y_val)

	n = len(y_val)
	p = len(parameters)

	if model_name == 'FIDA_1to1':
		y_fit = model_fida_1to1(x_val, parameters[0], parameters[1], parameters[2])

	elif model_name == 'FIDA_excess':
		y_fit = model_fida_excess(x_val, parameters[0], parameters[1], parameters[2], parameters[3])

	elif model_name == 'Hill_simple':
		y_fit = hill_simple(x_val, parameters[0], parameters[1], parameters[2])

	elif model_name == 'Hill':
		y_fit = hill_equation(x_val, parameters[0], parameters[1], parameters[2], parameters[3])

	elif model_name == '4PL':
		y_fit = fit_4pl(x_val, parameters[0], parameters[1], parameters[2], parameters[3])

	else:
		raise TypeError('Fitting model not recognized!')

	# Calculate R square based on the fitted values and the observed values
	r_square = r2_score(y_val, y_fit)

	# Calculate the adjusted R square
	term1 = 1 - r_square
	term2 = n - 1
	term3 = n - p - 1
	r_square_adj = 1 - (term1 * term2) / (term3)

	return r_square, r_square_adj


def chi_square_function(x_val, y_val, parameters, model_name, alpha_level):
	# Converting values to arrays
	x_val = np.asarray(x_val)
	y_val = np.asarray(y_val)

	if model_name == 'FIDA_1to1':
		y_fit = model_fida_1to1(x_val, parameters[0], parameters[1], parameters[2])

	elif model_name == 'FIDA_excess':
		y_fit = model_fida_excess(x_val, parameters[0], parameters[1], parameters[2], parameters[3])

	elif model_name == 'Hill_simple':
		y_fit = hill_simple(x_val, parameters[0], parameters[1], parameters[2])

	elif model_name == 'Hill':
		y_fit = hill_equation(x_val, parameters[0], parameters[1], parameters[2], parameters[3])

	elif model_name == '4PL':
		y_fit = fit_4pl(x_val, parameters[0], parameters[1], parameters[2], parameters[3])

	else:
		raise TypeError('Fitting model not recognized!')

	chi2 = 0
	DoF = len(y_val) - 1

	# Calculating the chi square critical value
	for a in range(len(y_val)):
		obs = y_val[a]
		pred = y_fit[a]
		chi2 = chi2 + (obs - pred) ** 2 / pred

	# Calculating the tail values for the given alpha level
	tail_right = scipy.stats.chi2.ppf(1 - alpha_level, df=DoF)
	tail_left = scipy.stats.chi2.ppf(alpha_level, df=DoF)

	p_val = 1 - scipy.stats.chi2.sf(chi2, DoF)
	if p_val < 0.001:
		p_val_str = '<0.001'
	elif 0.001 < p_val and 0.005 > p_val:
		p_val_str = '<0.005'

	elif 0.005 < p_val and 0.01 > p_val:
		p_val_str = '<0.01'

	elif 0.01 < p_val and 0.05 > p_val:
		p_val_str = '<0.05'

	else:
		p_val_str = '=' + str("{:.2f}".format(p_val))

	# Reformating all the values for outputing
	chi2 = "{:.3f}".format(chi2)
	tail_right = "{:.2f}".format(tail_right)
	tail_left = "{:.2f}".format(tail_left)

	return chi2, DoF, tail_left, tail_right, p_val_str


def calculate_and_plot_residuals(x_val, y_val, parameters, model_name):
	df = pd.DataFrame()
	df['xs'] = x_val
	df['ys'] = y_val

	if model_name == 'FIDA_1to1':
		y_fit = model_fida_1to1(df['xs'], parameters[0], parameters[1], parameters[2])
		df['y_fit'] = y_fit

	elif model_name == 'FIDA_excess':
		y_fit = model_fida_excess(df['xs'], parameters[0], parameters[1], parameters[2], parameters[3])
		df['y_fit'] = y_fit

	elif model_name == 'Hill_simple':
		y_fit = hill_simple(x_val, parameters[0], parameters[1], parameters[2])
		df['y_fit'] = y_fit

	elif model_name == 'Hill':
		y_fit = hill_equation(x_val, parameters[0], parameters[1], parameters[2], parameters[3])
		df['y_fit'] = y_fit

	elif model_name == '4PL':
		y_fit = fit_4pl(x_val, parameters[0], parameters[1], parameters[2], parameters[3])
		df['y_fit'] = y_fit

	residuals = df['ys'] - df['y_fit']

	return residuals


def rmse_calculation(x_val, y_val, parameters, model_name):
	import math

	df = pd.DataFrame()
	df['xs'] = x_val
	df['ys'] = y_val

	if model_name == 'FIDA_1to1':
		y_fit = model_fida_1to1(df['xs'], parameters[0], parameters[1], parameters[2])

	elif model_name == 'FIDA_excess':
		y_fit = model_fida_excess(df['xs'], parameters[0], parameters[1], parameters[2], parameters[3])

	elif model_name == 'Hill_simple':
		y_fit = hill_simple(x_val, parameters[0], parameters[1], parameters[2])

	elif model_name == 'Hill':
		y_fit = hill_equation(x_val, parameters[0], parameters[1], parameters[2], parameters[3])

	elif model_name == '4PL':
		y_fit = fit_4pl(x_val, parameters[0], parameters[1], parameters[2], parameters[3])

	MSE = np.square(np.subtract(y_val, y_fit)).mean()

	RMSE = math.sqrt(MSE)

	RMSE = "{:.3f}".format(RMSE)

	return RMSE


def sigma_list_for_scatter_fit(unique_x, x_val, y_val, x_val_list, fitting_mode):
	sigma_list = []

	# If fitting mode is local then use the standard errors as sigma values.
	if fitting_mode == 'Local':
		sigma = replicate_mean_error(unique_x, x_val, y_val)[1]
		sigma = np.array(sigma)
		sigma[sigma == 0] = 1
		sigma_list.append(sigma)

	# If fitting mode is global then we just fill arrays with ones since all data points
	# will be weighted equal.
	else:
		for a in x_val_list:
			list_length = len(a)
			sigma = np.ones(list_length)
			sigma_list.append(sigma)

	return sigma_list
