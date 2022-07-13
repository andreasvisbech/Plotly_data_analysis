# Import relevant modules
import statistics
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from scipy.stats import chisquare

# Import functions from other scripts
from plot_it_scripts.plotting_script import *


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
	unique_x = list(dict.fromkeys(
		xs))  # Create a list with unique concentration values. These should include ALL relevant x values.

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
	ID_list = user_input_dict['ID_list']
	sample_notes = user_input_dict['sample_notes']

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

	# If local fitting approach the data will be averaged and plotted with error bars. Else the full data set will be plotted.
	if fit_modes[sample_idx] == 'Local':
		y_mean, std_dev = replicate_mean_error(unique_x, xs, ys)
		plot_func(figure, graph_name, unique_x, y_mean, std_dev, plot_marker, x_title, y_title, subplot_row,
				  subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
				  user_input_dict, master_dict)
		plot_func(plot_figure, graph_name, unique_x, y_mean, std_dev, plot_marker, x_title, y_title, subplot_row,
				  subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
				  user_input_dict, master_dict)

	else:
		plot_func(figure, graph_name, xs, ys, np.zeros(len(xs)), plot_marker, x_title, y_title, subplot_row,
				  subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
				  user_input_dict, master_dict)
		plot_func(plot_figure, graph_name, xs, ys, np.zeros(len(xs)), plot_marker, x_title, y_title, subplot_row,
				  subplot_col, 'Scatter_error', sample_idx, param_dict, color_list, color_count,
				  user_input_dict, master_dict)

	# Extract boundaries for the fit from the excel sheet and turn to list of floats. If no interval is supplied code will automaticallt assign 0;0.
	fitting_interval = fit_intervals[sample_idx].split(';')
	for j in range(0, len(fitting_interval)):
		fitting_interval[j] = float(
			fitting_interval[j].replace(",", "."))  # Replace commas with proper dots to get proper numbers.

	# Running fitting function and extracting parameters.
	if fit_models[sample_idx] == 'Hill':
		scatter_fit(hill_equation, 'Hill', xs, ys, fit_modes[sample_idx], fitting_interval[0], fitting_interval[1],
					sample_idx, param_dict, color_list, master_dict, user_input_dict, plot_dict)

	elif fit_models[sample_idx].lower() == 'hill_simple':
		scatter_fit(hill_simple, 'Hill_simple', xs, ys, fit_modes[sample_idx], fitting_interval[0], fitting_interval[1],
					sample_idx, param_dict, color_list, master_dict, user_input_dict, plot_dict)

	elif fit_models[sample_idx].lower() == '4pl':
		scatter_fit(fit_4pl, '4PL', xs, ys, fit_modes[sample_idx], fitting_interval[0], fitting_interval[1], sample_idx,
					param_dict, color_list, master_dict, user_input_dict, plot_dict)

	elif fit_models[sample_idx].lower() in ['fida_1to1', '1to1', 'one2one']:
		scatter_fit(model_fida_1to1, 'FIDA_1to1', xs, ys, fit_modes[sample_idx], fitting_interval[0],
					fitting_interval[1], sample_idx, param_dict, color_list, master_dict, user_input_dict, plot_dict)

	elif fit_models[sample_idx].lower() in ['fida_excess', 'excess']:
		scatter_fit(model_fida_excess, 'FIDA_excess', xs, ys, fit_modes[sample_idx], fitting_interval[0],
					fitting_interval[1], sample_idx, param_dict, color_list, master_dict, user_input_dict, plot_dict)


def scatter_plot_without_fit(master_dict, sample_idx, user_input_dict, plot_dict, xs, ys, param_dict):
	# Appending dummy values to the master dict to avoid mispairing of cells in the output table
	master_dict['ID_list_new'].append(user_input_dict['ID_list'][sample_idx])
	master_dict['notes_list'].append(user_input_dict['sample_notes'][sample_idx])
	master_dict['model_list'].append(' ')
	master_dict['fit_parameters'].append(' ')
	master_dict['R_square'].append(' ')
	master_dict['KD_fit'].append(' ')

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

	plot_func(figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col, 'None',
			  sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
	plot_func(plot_figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col, 'None',
			  sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)


def replicate_mean_error(unique_x, redundant_x, y_val):
	"""
	The function takes a list of unique concentrations, the redundant concentration list and a list of y values corresponding to the redundant concentration list.
	The function then calculates appropriate mean values and the standard deviation on the values used for the mean calculation.

	:param unique_conc:
	:param redundant_conc:
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


def scatter_fit(model_function, model_name, x_val, y_val, fitting_mode, fitting_min, fitting_max, sample_idx,
				param_dict, color_list, master_dict, user_input_dict, plot_dict):
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

	# Loading local variables for analysis
	ID_list = user_input_dict['ID_list']
	sample_notes = user_input_dict['sample_notes']
	fit_models = user_input_dict['fit_models']

	out_dict = {}

	x_val_list = []
	y_val_list = []

	unique_x = list(dict.fromkeys(x_val))

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

	if fitting_mode == 'Local':

		# In local fitting mode the list of unique x values is the only x list needed.
		x_val_list.append(unique_x)

		# Calculate mean and std error on replicate values in the data
		y_mean = replicate_mean_error(unique_x, x_val, y_val)[0]
		y_val_list.append(y_mean)
		y_mean_err = replicate_mean_error(unique_x, x_val, y_val)[1]
		out_dict['y_mean'] = y_mean

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
			conc_indices = [k for k, x in enumerate(list(x_val)) if
							x == j]  # Getting indices of the specific concentration
			for k in range(len(conc_indices)):
				x_val_list[k].append(list(x_val)[conc_indices[k]])
				y_val_list[k].append(list(y_val)[conc_indices[k]])

	func_KD_list = []
	func_R2_list = []

	for c in range(len(x_val_list)):
		parameters, covariance = curve_fit(model_function, list(x_val_list[c]), list(y_val_list[c]), bounds=bound_param,
										   method='trf', maxfev=10000)

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
			plot_func(figure, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'None', 'line', x_title, y_title,
					  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
					  user_input_dict, master_dict)

			# Calcularing R squared value
			y_fit_small = hill_equation(x_val_list[c], parameters[0], parameters[1], parameters[2], parameters[3])
			method_name(c, func_KD_list, func_R2_list, parameters, y_fit_small, y_val_list, master_dict)

		elif model_name == 'Hill_simple':
			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = hill_simple(x_fit, parameters[0], parameters[1], parameters[2])
			plot_func(figure, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'None', 'line', x_title, y_title,
					  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
					  user_input_dict, master_dict)

			# Calcularing R squared value
			y_fit_small = hill_simple(np.asarray(x_val_list[c]), parameters[0], parameters[1], parameters[2])
			r_square = r2_score(y_val_list[c], y_fit_small)
			master_dict['R_square'].append("{:.3f}".format(r_square))
			func_R2_list.append(r_square)

			func_Bmin = "{:.3f}".format(parameters[0])
			func_Bmax = "{:.3f}".format(parameters[1])
			func_KD = "{:.3f}".format(parameters[2])
			func_KD_list.append(parameters[2])

			master_dict['KD_fit'].append(func_KD)
			master_dict['fit_parameters'].append('Bmin=' + str(func_Bmin) + ', Bmax=' + str(func_Bmax))

		elif model_name == '4PL':
			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = fit_4pl(x_fit, parameters[0], parameters[1], parameters[2], parameters[3])
			plot_func(figure, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'None', 'line', x_title, y_title,
					  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
					  user_input_dict, master_dict)

			# Calcularing R squared value
			y_fit_small = fit_4pl(x_val_list[c], parameters[0], parameters[1], parameters[2], parameters[3])
			method_name(c, func_KD_list, func_R2_list, parameters, y_fit_small, y_val_list, master_dict)

		elif model_name == 'FIDA_1to1':
			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = model_fida_1to1(x_fit, parameters[0], parameters[1], parameters[2])
			plot_func(figure, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'None', 'line', x_title, y_title,
					  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
					  user_input_dict, master_dict)

			# Calcularing R squared value
			y_fit_small = model_fida_1to1(np.asarray(x_val_list[c]), parameters[0], parameters[1], parameters[2])
			r_square = r2_score(y_val_list[c], y_fit_small)
			master_dict['R_square'].append("{:.3f}".format(r_square))
			func_R2_list.append(r_square)

			func_RI = "{:.3f}".format(parameters[0])
			func_RIA = "{:.3f}".format(parameters[1])
			func_KD = "{:.3f}".format(parameters[2])
			func_KD_list.append(parameters[2])

			master_dict['KD_fit'].append(func_KD)
			master_dict['fit_parameters'].append('RI=' + str(func_RI) + ', RIA=' + str(func_RIA))

		elif model_name == 'FIDA_excess':

			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = model_fida_excess(x_fit, parameters[0], parameters[1], parameters[2], parameters[3])
			plot_func(figure, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'None', 'line', x_title, y_title,
					  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
					  user_input_dict, master_dict)

			# Calcularing R squared value
			y_fit_small = model_fida_excess(np.asarray(x_val_list[c]), parameters[0], parameters[1], parameters[2],
											parameters[3])
			r_square = r2_score(y_val_list[c], y_fit_small)
			master_dict['R_square'].append("{:.3f}".format(r_square))
			func_R2_list.append(r_square)

			func_RI = "{:.3f}".format(parameters[0])
			func_RIA = "{:.3f}".format(parameters[1])
			func_KD = "{:.3f}".format(parameters[2])
			func_KD_list.append(parameters[2])
			func_CI = "{:.3f}".format(parameters[3])

			master_dict['KD_fit'].append(func_KD)
			master_dict['fit_parameters'].append(
				'RI=' + str(func_RI) + ', RIA=' + str(func_RIA) + ', CI=' + str(func_CI))

	# If global fitting (i.e. more than one KD and R^2 has been calculated for the data) add final result to table
	if fitting_mode == 'Global':
		master_dict['ID_list_new'].append(str(ID_list[sample_idx]))
		master_dict['notes_list'].append(sample_notes[sample_idx])
		master_dict['model_list'].append(fit_models[sample_idx] + ', ' + fitting_mode)
		master_dict['fit_parameters'].append(' ')

		global_KD = statistics.mean(func_KD_list)
		global_R2 = statistics.stdev(func_KD_list)
		master_dict['KD_fit'].append(
			"{:.3f}".format(global_KD) + ' ' + u"\u00B1" + ' ' + str("{:.3f}".format(global_R2)))

		master_dict['R_square'].append(' ')


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


def method_name(c, func_KD_list, func_R2_list, parameters, y_fit_small, y_val_list, master_dict):
	r_square = r2_score(y_val_list[c], y_fit_small)
	master_dict['R_square'].append("{:.3f}".format(r_square))
	func_R2_list.append(r_square)
	func_Bmin = "{:.3f}".format(parameters[0])
	func_Bmax = "{:.3f}".format(parameters[1])
	func_KD = "{:.3f}".format(parameters[2])
	func_KD_list.append(parameters[2])
	func_k_coop = "{:.3f}".format(parameters[3])
	master_dict['KD_fit'].append(func_KD)
	master_dict['fit_parameters'].append(
		'Bmin=' + str(func_Bmin) + ', Bmax=' + str(func_Bmax) + ', k_coop=' + str(func_k_coop))


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


None
