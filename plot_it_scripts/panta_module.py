# Import modules
from scipy.stats import binned_statistic
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

# Import functions from other scripts
# from plotting_script import *
from plot_it_scripts.plotting_script import *


def panta_main(df, x_id, y_id, sample_idx, plot_dict, user_input_dict, param_dict):
	# Creating local variables for plotting
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

	# Getting data and slicing if needed
	xs = df[x_id][pd.to_numeric(df[x_id], errors='coerce').notnull()]
	ys = df[y_id][pd.to_numeric(df[y_id], errors='coerce').notnull()]

	# Slicing the data so only data within the specified data interval is included.
	xs, ys = panta_data_slice(x_id, y_id, sample_idx, xs, ys, user_input_dict)

	# Plotting the data values into interactive plot and the static "plotting plot"
	plot_func(figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col, 'None',
			  sample_idx, param_dict, color_list, color_count, user_input_dict)
	plot_func(plot_figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col, 'None',
			  sample_idx, param_dict, color_list, color_count, user_input_dict)

	# Smoothing data and getting first and second derivative of the data. The last argument enables plotting of bin data,
	# interpolated data, first derivative and second derivative ('yes')
	xs_new, ys_interpol, ys_diff1, ys_diff2 = data_diff(xs, ys, 10, 81, 2, param_dict, plot_dict, user_input_dict,
														sample_idx)


def panta_data_slice(x_id, y_id, sample_idx, xs, ys, user_input_dict):
	data_interval = user_input_dict['data_interval']

	# Slicing the data so only data within the specified data interval is included.
	if data_interval[sample_idx] != 0:
		interval_var = data_interval[sample_idx].split(';')

		xsys_interval = pd.concat([xs, ys], axis=1)
		xsys_interval_slice = xsys_interval[
			(xsys_interval[x_id] >= float(interval_var[0])) & (xsys_interval[x_id] <= float(interval_var[1]))]

		xs = xsys_interval_slice[x_id]
		ys = xsys_interval_slice[y_id]

	return xs, ys


def data_diff(x_val, y_val, bin_width, savgol_window, savgol_pol, param_dict, plot_dict, user_input_dict, sample_idx):
	# Creating local variables for plotting
	figure = plot_dict['figure']
	graph_name = plot_dict['graph_names'][sample_idx]
	plot_marker = user_input_dict['plot_markers'][sample_idx]
	x_title = user_input_dict['x_titles'][sample_idx]
	y_title = user_input_dict['y_titles'][sample_idx]
	subplot_row = user_input_dict['subplot_row'][sample_idx]
	subplot_col = user_input_dict['subplot_col'][sample_idx]
	color_list = plot_dict['color_list']
	color_count = plot_dict['color_count']
	plot_intermediate = param_dict['panta_intermediate_plot']

	bin_size = len(y_val) / bin_width
	bin_y = binned_statistic(list(x_val), list(y_val), statistic='mean', bins=bin_size)[0]  # Create bin for y values
	bin_x = binned_statistic(list(x_val), list(x_val), statistic='mean', bins=bin_size)[0]  # Create bin for x values

	# Getting min and max values of the x values and y values
	xs_min = list(x_val)[0]
	xs_max = list(x_val)[-1]
	ys_min = list(y_val)[0]
	ys_max = list(y_val)[-1]

	# Update the bins to include first and last value of the original x_val and y_val.
	# This is to make sure the interpolation range can accomodate all values of the original data.
	bin_x = np.insert(bin_x, 0, xs_min)
	bin_x = np.insert(bin_x, len(bin_x), xs_max)
	bin_y = np.insert(bin_y, 0, ys_min)
	bin_y = np.insert(bin_y, len(bin_y), ys_max)

	# Create interpolation function based on the binned data
	f = interp1d(bin_x, bin_y)

	# Create new x_val and y_val using the interpolation function. This is basically to fill data between the bin points i.e. making data big again.
	# I also create new x values that are evenly spaced in the original x value range
	x_val_new = np.linspace(list(x_val)[0], list(x_val)[-1], len(x_val))
	y_val_new = f(list(x_val_new))

	# Smooth the interpolated data using savgol filter and take first derivative.
	y_val_filter = savgol_filter(y_val_new, savgol_window, savgol_pol, deriv=1)

	# Smooth the first derivative data and take the second derivative.
	y_val_filter2 = savgol_filter(y_val_filter, savgol_window, savgol_pol, deriv=1)

	if plot_intermediate == 'yes':
		plot_func(figure, graph_name + '_bin', bin_x, bin_y, 'None', 'line', x_title, y_title, subplot_row, subplot_col,
				  'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
		plot_func(figure, graph_name + '_interpolated', x_val_new, y_val_new, 'None', 'line', x_title, y_title,
				  subplot_row,
				  subplot_col, 'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
		plot_func(figure, graph_name + '_1st deriv', x_val_new, y_val_filter, 'None', 'line', x_title, y_title,
				  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count, user_input_dict,
				  master_dict)
		plot_func(figure, graph_name + '_2nd deriv', x_val_new, y_val_filter2, 'None', 'line', x_title, y_title,
				  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count, user_input_dict,
				  master_dict)

	return x_val_new, y_val_new, y_val_filter, y_val_filter2
