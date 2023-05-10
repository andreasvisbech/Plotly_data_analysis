# Import modules
from scipy.signal import savgol_filter
from sklearn.ensemble import IsolationForest

# Import relevant scripts
from plot_it_scripts.plotting_script import *

def savgol_1st_deriv(xs, ys, sample_idx, param_dict, master_dict, user_input_dict, plot_dict):

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

	# Getting the filter specs from the advanced box input
	savgol_window = param_dict['savgol_1st_deriv_window']
	savgol_pol = param_dict['savgol_1st_deriv_pol']

	# Check if the inputs have been changed from N/A. If they have been changed then apply the computation.
	# Set window to 1 and polynomium to 0 if no smoothing is needed.
	if savgol_window != 'N/A' and savgol_pol != 'N/A':

		y_val_filter = savgol_filter(ys, int(savgol_window), int(savgol_pol), deriv=1)

		plot_func(figure, graph_name + '_1st deriv', xs, y_val_filter, 'None', plot_marker, x_title, y_title,
				  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
				  user_input_dict, master_dict)
		plot_func(plot_figure, graph_name + '_1st deriv', xs, y_val_filter, 'None', plot_marker, x_title, y_title,
				  subplot_row, subplot_col, 'None', sample_idx, param_dict, color_list, color_count,
				  user_input_dict, master_dict)


def outlier_detect(xs, ys):

	from sklearn.neighbors import NearestNeighbors

	X = np.vstack((xs,ys)).T

	nbrs = NearestNeighbors(n_neighbors=2)

	nbrs.fit(X)

	distances, indexes = nbrs.kneighbors(X)

	outlier_index = np.where(distances.mean(axis=1) > 0.5)

	print(ys[outlier_index])

	None


def normalize_to_max(ys, idx, user_input_dict):

	flags = user_input_dict['python_misc'][idx].split(';')

	if 'normalize_to_max' in flags:

		# Get max value of the data
		ys_max = max(ys)

		# Update the y values by getting them as percent of the max value in the data.
		ys = (ys/ys_max)*100

	return ys


def find_nearest(arr, value):

	# The function takes an array and finds the value in the array closest to a specified value

	idx = np.abs(arr - value).argmin()

	return arr[idx]


def normalize_to_x(xs, ys, idx, user_input_dict):

	# Getting flags
	flags = user_input_dict['python_misc'][idx].split(';')

	# Make sure data is arrays
	xs = np.array(xs)
	ys = np.array(ys)

	for flag in flags:
		if 'baseline_x=' in flag:

			# The user specifies which x value should be used as the baseline point
			x_ref = flag.split('=')[1]
			x_ref = float(x_ref)

			# Finding the x value in the data closest to the user specified x value to be used as reference
			x_near = find_nearest(xs, x_ref)

			# Getting the y value that will be used for subtraction.
			y_ref = ys[xs==x_near][0]

			# Subtract the baseline value from the ys data
			ys = ys-y_ref

	return ys

def data_slice(sample_idx, xs, ys, user_input_dict):

	data_interval = user_input_dict['data_interval']

	# Make sure xs and ys are arrays
	xs = np.array(xs)
	ys = np.array(ys)

	# Slicing the data so only data within the specified data interval is included.
	if data_interval[sample_idx] != 0:

		interval_var = data_interval[sample_idx].split(';')
		interval_min = float(interval_var[0])
		interval_max = float(interval_var[1])

		# Getting the x_range as the x values that are within the specified inter
		xs_new = xs[np.logical_and(xs >= interval_min, xs <= interval_max)]
		ys_new = ys[np.logical_and(xs >= interval_min, xs <= interval_max)]

	else:
		xs_new = xs
		ys_new = ys

	return xs_new, ys_new

