# Import modules
import pandas as pd
from scipy.stats import binned_statistic
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

# Import functions from other scripts
# from plotting_script import *
from plot_it_scripts.plotting_script import *
from plot_it_scripts.scatter_module import linear_model
#from scatter_module import linear_model


def panta_main(df, x_id, y_id, sample_idx, plot_dict, user_input_dict, param_dict, master_dict):
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
			  sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
	plot_func(plot_figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col, 'None',
			  sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

	# Setting the bin width as 1% of the data size.
	bin_width = len(xs)/100

	# Setting the Savgol filter window
	savgol_length = int(len(xs)/10)

	# Smoothing data and getting first and second derivative of the data. The last argument enables plotting of bin data,
	# interpolated data, first derivative and second derivative ('yes')
	xs_new, ys_interpol, ys_diff1, ys_diff2 = data_diff(xs, ys, bin_width, savgol_length, 2, param_dict, plot_dict,
														user_input_dict, master_dict, sample_idx)

	# Creating a linear baseline based on the second derivative.
	# The last argument enables plotting of the baseline ('yes').
	auto_baseline_coeff = auto_baseline(xs_new, ys_interpol, ys_diff2, param_dict, plot_dict, user_input_dict,
										master_dict, sample_idx)

	auto_baseline_slope = auto_baseline_coeff[0]
	auto_baseline_intercept = auto_baseline_coeff[1]

	# Determine the onset of the peak.
	# The onset is calculated as the point where the interpolated data deviates x % from the calculated baseline.
	onset = onset_detect(xs_new, ys_interpol, auto_baseline_slope, auto_baseline_intercept,
						 param_dict['peak_onset_var_deg'], auto_baseline_coeff[2])
	master_dict['peak_onset'].append(onset)

	# Determine inflection points as vertex points on the first derivative data.
	ip_points = ip_detect(xs_new, ys_diff1, param_dict['min_peak_prominence'], onset)
	master_dict['inflection_points'].append(ip_points)

	vertex_max = vertex_detect(np.asarray(xs), np.asarray(ys), param_dict['min_peak_prominence'])
	vertex_max = ["{:.1f}".format(a) for a in vertex_max]
	master_dict['vertex_max'].append(vertex_max)

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


def data_diff(x_val, y_val, bin_width, savgol_window, savgol_pol, param_dict, plot_dict, user_input_dict,
			  master_dict, sample_idx):

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

	# Create new x_val and y_val using the interpolation function.
	# This is basically to fill data between the bin points i.e. making data big again.
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


def auto_baseline(x_val, y_val, y_val_diff2, param_dict, plot_dict, user_input_dict,
			  master_dict, sample_idx):

	# Creating local variables for plotting
	figure = plot_dict['figure']
	graph_name = plot_dict['graph_names'][sample_idx]
	x_title = user_input_dict['x_titles'][sample_idx]
	y_title = user_input_dict['y_titles'][sample_idx]
	subplot_row = user_input_dict['subplot_row'][sample_idx]
	subplot_col = user_input_dict['subplot_col'][sample_idx]
	color_list = plot_dict['color_list']
	color_count = plot_dict['color_count']
	plot_intermediate = param_dict['panta_intermediate_plot']

	# Converting x values and 2nd derivative data to lists.
	x_val_list = list(x_val)
	y_val_list = list(y_val_diff2)

	# Determine baseline start as the point on the 2nd derivative data where the following
	# 10 data points are all within 3% of the max value of that derivative i.e. very close to zero
	percent_cutoff = 0.03
	for a in range(len(y_val_list)):
		data_window = y_val_list[a:a + 10]
		counter = 0
		for b in range(len(data_window)):
			if abs(data_window[b]) < percent_cutoff * max(y_val_diff2):
				counter = counter + 1
			else:
				counter = counter + 0

		if counter == 10:
			baseline_start_idx = a
			break

		else:
			baseline_start_idx = 0

	# The end of the autobaseline
	for c in range(baseline_start_idx, len(y_val_list)):
		data_window = y_val_list[c:c + 10]
		counter = 0
		for d in range(len(data_window)):
			if abs(data_window[d]) > 0.1 * max(y_val_diff2) and c > baseline_start_idx:  # Make sure the index of baseline end is bigger than the start index
				counter = counter + 1
			else:
				counter = counter + 0
		if counter == 10:
			baseline_end_idx = c
			break

		else:
			baseline_end_idx = len(y_val_list) - 1

	baseline_x = x_val_list[baseline_start_idx:baseline_end_idx]
	baseline_y = list(y_val)[baseline_start_idx:baseline_end_idx]

	parameters, covariance = curve_fit(linear_model, baseline_x, baseline_y, method='trf', maxfev=10000)

	# y_pred = linear_model(np.asarray(baseline_x), parameters[0], parameters[1])
	# r_square = r2_score(baseline_y, y_pred)
	# print('R square: ' + str(r_square))

	y_pred = linear_model(x_val, parameters[0], parameters[1])

	if plot_intermediate == 'yes':
		plot_func(figure, graph_name + '_auto_baseline', x_val, y_pred, 'None', 'line', x_title, y_title, subplot_row,
				  subplot_col, 'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

	# Return coefficients for the linear fit
	return [parameters[0], parameters[1], baseline_x[-1]]


def onset_detect(x_val, y_val, slope, intercept, detect_limit, baseline_cutoff):
	x_val_list = list(x_val)
	y_val_list = list(y_val)

	for a in range(len(y_val_list)):
		data_obs = y_val_list[a]
		data_calc = linear_model(x_val_list[a], slope, intercept)

		# The onset is defined as as the point where the observed data deviates more and 0.5% from the linear baseline.
		# It is further required that the x value is higher than the baseline cutoff to avoid predicting onsets on noise in the beginning of the dataset.
		if data_obs > data_calc * (1 + (detect_limit / 100)) and x_val_list[a] > baseline_cutoff:
			onset_val = x_val_list[a]
			onset_val = "{:.1f}".format(onset_val)
			break

		else:
			onset_val = 0

	return onset_val


def ip_detect(x_val, y_val, min_peak_prominence, onset_cutoff):

	IP_list = []

	peaks_temp = vertex_detect(x_val, y_val, min_peak_prominence)

	for a in peaks_temp:
		if float(a) > float(onset_cutoff):
			my_IP = "{:.1f}".format(a)
			IP_list.append(my_IP)

	return IP_list


def vertex_detect(x_val, y_val, min_peak_prominence):

	# Getting indexes of peaks
	peaks, _ = find_peaks(y_val, prominence=min_peak_prominence)

	# Getting the corresponding x data
	peaks = x_val[peaks]

	return peaks


def vertex_detect2(x_val, y_val, window_size):

    vertex_min = []
    vertex_max = []

    x_values = x_val.tolist()
    y_values = y_val.tolist()

    if len(x_values) != len(y_values):
            print('ERROR! Different length of arrays')
    else:
        for a in range(window_size, len(x_values)-window_size):
            window_values = y_values[a-window_size:a+window_size]
            if y_values[a] == max(window_values):
                vertex_max.append("{:.1f}".format(x_values[a]))
            elif y_values[a] == min(window_values):
                vertex_min.append("{:.1f}".format(x_values[a]))

    return [vertex_min, vertex_max]