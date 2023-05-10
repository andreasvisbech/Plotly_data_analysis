# import relevant modules


# Import functions from other scripts
from plot_it_scripts.plotting_script import *


# from plotting_script import *

def bioanalyzer_main(df, xs, ys, sample_idx, x_id, y_id, param_dict, master_dict, user_input_dict, plot_dict):

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

	# Plotting the data values
	plot_func(figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col,
			  'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
	plot_func(plot_figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col,
			  'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

def bioanalyzer_clean_input(pd_column):

	data = pd_column.values.tolist()

	# Create empty lists for populating with data values
	xs = []
	ys = []

	for i in data:

		if str(i).count(',') >= 2:

			# Find indexes of commas
			idx_list = []
			for idx, j in enumerate(i):
				if j == ',':
					idx_list.append(idx)

			# If two commas are present we want to split according to the first.
			if str(i).count(',') == 2:
				split_idx = int(idx_list[0])

			# If three commas are present we want to split according to the second
			elif str(i).count(',') == 3:
				split_idx = int(idx_list[1])

			# Getting the x values as left part of the string
			x = i[0:split_idx].replace(',' , '.')
			xs.append(float(x))

			# Getting the y values as the right part of the string
			y = i[split_idx+1:len(i)].replace(',' , '.')
			ys.append(float(y))

	# Convert the data lists to arrays
	xs = np.array(xs)
	ys = np.array(ys)

	bioanalyzer_data_QC(xs)

	return xs, ys

def bioanalyzer_data_QC(xs):

	# Purpose of function is to check if some x values have disappeared when importing from BioAnalyzer
	# software to the excel sheet

	# create a list for storing wronly imported x values
	wrong_x_list = []

	steps = np.diff(xs)

	for idx, i in enumerate(steps):

		if round(i, 3) != 0.05:

			wrong_x = (xs[idx] + xs[idx+1])/2

			wrong_x_list.append(wrong_x)

	if len(wrong_x_list) > 0:
		print('Please check the following x values: ' + str(wrong_x_list))








