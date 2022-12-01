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






