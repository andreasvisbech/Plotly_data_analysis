# import relevant modules

# Import functions from other scripts
from plot_it_scripts.plotting_script import *


# from plotting_script import *

def bar_main(df, x_id, y_id, error_id, sample_idx, plot_dict, user_input_dict, param_dict, master_dict):
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

	# There might be key error if the user has only specified one x column
	try:
		xs = df[x_id].dropna()
		ys = df[y_id][pd.to_numeric(df[y_id], errors='coerce').notnull()]
		# xs = df[x_id].dropna()
		# ys = df[y_id].dropna()
	except KeyError:
		None

	# Check if error column exists
	if error_id in df.columns:

		error_values = df[error_id][pd.to_numeric(df[error_id], errors='coerce').notnull()]

		plot_func(figure, graph_name, xs, ys, error_values, plot_marker, x_title, y_title, subplot_row, subplot_col,
				  'Bar_error', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
		plot_func(plot_figure, graph_name, xs, ys, error_values, plot_marker, x_title, y_title, subplot_row,
				  subplot_col, 'Bar_error', sample_idx, param_dict, color_list, color_count, user_input_dict,
				  master_dict)

	else:
		plot_func(figure, graph_name, xs, ys, 'N/A', plot_marker, x_title, y_title, subplot_row, subplot_col, 'Bar',
				  sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
		plot_func(plot_figure, graph_name, xs, ys, 'N/A', plot_marker, x_title, y_title, subplot_row, subplot_col,
				  'Bar', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
