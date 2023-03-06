# Import relevant modules
import plotly.graph_objects as go
import pandas as pd
import matplotlib.colors
import numpy as np

global used_graph_names
used_graph_names = []

global used_color_IDs
used_colors = []

global secondary_y_counter



# color_count_global = 0

def plot_func(figure, graph_name, x_val, y_val, std_dev, marker, x_title, y_title, subplot_row, subplot_col, comment, i,
			  param_dict, color_list, color_count, user_input_dict, master_dict):
	"""
    The plotting function adds trace to the specified figure. 'graph_name' is the name of the added trace.
    'x_val' and 'y_val' are the plotted data. 'marker' specifies line or dots. 'x_title' and 'y_title' is the
    axis names. 'subplot_row' and 'subplpot_col' specify which subplot to add to. 'comment' is used to specify
    special plotting cases e.g. baselines.

    :param figure:
    :param graph_name:
    :param x_val:
    :param y_val:
    :param marker:
    :param x_title:
    :param y_title:
    :param subplot_row:
    :param subplot_col:
    :param comment:
    :param i:
    :param param_dict:
    :param color_list: A list of colors in hex format
    :return:
    """

	# Set color for the trace
	user_input_dict['plot_color'] = color_list[color_count]

	# Define template for hover label
	my_hover_template = graph_name + '<extra></extra>' + '<br>x: %{x}' + '<br>y: %{y}<br>' + 'Note: ' + \
						user_input_dict['sample_notes'][i]

	# Determining the "number id" for the individual subplot. This can be used for targeting customization
	# such as axis labels to this subplot.
	# Note this code used to be defined after adding traces
	subplot_id = ((subplot_row - 1) * user_input_dict['subplot_col_count']) + subplot_col

	#ax_id_x, ax_id_y = get_ax_id(subplot_id, user_input_dict)
	#if subplot_id == 1:
	#	ax_id = ''
	#else:
	#	ax_id = str(subplot_id)

	# Customize the plot if specified by user
	# Note this code used to be defined after adding the trace
	#if user_input_dict['python_misc'][i] != 'None':
	plot_customize(figure, i, user_input_dict, param_dict, subplot_id)

	# Setting titles on the subplots
	# Note this code used to be defined after adding traces.
	#figure['layout']['xaxis' + ax_id_x]['title'] = x_title
	#figure['layout']['yaxis' + ax_id_y]['title'] = y_title

	# Keep track of which graph names have already been used.
	# If the graph name has already been used we don't include it in the legend again
	used_graph_names.append(graph_name)
	if used_graph_names.count(graph_name) > 1 and figure != 'plot_figure':
		legend_show = False
	else:
		legend_show = True

	if comment == 'None':
		if marker.lower() in ['lines', 'line']:
			figure.add_trace(
				go.Scatter(
					name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show,
					hovertemplate=my_hover_template,
					line=dict(width=param_dict['graph width'], color=user_input_dict['plot_color'])),
				row=subplot_row, col=subplot_col, secondary_y=user_input_dict['secondary_y']
			)

		elif marker.lower() in ['dot', 'dots', 'marker', 'markers']:
			figure.add_trace(
				go.Scatter(
					name=graph_name, x=x_val, y=y_val, mode='markers', legendgroup=graph_name, showlegend=legend_show,
					hovertemplate=my_hover_template,
					marker=dict(size=param_dict['marker_size'], color=user_input_dict['plot_color'])),
				row=subplot_row, col=subplot_col
			)

	elif comment == 'Scatter_error':

		rgb_color = matplotlib.colors.to_rgb(user_input_dict['plot_color'])
		rgb_color = 'rgba(' + str(rgb_color[0]) + ',' + str(rgb_color[1]) + ',' + str(rgb_color[2]) + ', 0.3)'

		y_val_upper = np.add(y_val, std_dev)
		y_val_lower = np.subtract(y_val, std_dev)

		if marker.lower() in ['lines', 'line']:

			if param_dict['error_marker'] == 'Error bands':

				figure.add_trace(
					go.Scatter(
						name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show,
						line=dict(width=param_dict['graph width'], color=user_input_dict['plot_color'])),
					row=subplot_row, col=subplot_col
				)

				figure.add_trace(
					go.Scatter(
						name=graph_name + '_upper', x=x_val, y=y_val_upper, mode='lines', legendgroup=graph_name,
						showlegend=False, line=dict(width=0, color=user_input_dict['plot_color'])),
					row=subplot_row, col=subplot_col
				)

				figure.add_trace(
					go.Scatter(
						name=graph_name + '_lower', x=x_val, y=y_val_lower, mode='lines', legendgroup=graph_name,
						showlegend=False, fillcolor=rgb_color, fill='tonexty',
						line=dict(width=0, color=user_input_dict['plot_color'])),
					row=subplot_row, col=subplot_col
				)


			else:
				figure.add_trace(
					go.Scatter(
						name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show,
						hovertemplate=my_hover_template,
						error_y=dict(type='data', array=std_dev, visible=True),
						line=dict(width=param_dict['graph width'], color=user_input_dict['plot_color'])),
					row=subplot_row, col=subplot_col
				)

		elif marker.lower() in ['dot', 'dots', 'marker', 'markers']:

			if param_dict['error_marker'] == 'Error bands':

				figure.add_trace(
					go.Scatter(
						name=graph_name, x=x_val, y=y_val, mode='markers', legendgroup=graph_name,
						showlegend=legend_show,
						hovertemplate=my_hover_template,
						marker=dict(size=param_dict['marker_size'], color=user_input_dict['plot_color'])),
					row=subplot_row, col=subplot_col
				)

				figure.add_trace(
					go.Scatter(
						name=graph_name + '_upper', x=x_val, y=y_val_upper, mode='lines', legendgroup=graph_name,
						showlegend=False,
						hovertemplate=my_hover_template,
						line=dict(width=0, color=user_input_dict['plot_color'])),
					row=subplot_row, col=subplot_col
				)

				figure.add_trace(
					go.Scatter(
						name=graph_name + '_lower', x=x_val, y=y_val_lower, mode='lines', legendgroup=graph_name,
						showlegend=False, fillcolor=rgb_color, fill='tonexty',
						hovertemplate=my_hover_template,
						line=dict(width=0, color=user_input_dict['plot_color'])),
					row=subplot_row, col=subplot_col
				)

			else:

				figure.add_trace(
					go.Scatter(
						name=graph_name, x=x_val, y=y_val, mode='markers', legendgroup=graph_name,
						meta='scatter',
						showlegend=legend_show,
						hovertemplate=my_hover_template,
						error_y=dict(type='data', array=std_dev, visible=True),
						marker=dict(size=param_dict['marker_size'], color=user_input_dict['plot_color'])),
					row=subplot_row, col=subplot_col
				)

	elif comment == 'scatter_residuals':

		figure.add_trace(
			go.Scatter(
				name=graph_name, x=x_val, y=y_val, mode='markers', marker_symbol='x',
				legendgroup=graph_name, showlegend=legend_show,
				hovertemplate=my_hover_template,
				marker=dict(size=param_dict['marker_size'], color=user_input_dict['plot_color'])),
			row=subplot_row, col=subplot_col
		)

	elif comment == 'Bar':
		figure.add_trace(
			go.Bar(name=graph_name, x=x_val, y=y_val, legendgroup=graph_name, showlegend=legend_show,
				   hovertemplate=my_hover_template,
				   marker_line=dict(width=param_dict['bar_edge_width'], color='black'),
				   marker=dict(color=user_input_dict['plot_color'])),
			row=subplot_row, col=subplot_col
		)

	elif comment == 'Boxplot':
		figure.add_trace(
			go.Box(name=graph_name, y=y_val)
		)

	elif comment == 'Boxplot_all':
		figure.add_trace(
			go.Box(name=graph_name, y=y_val, boxpoints='all')
		)

	elif comment == 'Bar_error':
		figure.add_trace(
			go.Bar(
				name=graph_name, x=x_val, y=y_val, legendgroup=graph_name, showlegend=legend_show,
				hovertemplate=my_hover_template,
				error_y=dict(type='data', array=std_dev, visible=True),
				marker_line=dict(width=param_dict['bar_edge_width'], color='black'),
				marker=dict(color=user_input_dict['plot_color'])),
			row=subplot_row, col=subplot_col
		)

	elif comment == 'AKTA_baseline':
		figure.add_trace(
			go.Scatter(name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show,
					   hovertemplate=my_hover_template,
					   line=dict(width=param_dict['graph width'], color='rgb(255,0,0)')),
			row=subplot_row, col=subplot_col, secondary_y=user_input_dict['secondary_y']
		)

	elif comment == 'AKTA_fraction':
		figure.add_trace(
			go.Scatter(
				name=graph_name, x=x_val, y=y_val, mode='lines', fill='tonexty', legendgroup=graph_name,
				showlegend=legend_show,
				hovertemplate=my_hover_template,
				line=dict(width=param_dict['graph width'], color=user_input_dict['plot_color'])),
			row=subplot_row, col=subplot_col
		)

		# 'tozeroy'

	elif comment == 'Octet':
		# Re-define template for hover label
		my_hover_template = graph_name + '<extra></extra>' + '<br>x: %{x}' + '<br>y: %{y}<br>' + 'Note: ' + \
							user_input_dict['sensor']

		figure.add_trace(
			go.Scatter(
				name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show,
				hovertemplate=my_hover_template,
				line=dict(width=param_dict['graph width'], color=user_input_dict['plot_color'])),
			row=subplot_row, col=subplot_col
		)

	elif comment == 'Octet_fit':
		# Re-define template for hover label
		my_hover_template = graph_name + '<extra></extra>' + '<br>x: %{x}' + '<br>y: %{y}<br>' + 'Note: ' + \
							user_input_dict['sensor'] + '_fit'

		figure.add_trace(
			go.Scatter(
				name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show,
				hovertemplate=my_hover_template,
				line=dict(width=param_dict['graph width'], color='rgb(255,0,0)')),
			row=subplot_row, col=subplot_col
		)



	# Determining the "number id" for the individual subplot. This can be used for targeting customization
	# such as axis labels to this subplot.
	#subplot_id = ((subplot_row - 1) * user_input_dict['subplot_col_count']) + subplot_col
	#if subplot_id == 1:
	#	ax_id = ''
	#else:
	#	ax_id = str(subplot_id)

	#figure['layout']['xaxis' + ax_id]['title'] = x_title
	#figure['layout']['yaxis' + ax_id]['title'] = y_title

	#print(figure)

	figure.update_layout(
		template=param_dict['plot_template'],
		xaxis=dict(title_font=dict(size=param_dict['xaxis_title_font_size']),
				   tickfont_size=param_dict['xaxis_ticks_font_size']),
		yaxis=dict(title_font=dict(size=param_dict['yaxis_title_font_size']),
				   tickfont_size=param_dict['yaxis_ticks_font_size'])
	)

	figure.update_layout(
		hoverlabel=dict(
			font_size=12,
			font_family="Arial"
		)
	)

	#if user_input_dict['python_misc'][i] != 'None':
	#	plot_customize(figure, user_input_dict['python_misc'][i], ax_id)



def plot_customize(figure, sample_idx, user_input_dict, param_dict, subplot_id):

	x_title = user_input_dict['x_titles'][sample_idx]
	y_title = user_input_dict['y_titles'][sample_idx]
	x_title_font_size = param_dict['xaxis_title_font_size']
	y_title_font_size = param_dict['yaxis_title_font_size']

	flags = user_input_dict['python_misc'][sample_idx]

	flag_list = flags.split(';')

	# Pre-define some of the variables in case the length of the flag list is zero
	user_input_dict['secondary_y'] = False

	ax_id_x, ax_id_y = get_ax_id(subplot_id, user_input_dict)

	# We first go through the flag list the first time to see if secondary y is defined.
	# This is needed because the ax ids must be updated accordingly.
	for a in range(len(flag_list)):
		if flag_list[a].lower() == 'secondary_y':
			user_input_dict['secondary_y'] = True
			user_input_dict['second_y_traces'].append(subplot_id)

			if ax_id_y != '':
				ax_id_y = str(int(ax_id_y)+1)
			else:
				ax_id_y = '2'

	# We go through the flags again looking for other inputs.
	for b in range(len(flag_list)):
		if flag_list[b].lower() == 'logx':
			figure['layout']['xaxis' + ax_id_x]['type'] = 'log'
			figure['layout']['xaxis' + ax_id_x]['dtick'] = 1

		if flag_list[b].lower() == 'logy':
			figure['layout']['yaxis' + ax_id_y]['type'] = 'log'
			figure['layout']['yaxis' + ax_id_y]['dtick'] = 1

		if 'color' in flag_list[b].lower():
			# If a color is defined by the user this will be used over the default color.
			new_color = flag_list[b].lower().split('=')[1]
			user_input_dict['plot_color'] = new_color


	# Setting titles on the subplots
	# Note this code used to be defined after adding traces.
	figure['layout']['xaxis' + ax_id_x]['title'] = {'text': x_title, 'font': {'size':x_title_font_size}}
	figure['layout']['yaxis' + ax_id_y]['title'] = {'text': y_title, 'font': {'size':y_title_font_size}}


def plotly_buttons(plot_dict):
	# Getting figure from plotting dictionary
	figure = plot_dict['figure']

	figure.update_layout(
		updatemenus=[
			dict(
				type="buttons",
				direction="left",
				buttons=list([dict(args=["type","scatter"], label="Graphs", method="restyle"),
							  dict(args=["type", "table"], label="Stats", method="restyle")]),
				pad={"r": 10, "t": 10}, showactive=True, x=0.9, xanchor="left", y=1.1, yanchor="top"), ])


def table_plot(plot_dict, col_names_list, col_values_list, user_input_dict, param_dict):
	"""
    The function is intended for making the table going into the html output.
    Variable "figure" specifies the figure that the table are added to.
    Variable "col_names_list" specifies a list the the column headers that should go into the table.
    Variable "col_values_list" is a list of lists. Each list contains data for an individual column

    :param col_names_list:
    :param col_values_list:
    param_dict['plot_fig_xmin']
    :return:
    """

	# Read local variables into function
	figure = plot_dict['figure']
	ID_list = user_input_dict['ID_list']
	python_misc = user_input_dict['python_misc']

	# Getting the list of colors to use in the table. If the color list has no colors we add white.
	fill_color_list = plot_dict['table_color_list']
	if not fill_color_list:
		fill_color_list = ['#ffffff']

	table_pos_list = []

	# Defining a dataframe with the first column
	data_func = {col_names_list[0]: col_values_list[0]}
	df_func = pd.DataFrame(data=data_func)

	# Go over the lists loaded into the function and append them to the dataframe.
	for a in range(1, len(col_names_list)):
		df_func[col_names_list[a]] = col_values_list[a]

	# Check if the number of samples given by user is the same as given to the function. Deviations can occur e.g. in
	# global fittings where the code is generating extra traces.
	if len(ID_list) == len(col_values_list[0]):
		# Code for extracting the user defined table coordinates i.e. what line in the table the data should be inputted on.
		# First loop iterates through ID_list i.e. the samples in the data sheet and extracts the flags separated by semi colon.
		# The code goes over each flag and if table_pos= flag is present it extracts the integer following the equal sign.

		for b in range(len(ID_list)):
			func_var = python_misc[b].split(';')

			# If user has specified table_pos in the python_misc column the entry will be extracted into dummy list.
			func_list = []
			for c in range(len(func_var)):
				if 'table_pos=' in func_var[c]:
					func_list.append(func_var[c])

			# If no table position has been specified the code will just use the row index from excel sheet.
			if len(func_list) == 0:
				table_pos_list.append(b + 1)
			else:
				func_coord = func_var[c].split('=')[1]
				func_coord = int(func_coord)
				table_pos_list.append(func_coord)

		# Add the table position list to the dataframe
		df_func['table_coords'] = table_pos_list

		# Sort the dataframe by the table coordiates to get the desired order
		df_func2 = df_func.sort_values(by='table_coords', ascending=True)
		del df_func2['table_coords']

		if param_dict['table_coloring'] == 'yes':

			figure.add_trace(
				go.Table(
					header=dict(values=list(df_func2.columns), align='left'),
					cells=dict(values=df_func2.transpose().values.tolist(), align='left',
							   height=50, fill_color=[fill_color_list]*len(col_names_list)))
			)

		else:

			figure.add_trace(
				go.Table(
					header=dict(values=list(df_func2.columns), align='left'),
					cells=dict(values=df_func2.transpose().values.tolist(), align='left', height=50),
					meta='table')
			)

	else:

		if param_dict['table_coloring'] == 'yes':

			figure.add_trace(
				go.Table(
					header=dict(values=list(df_func.columns), align='left'),
					cells=dict(values=df_func.transpose().values.tolist(), align='left',
							   height=50, fill_color=[fill_color_list]*len(col_names_list)))
			)

		else:

			figure.add_trace(
				go.Table(
					header=dict(values=list(df_func.columns), align='left'),
					cells=dict(values=df_func.transpose().values.tolist(), align='left', height=50))
			)

def get_ax_id(subplot_id, user_input_dict):

	# The function gets the internal ids that is used for navigating the plot layout

	# We need to keep track of the number of secondary y axis have been added because this affects the axis ids.
	# The index is being added both for the interactive plot and the static plot so we need to remove the duplicates
	sec_y_count = user_input_dict['second_y_traces']
	sec_y_count = list(dict.fromkeys(sec_y_count))
	sec_y_count = len(sec_y_count)

	#print(sec_y_count)

	if subplot_id == 1 and sec_y_count == 0:
		x_ax_id = ''
		y_ax_id = ''
	else:
		x_ax_id = str(subplot_id)
		y_ax_id = str((subplot_id*2)-1)
		#y_ax_id = str(subplot_id + sec_y_count)

	return x_ax_id, y_ax_id
