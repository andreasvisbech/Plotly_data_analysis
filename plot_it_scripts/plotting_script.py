# Import relevant modules
import plotly.graph_objects as go
import pandas as pd
import matplotlib.colors
import numpy as np

global used_graph_names
used_graph_names = []

global used_color_IDs
used_colors = []


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
	# global used_graph_names, color_count_global

	# Define template for hover label
	my_hover_template = graph_name + '<extra></extra>' + '<br>x: %{x}' + '<br>y: %{y}<br>' + 'Note: ' + \
						user_input_dict['sample_notes'][i]


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
					line=dict(width=param_dict['graph width'], color=color_list[color_count])),
				row=subplot_row, col=subplot_col
			)

		elif marker.lower() in ['dot', 'dots', 'marker', 'markers']:
			figure.add_trace(
				go.Scatter(
					name=graph_name, x=x_val, y=y_val, mode='markers', legendgroup=graph_name, showlegend=legend_show,
					hovertemplate=my_hover_template,
					marker=dict(size=param_dict['marker_size'], color=color_list[color_count])),
				row=subplot_row, col=subplot_col
			)

	elif comment == 'Scatter_error':

		rgb_color = matplotlib.colors.to_rgb(color_list[color_count])
		rgb_color = 'rgba(' + str(rgb_color[0]) + ',' + str(rgb_color[1]) + ',' + str(rgb_color[2]) + ', 0.3)'

		y_val_upper = np.add(y_val, std_dev)
		y_val_lower = np.subtract(y_val, std_dev)

		if marker.lower() in ['lines', 'line']:

			if param_dict['error_marker'] == 'Error bands':

				figure.add_trace(
					go.Scatter(
						name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show,
						line=dict(width=param_dict['graph width'], color=color_list[color_count])),
					row=subplot_row, col=subplot_col
				)

				figure.add_trace(
					go.Scatter(
						name=graph_name + '_upper', x=x_val, y=y_val_upper, mode='lines', legendgroup=graph_name,
						showlegend=False, line=dict(width=0, color=color_list[color_count])),
					row=subplot_row, col=subplot_col
				)

				figure.add_trace(
					go.Scatter(
						name=graph_name + '_lower', x=x_val, y=y_val_lower, mode='lines', legendgroup=graph_name,
						showlegend=False, fillcolor=rgb_color, fill='tonexty',
						line=dict(width=0, color=color_list[color_count])),
					row=subplot_row, col=subplot_col
				)


			else:
				figure.add_trace(
					go.Scatter(
						name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show,
						hovertemplate=my_hover_template,
						error_y=dict(type='data', array=std_dev, visible=True),
						line=dict(width=param_dict['graph width'], color=color_list[color_count])),
					row=subplot_row, col=subplot_col
				)

		elif marker.lower() in ['dot', 'dots', 'marker', 'markers']:

			if param_dict['error_marker'] == 'Error bands':

				figure.add_trace(
					go.Scatter(
						name=graph_name, x=x_val, y=y_val, mode='markers', legendgroup=graph_name,
						showlegend=legend_show,
						hovertemplate=my_hover_template,
						marker=dict(size=param_dict['marker_size'], color=color_list[color_count])),
					row=subplot_row, col=subplot_col
				)

				figure.add_trace(
					go.Scatter(
						name=graph_name + '_upper', x=x_val, y=y_val_upper, mode='lines', legendgroup=graph_name,
						showlegend=False,
						hovertemplate=my_hover_template,
						line=dict(width=0, color=color_list[color_count])),
					row=subplot_row, col=subplot_col
				)

				figure.add_trace(
					go.Scatter(
						name=graph_name + '_lower', x=x_val, y=y_val_lower, mode='lines', legendgroup=graph_name,
						showlegend=False, fillcolor=rgb_color, fill='tonexty',
						hovertemplate=my_hover_template,
						line=dict(width=0, color=color_list[color_count])),
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
						marker=dict(size=param_dict['marker_size'], color=color_list[color_count])),
					row=subplot_row, col=subplot_col
				)

	elif comment == 'scatter_residuals':

		figure.add_trace(
			go.Scatter(
				name=graph_name, x=x_val, y=y_val, mode='markers', marker_symbol='x',
				legendgroup=graph_name, showlegend=legend_show,
				hovertemplate=my_hover_template,
				marker=dict(size=param_dict['marker_size'], color=color_list[color_count])),
			row=subplot_row, col=subplot_col
		)

	elif comment == 'Bar':
		figure.add_trace(
			go.Bar(name=graph_name, x=x_val, y=y_val, legendgroup=graph_name, showlegend=legend_show,
				   hovertemplate=my_hover_template,
				   marker=dict(color=color_list[color_count])),
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
				marker=dict(color=color_list[color_count])),
			row=subplot_row, col=subplot_col
		)

	elif comment == 'AKTA_baseline':
		figure.add_trace(
			go.Scatter(name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show,
					   hovertemplate=my_hover_template,
					   line=dict(width=param_dict['graph width'], color='rgb(255,0,0)')),
			row=subplot_row, col=subplot_col
		)

	elif comment == 'AKTA_fraction':
		figure.add_trace(
			go.Scatter(
				name=graph_name, x=x_val, y=y_val, mode='lines', fill='tozeroy', legendgroup=graph_name,
				showlegend=legend_show,
				hovertemplate=my_hover_template,
				line=dict(width=param_dict['graph width'], color=color_list[color_count])),
			row=subplot_row, col=subplot_col
		)

	# Determining the "number id" for the individual subplot. This can be used for targeting customization
	# such as axis labels to this subplot.
	subplot_id = ((subplot_row - 1) * user_input_dict['subplot_col_count']) + subplot_col
	if subplot_id == 1:
		ax_id = ''
	else:
		ax_id = str(subplot_id)

	figure['layout']['xaxis' + ax_id]['title'] = x_title
	figure['layout']['yaxis' + ax_id]['title'] = y_title

	#figure.for_each_xaxis(lambda axis: axis.title.update(font=dict(size=param_dict['xaxis_title_font_size'])))
	#figure.for_each_yaxis(lambda axis: axis.title.update(font=dict(size=param_dict['yaxis_title_font_size'])))

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

	if user_input_dict['python_misc'][i] != 'None':
		plot_customize(figure, user_input_dict['python_misc'][i], ax_id)

def plot_customize(figure, flags, ax_id):
	in_list = flags.split(';')

	for a1 in range(len(in_list)):
		if in_list[a1] == 'logx':
			figure['layout']['xaxis' + ax_id]['type'] = 'log'
			figure['layout']['xaxis' + ax_id]['dtick'] = 1

		elif in_list[a1] == 'logy':
			figure['layout']['yaxis' + ax_id]['type'] = 'log'
			figure['layout']['yaxis' + ax_id]['dtick'] = 1


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
				pad={"r": 10, "t": 10}, showactive=True, x=0.11, xanchor="left", y=1.1, yanchor="top"), ])

def table_plot(plot_dict, col_names_list, col_values_list, user_input_dict):
	"""
    The function is intended for making the table going into the html output.
    Variable "figure" specifies the figure that the table are added to.
    Variable "col_names_list" specifies a list the the column headers that should go into the table.
    Variable "col_values_list" is a list of lists. Each list contains data for an individual column

    :param figure:
    :param col_names_list:
    :param col_values_list:
    :param ID_list:
    :return:
    """

	# Read local variables into function
	figure = plot_dict['figure']
	ID_list = user_input_dict['ID_list']
	python_misc = user_input_dict['python_misc']

	fill_color_list = plot_dict['table_color_list']

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

		if len(fill_color_list) > 0:

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

		if len(fill_color_list) > 0:

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

