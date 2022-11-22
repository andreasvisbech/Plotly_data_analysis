# Import modules
import pandas as pd
import matplotlib.colors
from plotly.subplots import make_subplots
import random


# Define functions
def data_clean(data):
	"""
	Split columns where semicolon separator is used to indicate that the column
	should be used in more than one data analysis

	:param data:
	:return:
	"""
	for col in data.columns:
		if str(col).count(';'):  # Check if semicolon separators are present in column names

			col_data = data[col]  # Store the column data

			col_list = col.split(';')
			for i in range(len(col_list)):
				new_frame = pd.DataFrame({str(col_list[i]): col_data})
				data = pd.concat([data, new_frame], axis=1)

	# Check for duplicate columns
	for col in data.columns:
		for a in range(20):
			col_check = str(col) + '.' + str(int(a))
			#col_check = str(col + '.' + str(int(a)))
			#print(data.columns.tolist().count(col_check))
			if data.columns.tolist().count(col_check) > 0 and ('x' or 'y') in col_check:
				raise ValueError('Column ' + str(col) + ' appears more than once!')

	return data


def load_user_input(df):
	# The function collects the user input from the excel sheet and outputs it into a dictionary

	out_dict = {}

	ID_list = df['IDs'].tolist()
	ID_list = [x for x in ID_list if str(x) != 'nan']
	out_dict['ID_list'] = ID_list

	data_interval = df['Data_interval'].fillna(0).replace(',',
														  '.').tolist()  # Get the data intervals for slicing the data to only analyse sections of full dataset.
	data_interval = data_interval[0:len(ID_list)]
	out_dict['data_interval'] = data_interval

	plot_markers = df['Plot_markers'].fillna('line').tolist()
	out_dict['plot_markers'] = plot_markers

	subplot_coord = df['Sub_plot'].fillna('1;1').tolist()
	out_dict['subplot_coord'] = subplot_coord

	sample_notes = df['Notes'].fillna('None').tolist()
	out_dict['sample_notes'] = sample_notes

	python_misc = df['Python_misc'].fillna('None').tolist()
	out_dict['python_misc'] = python_misc

	# Getting the axis titles for each sample
	axis_title = df['Axis_titles'].fillna(';').tolist()
	x_titles = []
	y_titles = []
	for id_b in range(len(axis_title)):
		x_titles.append(axis_title[id_b].split(';')[0])
		y_titles.append(axis_title[id_b].split(';')[1])

	out_dict['x_titles'] = x_titles
	out_dict['y_titles'] = y_titles

	# Get the data intervals for slicing the data to only analyse sections of full dataset.
	data_interval = df['Data_interval'].fillna(0).replace(',', '.').tolist()
	data_interval = data_interval[0:len(ID_list)]
	out_dict['data_interval'] = data_interval

	plot_markers = df['Plot_markers'].fillna('line').tolist()
	out_dict['plot_markers'] = plot_markers

	sample_notes = df['Notes'].fillna('None').tolist()
	out_dict['sample_notes'] = sample_notes

	fit_models = df['Fit_model'].fillna('None').tolist()
	out_dict['fit_models'] = fit_models

	fit_modes = df['Fit_approach'].fillna('Local').tolist()
	out_dict['fit_modes'] = fit_modes

	fit_intervals = df['Fitting_interval'].fillna('0;0').tolist()
	out_dict['fit_intervals'] = fit_intervals

	# Getting user info on the subplotting setup
	subplot_row = []
	subplot_col = []
	subplot_coord = df['Sub_plot'].fillna('1;1').tolist()
	for id_a in range(len(subplot_coord)):
		# The first subplot coordinate is appended to subplot row list as integer
		subplot_row.append(int(subplot_coord[id_a].split(';')[0]))
		# The first subplot coordinate is appended to subplot row list as integer
		subplot_col.append(int(subplot_coord[id_a].split(';')[1]))
	out_dict['subplot_row'] = subplot_row
	out_dict['subplot_col'] = subplot_col
	subplot_row_count = max(subplot_row)
	subplot_col_count = max(subplot_col)
	out_dict['subplot_row_count'] = subplot_row_count
	out_dict['subplot_col_count'] = subplot_col_count

	# Get listed extinction coefficient from excel sheet and replace empty values with zero. Also make sure the the right decimal seperator is used in case people use e.g. Danish excel
	ext_coeff = df['AKTA_extinc_coeff'].fillna(0).replace(',', '.').tolist()
	# ext_coeff = ext_coeff[0:len(ID_list)]  # Restrict the list so it has the same length as ID list. Basically just removing nan values
	out_dict['ext_coeff'] = ext_coeff

	# Get the volumes loaded. This is the total volume of supernatant loaded during experiments. Replace empty values with zero
	volume_load = df['AKTA_volume_load'].fillna(0).replace(',', '.').tolist()
	# volume_load = volume_load[0:len(ID_list)]  # Restrict the list so it has the same length as ID list
	out_dict['akta volume loaded'] = volume_load

	return out_dict


def define_master_dict():
	master_dict = {}

	master_dict['ID_list'] = []
	master_dict['ID_list_new'] = []
	master_dict['notes_list'] = []
	master_dict['model_list'] = []
	master_dict['fit_parameters'] = []
	master_dict['R_square'] = []
	master_dict['func_R2_list'] = []
	master_dict['RMSE'] = []
	master_dict['chi_square'] = []
	master_dict['Fisher_test'] = []
	master_dict['KD_fit'] = []
	master_dict['func_KD_list'] = []
	master_dict['peak_id'] = []
	master_dict['sample_areas_tot'] = []
	master_dict['baseline_area_tot'] = []
	master_dict['fraction_retentions'] = []
	master_dict['fraction_areas'] = []
	master_dict['fraction_baseline'] = []
	master_dict['fraction_calculation'] = []
	master_dict['fraction_concentrations'] = []
	master_dict['fraction_yield'] = []
	master_dict['culture_yield'] = []
	master_dict['peak_onset'] = []
	master_dict['inflection_points'] = []
	master_dict['vertex_max'] = []
	master_dict['octet_sensors'] = []
	master_dict['octet_sensor_conc'] = []
	master_dict['octet_ka'] = []
	master_dict['octet_ka_err'] = []
	master_dict['octet_kd'] = []
	master_dict['octet_kd_err'] = []
	master_dict['octet_kinetic_KD'] = []
	master_dict['octet_kinetic_KD_err'] = []
	master_dict['octet_R2_full'] = []
	master_dict['octet_KD_SS'] = []
	master_dict['octet_R2_SS'] = []

	return master_dict


def define_plot_dict():
	plot_dict = {}

	plot_dict['figure'] = ''
	plot_dict['plot_figure'] = ''
	plot_dict['graph_names'] = []
	plot_dict['color_count'] = -1
	plot_dict['table_color_list'] = []

	return plot_dict


def define_color_list(user_input_dict, param_dict):

	color_scheme = param_dict['color_palette']

	if color_scheme == 'Default (10 color)':
		# Loading a list with colors for plotting. The colors come from https://plotly.com/python/discrete-color/.
		color_list = ['#1F77B4', '#FF7F0E', '#2CA02C', '#9467BD', '#FECB52', '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22',
					  '#17BECF']

	elif color_scheme == '15 color':
		color_list = ['#9352a8', '#f7760b', '#d51c3c', '#c8b18b', '#23eaa5', '#f483cd', '#276cbd',
					  '#2ca02c', '#ebdd21', '#8b1c0e', '#a7dc26', '#673f0b', '#00FFFF',
					 '#8a8489', '#20B2AA']

	elif color_scheme == '18 color':
		color_list =['#9352a8', '#f7760b', '#d51c3c', '#c8b18b', '#23eaa5', '#f483cd', '#276cbd', '#f59080',
					 '#61419c', '#b83773', '#ebdd21', '#8b1c0e', '#a7dc26', '#673f0b', '#e83b1b', '#495b22',
					 '#8a8489', '#99c6f9']

	elif color_scheme == '4 color (color blind safe':
		color_list = ['#D41159', '#40B0A6', '#D35FB7', '#52D1EA']

	# If the number of samples equals the length of the color list the script will append one extra color to avoid the same
	# color comparisons in subplots.
	if len(user_input_dict['ID_list']) % len(color_list) == 0:
		color_list.append('#B6E880')

	return color_list


def color_selector(plot_dict, graph_name, used_graph_names):
	"""
	The function is used to set the color-determining index used in the plotting function.
	If more graphs are plotted than there are colors in the color_list then the counting will start over.

	:param plot_dict:
	:return:
	"""

	color_count = plot_dict['color_count']
	color_list = plot_dict['color_list']

	#if used_graph_names.count(graph_name) > 1:
	#	None
	#else:
	if used_graph_names.count(graph_name) < 1:
		if color_count < len(color_list) - 1:
			color_count += 1
		else:
			color_count = 0

	return color_count


def html_plot_out(figure, output_file_name):
	figure.write_html(output_file_name)


def write_output_table(name_list, data_list):
	# Making an empty dataframe
	pd_out = pd.DataFrame({})

	# Make sure lengths of the two lists are the same
	if len(name_list) == len(data_list):

		for a in range(len(name_list)):
			pd_out = pd.concat([pd_out, pd.DataFrame({name_list[a]: data_list[a]}).reset_index(drop=True)],
							   ignore_index=False, axis=1)

	else:
		print('Output table could not be made')

	return pd_out


def table_color_list_manager(color,list):
	# Add the appropriate color to the table coloring list

	rgb_color = matplotlib.colors.to_rgb(color)

	rgb_color = 'rgba(' + \
				str(rgb_color[0] * 255) + ',' + \
				str(rgb_color[1] * 255) + ',' + \
				str(rgb_color[2] * 255) + ', 0.15)'

	list.append(rgb_color)


def create_subplot_function(ID_list, user_input_dict):

	figure = make_subplots(
		rows=user_input_dict['subplot_row_count'],
		cols=user_input_dict['subplot_col_count'],
		vertical_spacing=0.07,
		horizontal_spacing=0.07)

	# Positioning the subplot titles
	#for c in figure['layout']['annotations']:
	#	c['xref'] = 'x' + str(ax_id)
	#	c['x'] = 0.1
		#c['yref'] = 'y' + str(ax_id)
		#c['y'] = 0.1

	return figure

def create_subplot_octet():

	figure = make_subplots(rows=4,
						   cols=4,
						   vertical_spacing=0.1,
						   horizontal_spacing=0.04,
						   specs=[[{'rowspan': 2, 'colspan': 3}, None, None, {'rowspan': 2}],
								  [{}, {}, {}, {}],
								  [{'colspan': 2}, None, {'rowspan': 2}, {'rowspan': 2}],
								  [{'colspan': 2}, None, None, None]])

	return figure


def log_file_out(param_dict, file_name):

	log_pd = pd.DataFrame.from_dict(param_dict, orient = 'index')

	log_pd.to_csv(file_name + '_log.txt', sep=',')

def quote(file):
	df = pd.read_csv(file)

	author_list = df['Author'].tolist()
	quote_list = df['Quote'].tolist()

	number = random. randint(0, len(author_list))

	print(' ')
	print('Thank you for using the script! And remember: "' + str(quote_list[number]) + '"' + ' - ' + str(author_list[number]))
	print(' ')


def write_plot_fig_out(name, plot_fig, param_dict):

	if param_dict['plot_fig_width'] != 'N/A' and param_dict['plot_fig_height'] != 'N/A':

		plot_fig_width = float(param_dict['plot_fig_width'])
		plot_fig_height = float(param_dict['plot_fig_height'])

		plot_fig.update_layout(
			autosize=False,
			width=plot_fig_width,
			height=plot_fig_height)

	# If specified by user manually set the x axis range,
	if param_dict['plot_fig_xmin'] != 'N/A' and param_dict['plot_fig_xmax'] != 'N/A':

		x_min = float(param_dict['plot_fig_xmin'])
		x_max = float(param_dict['plot_fig_xmax'])

		plot_fig.update_xaxes(range=[x_min, x_max])

	# If specified by user manually set the y axis range.
	if param_dict['plot_fig_ymin'] != 'N/A' and param_dict['plot_fig_ymax'] != 'N/A':

		y_min = float(param_dict['plot_fig_ymin'])
		y_max = float(param_dict['plot_fig_ymax'])

		plot_fig.update_yaxes(range=[y_min, y_max])

	plot_fig.write_image(name + '_Output.svg')

def get_analysis_type(input):

	if input in ['Scatter', 'scatter']:
		analysis_type = 'scatter'
	elif input in ['AKTA', 'akta', 'Akta', 'FPLC', 'fplc']:
		analysis_type = 'fplc'
	elif input in ['FIDA', 'fida', 'Fida']:
		analysis_type = 'fida'
	elif input in ['panta', 'Panta', 'PANTA']:
		analysis_type = 'panta'
	elif input in ['Octet', 'OCTET', 'octet']:
		analysis_type = 'octet'
	elif input in ['Bar', 'BAR', 'bar']:
		analysis_type = 'bar'


	return analysis_type