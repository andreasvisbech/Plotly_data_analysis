'''
Author: Andreas Visbech Madsen
'''

__version__ = "3.2.0"

# Load modules
import argparse
import sys

# Load functions from related scripts
from plot_it_scripts.advanced_box import create_advanced_box, default_param_dict
from plot_it_scripts.main_functions import *

# Initialize parser
parser = argparse.ArgumentParser()
parser.add_argument("-s", "--plot_type", help="plot type", required=True)
parser.add_argument("-i", "--input_file", help="Excel file", required=True)
parser.add_argument("-out", "--output", help="y/n to output", action='store_true')
parser.add_argument('-v', '--version', action='version', version=__version__)
parser.add_argument("-plot", "--plotting", help="y/n to output", action='store_true')
parser.add_argument("-log_out", "--log_output", help="y/n to log file output", action='store_true')
parser.add_argument("-log_in", "--log_input", help="y/n to log file output")
parser.add_argument("-advanced", "--advanced_option_box", help="y/n to advanced option box", action='store_true')
args = parser.parse_args()

# Figure out the type of analysis
analysis_type = get_analysis_type(args.plot_type)

if args.advanced_option_box:
	param_dict = create_advanced_box(analysis_type)

else:
	param_dict = default_param_dict()

# If the user has supplied a log file using the log_in flag then the param_dict will be overwritten with the user input
if args.log_input:
	log_file = read_log_file_in(args.log_input, param_dict)

# Define globals
global used_graph_names

# Load data from excel sheet and prep it for analysis
print('Loading the data')
df = pd.concat(pd.read_excel(args.input_file, sheet_name=None), ignore_index=True)
print('Cleaning the data')
df = data_clean(df)

# Getting user inputs from excel sheet and separating it into appropriate variables
user_input_dict = load_user_input(df)
ID_list = user_input_dict['ID_list']

# Defining a master dictionary for storing data
master_dict = define_master_dict()
master_dict['columns'] = df.columns

# Defining figure for making the plots in. 
fig = create_subplot_function(user_input_dict)
plot_fig = create_subplot_function(user_input_dict)

# Define color list for coloring graphs
color_list = define_color_list(user_input_dict, param_dict)

# Defining plotting dict
plot_dict = define_plot_dict()
plot_dict['figure'] = fig
plot_dict['plot_figure'] = plot_fig
plot_dict['color_list'] = color_list

if analysis_type == 'scatter':
	from plot_it_scripts.scatter_module_v2 import *

	# Go over each sample in the excel sheet
	for i in range(len(ID_list)):

		print('Analysing data: ' + str(ID_list[i]))

		plot_dict['graph_names'].append(ID_list[i])

		x_id = 'x' + str(i + 1)
		y_id = 'y' + str(i + 1)
		ignore_id = 'ignore' + str(i + 1)
		graph_name = ID_list[i]

		error_id = 'error' + str(i + 1)
		plot_dict['error_id'] = error_id

		# Updating the color counter to ensure the graphs are different colors.
		plot_dict['color_count'] = color_selector(plot_dict, graph_name, used_graph_names)

		# Check if the "ignore column" exists i.e. if user has specified dropping specific values of data.
		if ignore_id in df.columns:
			df = df[df[ignore_id].fillna(' ').str.isalpha() == False]

		# Check if user has supplied error values
		if error_id in df.columns:
			plot_dict['errors'] = df[error_id][pd.to_numeric(df[error_id], errors='coerce').notnull()]

		# Getting x values, y values and a non-redundant list of x values
		xs, ys, unique_x = scatter_data_slice(df, i, x_id, y_id, user_input_dict)

		# Check if the data should just be plotted or if user is trying to fit.
		# This is checked by seing if user has inputed fitting interval and fitting model
		if user_input_dict['fit_intervals'][i] != '0;0' and user_input_dict['fit_models'][i] != 'None':
			scatter_plot_with_fit(i, user_input_dict, plot_dict, param_dict, master_dict, xs, ys, unique_x)

		else:
			scatter_plot_without_fit(master_dict, i, user_input_dict, plot_dict, xs, ys, param_dict)

		# Updating the color counter to ensure the graphs are different colors.
		plot_dict['color_count'] = color_selector(plot_dict, graph_name, used_graph_names)

		# Adding table to the interactive plot
		table_plot(plot_dict, [
			'Samples',
			'Sample notes',
			'Model list',
			'Fitted KD/EC50',
			'Goodness of fit (R<sup>2</sup>)',
			'RMSE value (beta)',
			'Fit parameters'], [
					   master_dict['ID_list_new'],
					   master_dict['notes_list'],
					   master_dict['model_list'],
					   master_dict['KD_fit'],
					   master_dict['R_square'],
					   master_dict['RMSE'],
					   master_dict['fit_parameters']], user_input_dict, param_dict)

	# Adding interactive buttons to the plotly plot
	plotly_buttons(plot_dict)

elif analysis_type == 'fplc':

	sys.path.append('./plot_it_scripts/')

	from plot_it_scripts.akta_module import *

	# Go over each sample in the excel sheet
	for i in range(len(ID_list)):

		print('Analysing data: ' + str(ID_list[i]))

		plot_dict['graph_names'].append(ID_list[i])
		#master_dict['ID_list_new'].append(ID_list[i])
		#master_dict['notes_list'].append(user_input_dict['sample_notes'][i])

		x_id = 'x' + str(i + 1)
		y_id = 'y' + str(i + 1)
		graph_name = ID_list[i]

		# Updating the color counter to ensure the graphs are different colors.
		plot_dict['color_count'] = color_selector(plot_dict, graph_name, used_graph_names)
		
		# Extracting the raw x and y values from the excel sheet
		xs = df[x_id][pd.to_numeric(df[x_id], errors='coerce').notnull()]
		ys = df[y_id][pd.to_numeric(df[y_id], errors='coerce').notnull()]

		# Slicing the data so only data within the specified data interval is included.
		xs, ys = akta_data_slice(x_id, y_id, i, xs, ys, user_input_dict)

		# If specified by user all the negative y values will be replaced by zeros.
		# Can be relevant since the negative values negatively impact the area calculations
		if param_dict['AKTA neg val handle'] == 'yes':
			ys[ys < 0] = 0

		akta_main_func(df, xs, ys, i, x_id, y_id, param_dict, master_dict, user_input_dict, plot_dict)

	# Adding table to the interactive plot
	table_plot(plot_dict, [
		'Sample ID',
		'Sample notes',
		'Peak ID',
		'Total sample area',
		'Total baseline area',
		'Retention time/volume (beta)',
		'Peak area (raw)',
		'Peak area (%)',
		'Fraction yield [mg]',
		'Culture yield [ug/mL]'],
		[
			master_dict['ID_list_new'],
			master_dict['notes_list'],
			master_dict['peak_id'],
			master_dict['sample_areas_tot'],
			master_dict['baseline_area_tot'],
			master_dict['fraction_retentions'],
			master_dict['fraction_calculation'],
			master_dict['fraction_AUC_of_total'],
			master_dict['fraction_yield'],
			master_dict['culture_yield']], user_input_dict, param_dict)

	# Adding interactive buttons to the plotly plot
	plotly_buttons(plot_dict)

	# Creating an output file with the raw values used for the plotting
	pd_out = write_output_table(['ID_list', ], [master_dict['ID_list']])

elif analysis_type == 'fida':

	from plot_it_scripts.scatter_module_v2 import *
	from plot_it_scripts.data_manipulation import *

	# Go over each sample in the excel sheet
	for i in range(len(ID_list)):

		print('Analysing data: ' + str(ID_list[i]))

		plot_dict['graph_names'].append(ID_list[i])
		#master_dict['notes_list'].append(user_input_dict['sample_notes'][i])

		x_id = 'x' + str(i + 1)
		y_id = 'y' + str(i + 1)
		ignore_id = 'ignore' + str(i+1)
		graph_name = ID_list[i]

		# Updating the color counter to ensure the graphs are different colors.
		plot_dict['color_count'] = color_selector(plot_dict, graph_name, used_graph_names)

		# Check if the "ignore column" exists i.e. if user has specified dropping specific values of data.
		if ignore_id in df.columns:
			df = df[df[ignore_id].fillna(' ').str.isalpha() == False]

		# Extract concentration float values and float apparent Rh values from the text FIDA output
		# Extract the values into lists and make into dataframe
		xs, ys = fida_text2floats(x_id, y_id, df)
		df_new = {x_id: xs, y_id: ys}
		df_new = pd.DataFrame(df_new)

		# Getting x values, y values and a non-redundant list of x values
		xs, ys, unique_x = scatter_data_slice(df_new, i, x_id, y_id, user_input_dict)

		# Check if the data should just be plotted or if user is trying to fit.
		# This is checked by seing if user has inputed fitting interval and fitting model
		if user_input_dict['fit_intervals'][i] != '0;0' and user_input_dict['fit_models'][i] != 'None':
			scatter_plot_with_fit(i, user_input_dict, plot_dict, param_dict, master_dict, xs, ys, unique_x)

		else:
			scatter_plot_without_fit(master_dict, i, user_input_dict, plot_dict, xs, ys, param_dict)

	# Adding table to the interactive plot
	table_plot(plot_dict,[
		'Samples',
		'Sample notes',
		'Model list',
		'Fitted KD/EC50',
		'Goodness of fit (R<sup>2</sup>)',
		'RMSE value (beta)',
		'Fit parameters'],	[
		master_dict['ID_list_new'],
		master_dict['notes_list'],
		master_dict['model_list'],
		master_dict['KD_fit'],
		master_dict['R_square'],
		master_dict['RMSE'],
		master_dict['fit_parameters']], user_input_dict, param_dict)

	# Adding interactive buttons to the plotly plot
	plotly_buttons(plot_dict)

elif analysis_type == 'bar':

	from plot_it_scripts.bar_module import *
	#from plot_it_scripts.bar_module import *

	for i in range(len(ID_list)):
		print('Analysing data: ' + str(ID_list[i]))

		plot_dict['graph_names'].append(ID_list[i])
		master_dict['notes_list'].append(user_input_dict['sample_notes'][i])

		x_id = 'x' + str(i + 1)
		y_id = 'y' + str(i + 1)
		error_id = 'error' + str(i + 1)
		graph_name = ID_list[i]

		# Updating the color counter to ensure the graphs are different colors.
		plot_dict['color_count'] = color_selector(plot_dict, graph_name, used_graph_names)

		bar_main(df, x_id, y_id, error_id, i, plot_dict, user_input_dict, param_dict, master_dict)

elif analysis_type == 'octet':

	from plot_it_scripts.octet_module import *

	# Re-define a specific plot for the octet.
	fig = create_subplot_octet()
	#plot_fig = create_subplot_octet()
	plot_dict['figure'] = fig
	plot_dict['plot_figure'] = plot_fig

	# If no parameter log file has been supplied the user is prompted for input.
	if args.log_input:
		None
	else:
		param_dict = octet_box(param_dict)

	# Go over each sample in the excel sheet
	for i in range(len(ID_list)):
		print('Analysing data: ' + str(ID_list[i]))

		plot_dict['graph_names'].append(ID_list[i])

		graph_name = ID_list[i]

		# Updating the color counter to ensure the graphs are different colors.
		plot_dict['color_count'] = color_selector(plot_dict, graph_name, used_graph_names)

		octet_main(i, user_input_dict, plot_dict, param_dict, master_dict, df)

	# Adding table to the interactive plot
	table_plot(plot_dict, [
		'Samples',
		'Sensor',
		'Concentrations (' + param_dict['octet_conc_unit'] + ')',
		'ka (M<sup>-1</sup>s<sup>-1</sup>)',
	    # 'ka stderr',
		'kd (s<sup>-1</sup>)',
		# 'kd stderr',
		'Kinetics KD (' + param_dict['octet_conc_unit'] + ')',
		'Kinetic KD error',
		'R2 full',
		'Steady-state KD (' + param_dict['octet_conc_unit'] + ')',
		'Steady-state R2'], [
		master_dict['ID_list_new'],
		master_dict['octet_sensors'],
		master_dict['octet_sensor_conc'],
		master_dict['octet_ka'],
		# master_dict['octet_ka_err'],
		master_dict['octet_kd'],
		# master_dict['octet_kd_err'],
		master_dict['octet_kinetic_KD'],
		master_dict['octet_kinetic_KD_err'],
		master_dict['octet_R2_full'],
	    master_dict['octet_KD_SS'],
	    master_dict['octet_R2_SS']], user_input_dict, param_dict)

	# Adding interactive buttons to the plotly plot
	plotly_buttons(plot_dict)

elif analysis_type == 'panta':

	from plot_it_scripts.panta_module import *

	# Go over each sample in the excel sheet
	for i in range(len(ID_list)):
		print('Analysing data: ' + str(ID_list[i]))

		plot_dict['graph_names'].append(ID_list[i])
		master_dict['notes_list'].append(user_input_dict['sample_notes'][i])

		x_id = 'x' + str(i + 1)
		y_id = 'y' + str(i + 1)
		graph_name = ID_list[i]

		# Updating the color counter to ensure the graphs are different colors.
		plot_dict['color_count'] = color_selector(plot_dict, graph_name, used_graph_names)

		panta_main(df, x_id, y_id, i, plot_dict, user_input_dict, param_dict, master_dict)

	# Adding table to the interactive plot
	table_plot(plot_dict, ['Samples',
						   'Sample notes',
						   'Onset (beta)',
						   'Inflection point (beta)',
						   'Data vertex points (beta)'],
			   [user_input_dict['ID_list'],
				master_dict['notes_list'],
				master_dict['peak_onset'],
				master_dict['inflection_points'],
				master_dict['vertex_max']], user_input_dict, param_dict)

	# Adding interactive buttons to the plotly plot
	plotly_buttons(plot_dict)

if args.log_output == True:
	output_file_name = str(args.input_file[:len(args.input_file) - 5])
	log_file_out(param_dict, output_file_name)

# If specified by user output the plotting figure as an svg file    
if args.plotting == True:
	output_file_name = str(args.input_file[:len(args.input_file) - 5])

	write_plot_fig_out(output_file_name, plot_fig, param_dict)

# Writing output file
output_file_name = str(args.input_file[:len(args.input_file) - 5]) + '_Output.html'
html_plot_out(fig, output_file_name)

quote("quotes.csv")