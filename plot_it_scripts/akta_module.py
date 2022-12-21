# Import modules
from sklearn.metrics import auc
import peakutils

# Import functions from other scripts
from plot_it_scripts.plotting_script import *
from plot_it_scripts.main_functions import *
from plot_it_scripts.data_manipulation import *



def akta_data_slice(x_id, y_id, sample_idx, xs, ys, user_input_dict):
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


def akta_main_func(df, xs, ys, sample_idx, x_id, y_id, param_dict, master_dict, user_input_dict, plot_dict):

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

	# user specified values loaded
	ext_coeff = user_input_dict['ext_coeff']
	volume_load = user_input_dict['akta volume loaded']

	# Check the user input for the baseline. Single number input will give flat baseline.
	# Two number input will be used as boundaries for baseline
	if str(df['AKTA_baseline'][sample_idx]).count(';') == 0:
		baseline_coord = [float(str(df['AKTA_baseline'][sample_idx]).replace(',', '.'))]
	elif df['AKTA_baseline'][sample_idx].count(';') > 0:
		baseline_coord = list(df['AKTA_baseline'][sample_idx].replace(',', '.').split(';'))

	baseline_values = akta_baseline(baseline_coord, xs, ys, param_dict['AKTA baseline mode'], param_dict)

	# Calculate the total AUC from entire sample
	AUC_sample_tot = auc(xs, ys)
	#AUC_sample_tot = "{:.3f}".format(AUC_sample_tot)
	#master_dict['sample_areas_tot'].append(AUC_sample_tot)

	# Calculate AUC for the baseline
	AUC_baseline_tot = auc(baseline_values[0], baseline_values[1])
	#AUC_baseline_tot = "{:.3f}".format(AUC_baseline_tot)
	#master_dict['baseline_area_tot'].append(AUC_baseline_tot)

	# Plotting the data values into interactive plot and the stativ "plotting plot"
	plot_func(figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col,
			  'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
	plot_func(plot_figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col,
			  'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

	# Run 1st derivative function. The derivative can maybe help resolve peak borders.
	savgol_1st_deriv(xs, ys, sample_idx, param_dict, master_dict, user_input_dict, plot_dict)

	# Plotting the baseline values into interactive plot and the stativ "plotting plot"
	plot_func(figure, graph_name + '_baseline', baseline_values[0], baseline_values[1], 'None', 'lines',
			  x_title, y_title, subplot_row, subplot_col, 'AKTA_baseline', sample_idx, param_dict, 'N/A', 'N/A',
			  user_input_dict, master_dict)
	plot_func(plot_figure, graph_name + '_baseline', baseline_values[0], baseline_values[1], 'None', 'lines',
			  x_title, y_title, subplot_row, subplot_col, 'AKTA_baseline', sample_idx, param_dict, 'N/A', 'N/A',
			  user_input_dict, master_dict)

	# Check if the user has specified a peak. If this is the case minimum one ; must be placed
	if str(df['AKTA_fraction'][sample_idx]).count(';') > 0:

		peak_count = 0

		# Check if more than one peak is specified.
		# If only one peak is available i.e. no "|" is given the script just uses the cell value.
		if str(df['AKTA_fraction'][sample_idx]).count('|') > 0:
			peak_list = list(df['AKTA_fraction'][sample_idx].replace(',', '.').split('|'))
		else:
			peak_list = [df['AKTA_fraction'][sample_idx].replace(',', '.')]

		# Create a list for storing the AUC values for the peaks.
		peak_AUC_list = []

		for my_peak in peak_list:

			peak_count = peak_count + 1

			# Write values calculated above into master dict
			master_dict['ID_list_new'].append(graph_name)
			master_dict['notes_list'].append(user_input_dict['sample_notes'][sample_idx])
			master_dict['peak_id'].append('Peak' + str(peak_count))
			master_dict['sample_areas_tot'].append(round(AUC_sample_tot,3))
			#master_dict['sample_areas_tot'].append(AUC_sample_tot)
			master_dict['baseline_area_tot'].append(round(AUC_baseline_tot,3))

			# Get fraction boundaries from excel and turn them to floats
			fraction = my_peak.split(';')
			#fraction = list(df['AKTA_fraction'][sample_idx].replace(',', '.').split(';'))
			for a in range(0, len(fraction)):
				fraction[a] = float(fraction[a])
			frac_volume = abs(fraction[1] - fraction[0])

			# Slicing AKTA data to only get data points that are part of the fraction interval
			xsys = pd.concat([xs, ys], axis=1)
			xsys_slice = xsys[(xsys[x_id] >= fraction[0]) & (xsys[x_id] <= fraction[1])]

			# Slicing baseline data to only get data points that are part of fraction interval
			# Create a dataframe from the lists containing the baseline values
			xsys_baseline = pd.DataFrame({str(x_id): baseline_values[0], str(y_id): baseline_values[1]})
			xsys_baseline_slice = xsys_baseline[(xsys_baseline[x_id] >= fraction[0]) & (xsys_baseline[x_id] <= fraction[1])]

			# The retention time/volume is defined as the x value where the y value hits max
			Retention_frac = float(xsys_slice[x_id][xsys_slice[y_id] == max(xsys_slice[y_id])])

			# Calculate AUC for the specified fraction
			AUC_frac = auc(xsys_slice[x_id], xsys_slice[y_id])

			# Calculate AUC of the baseline in the fraction
			AUC_frac_baseline = auc(xsys_baseline_slice[x_id], xsys_baseline_slice[y_id])

			# Do calculation on the total fraction AUC by subtracting fraction AUC and fraction baseline AUC
			AUC_calculation = AUC_frac - AUC_frac_baseline
			peak_AUC_list.append(AUC_calculation)

			# Calculate fraction yield using the user specified extinction coefficient and the path length. Convert to mg to dividing by 1000
			frac_yield = AUC_calculation / (float(ext_coeff[sample_idx]) * param_dict['AKTA_pathlength'])
			frac_yield = frac_yield / 1000

			# Calculating the concentration in the fraction
			fraction_conc = frac_yield / frac_volume

			# Calculating the culture yield based on the fraction yield and the volume loaded.
			culture_yield = (frac_yield / volume_load[sample_idx] * 1000)

			# Putting all the calculated values into master dict in a re-formatted way
			Retention_frac = "{:.3f}".format(Retention_frac)
			master_dict['fraction_retentions'].append(Retention_frac)
			AUC_frac = "{:.3f}".format(AUC_frac)
			master_dict['fraction_areas'].append(AUC_frac)
			AUC_frac_baseline = "{:.3f}".format(AUC_frac_baseline)
			master_dict['fraction_baseline'].append(AUC_frac_baseline)
			AUC_calculation = "{:.3f}".format(AUC_calculation)
			master_dict['fraction_calculation'].append(AUC_calculation)
			frac_yield = "{:.3f}".format(frac_yield)
			master_dict['fraction_yield'].append(frac_yield)
			culture_yield = "{:.3f}".format(culture_yield)
			master_dict['culture_yield'].append(culture_yield)
			fraction_conc = "{:.3f}".format(fraction_conc)
			master_dict['fraction_concentrations'].append(fraction_conc)

			# Plotting the fraction
			plot_func(figure, graph_name + '_peak' + str(peak_count), xsys_slice[x_id], xsys_slice[y_id], 'N/A', 'lines', x_title,
					  y_title, subplot_row, subplot_col, 'AKTA_fraction', sample_idx, param_dict, color_list,
					  color_count, user_input_dict, master_dict)
			plot_func(plot_figure, graph_name + '_peak' + str(peak_count), xsys_slice[x_id], xsys_slice[y_id], 'N/A', 'lines', x_title,
					  y_title, subplot_row, subplot_col, 'AKTA_fraction', sample_idx, param_dict, color_list,
					  color_count, user_input_dict, master_dict)

			# Add the appropriate color to the table coloring list
			table_color_list_manager(color_list[color_count], plot_dict['table_color_list'])

		peak_frac_calculation(peak_AUC_list, AUC_sample_tot, AUC_baseline_tot, master_dict)

	# If no fraction is specified 'N/A' values are appended to the master dict.
	# This is to make sure value positions in the table are not shifted
	else:
		master_dict['fraction_retentions'].append('N/A')
		master_dict['fraction_areas'].append('N/A')
		master_dict['fraction_baseline'].append('N/A')
		master_dict['fraction_calculation'].append('N/A')
		master_dict['fraction_yield'].append('N/A')
		master_dict['culture_yield'].append('N/A')
		master_dict['fraction_concentrations'].append('N/A')


def akta_plotly_table(plot_dict, master_dict, user_input_dict):
	# Getting the figure from the plotting dictionary
	figure = plot_dict['figure']

	figure.add_trace(
		go.Table(
			header=dict(values=[
				'Sample ID',
				'Sample notes',
				'Total sample area',
				'Total baseline area',
				'Retention time/volume (beta)',
				'Area of fraction',
				'Baseline area (fraction)',
				'Area used for calculation',
				'Fraction yield [mg]',
				'Culture yield [ug/mL]'],
				align='left'),
			cells=dict(values=[
				user_input_dict['ID_list'],
				master_dict['notes_list'],
				master_dict['sample_areas_tot'],
				master_dict['baseline_area_tot'],
				master_dict['fraction_retentions'],
				master_dict['fraction_areas'],
				master_dict['fraction_baseline'],
				master_dict['fraction_calculation'],
				master_dict['fraction_yield'],
				master_dict['culture_yield']],
				align='left', height=50)))


def akta_baseline(baseline, x_val, y_val, mode, param_dict):
	baseline_values = [[], []]

	for k in range(0, len(baseline)):
		baseline[k] = float(baseline[k])

	if len(baseline) == 1:
		baseline_values = [x_val, [baseline[0]] * len(x_val)]

	elif len(baseline) == 2:

		# Calculate flat baseline similar to Unicorn software
		baseline_x1 = float(x_val.iloc[(x_val - baseline[0]).abs().argsort().iloc[:1]])
		baseline_y1 = float(y_val.iloc[(x_val - baseline[0]).abs().argsort().iloc[:1]])
		baseline_x2 = float(x_val.iloc[(x_val - baseline[1]).abs().argsort().iloc[:1]])
		baseline_y2 = float(y_val.iloc[(x_val - baseline[1]).abs().argsort().iloc[:1]])

		if mode == 'linear':
			a = (baseline_y2 - baseline_y1) / (baseline_x2 - baseline_x1)  # Calculate slope of linar baseline curve
			b = baseline_y1 - (a * baseline_x1)  # Calculate intercept of linear baseline curve

			baseline_values = [x_val, linear_model(x_val, a, b)]

		elif mode == 'peakutils':
			x_val_baseline = x_val[(x_val >= baseline_x1) & (x_val <= baseline_x2)]
			y_val_baseline = y_val[(x_val >= baseline_x1) & (
					x_val <= baseline_x2)]  # Getting all y values that are within the baseline interval

			baseline_values = [x_val_baseline, peakutils.baseline(y_val_baseline, deg=param_dict['AKTA baseline deg'])]

	return baseline_values

def peak_frac_calculation(peak_AUC_list, AUC_sample_tot, AUC_baseline_tot, master_dict):

	# The function calculates the relative AUC of the user specified peak compared to the total peak area.
	# It compares the individal peaks to the total peak area as well as the total sample area.

	# Convert list to array
	peak_AUC_list = np.array(peak_AUC_list)

	# Get the sum of the peak AUC values
	peak_sum = np.sum(peak_AUC_list)

	# Calculate the total AUC in the sample by subtracting total baseline AUC from total sample AUC
	AUC_tot = AUC_sample_tot - AUC_baseline_tot

	for i in peak_AUC_list:

		# Calculate the relationship between the peak AUC and the total AUC of peak.
		var1 = (i/peak_sum)*100

		# Calculate the relationship between the peak AUC and the total AUC of the sample.
		var2 = (i/AUC_tot)*100

		# Append the calculated values to the table.
		master_dict['fraction_AUC_of_total'].append('Of peak=' + str(round(var1, 1)) + '<br>' +
													'Of sample='+ str(round(var2, 1)))


def linear_model(x, a, b):
	y = x * a + b
	return y

