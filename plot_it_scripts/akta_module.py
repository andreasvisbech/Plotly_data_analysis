# Import modules
from sklearn.metrics import auc
import peakutils
from lmfit.models import GaussianModel, SkewedGaussianModel, ExponentialGaussianModel, LorentzianModel, \
							SplitLorentzianModel, VoigtModel, PseudoVoigtModel, PowerLawModel, Pearson7Model

# Import functions from other scripts
#from plot_it_scripts.plotting_script import *
from plot_it_scripts.main_functions import *
from plot_it_scripts.data_manipulation import *



def akta_data_slice(x_id, y_id, sample_idx, xs, ys, user_input_dict):

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
	baseline_bounds = user_input_dict['baseline_bounds']

	if baseline_bounds[sample_idx] != 'None':
		if str(baseline_bounds[sample_idx]).count(';') == 0:
			baseline_coord = str(baseline_bounds[sample_idx]).replace(',', '.')
			baseline_coord = [float(baseline_coord)]
		elif str(baseline_bounds[sample_idx]).count(';') == 1:
			baseline_coord = str(baseline_bounds[sample_idx]).replace(',', '.').split(';')
		else:
			raise ValueError('Something seems wrong with the specified baseline interval..?')

	else:
		baseline_coord = [0]

	# Calculate values for the baseline.
	baseline_xs, baseline_ys = akta_baseline(baseline_coord, xs, ys, param_dict['AKTA baseline mode'], param_dict)

	# Only plot the baseline if a baseline has been specified
	if str(baseline_bounds[sample_idx]) != 'None':
		plot_func(figure, graph_name + '_baseline', baseline_xs, baseline_ys, 'None', 'lines',
				  x_title, y_title, subplot_row, subplot_col, 'AKTA_baseline', sample_idx, param_dict,
				  color_list, color_count, user_input_dict, master_dict)

		plot_func(plot_figure, graph_name + '_baseline', baseline_xs, baseline_ys, 'None', 'lines',
				  x_title, y_title, subplot_row, subplot_col, 'AKTA_baseline', sample_idx, param_dict,
				  color_list, color_count, user_input_dict, master_dict)

	# Calculate the total AUC from entire sample
	AUC_sample_tot = auc(xs, ys)

	# Calculate AUC for the baseline
	AUC_baseline_tot = auc(baseline_xs, baseline_ys)

	# Plotting the data values
	plot_func(figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col,
			  'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
	plot_func(plot_figure, graph_name, xs, ys, 'None', plot_marker, x_title, y_title, subplot_row, subplot_col,
			  'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

	# Run 1st derivative function. The derivative can maybe help resolve peak borders.
	savgol_1st_deriv(xs, ys, sample_idx, param_dict, master_dict, user_input_dict, plot_dict)

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
			master_dict['baseline_area_tot'].append(round(AUC_baseline_tot,3))

			# Get fraction boundaries from excel and turn them to floats
			fraction = my_peak.split(';')
			for a in range(0, len(fraction)):
				fraction[a] = float(fraction[a])
			frac_volume = abs(fraction[1] - fraction[0])

			# Slicing AKTA data to only get data points that are part of the fraction interval
			xsys = pd.DataFrame({x_id: xs, y_id: ys})
			xsys_slice = xsys[(xsys[x_id] >= fraction[0]) & (xsys[x_id] <= fraction[1])]

			# Slicing baseline data to only get data points that are part of fraction interval
			# Create a dataframe from the lists containing the baseline values
			xsys_baseline = pd.DataFrame({str(x_id): baseline_xs, str(y_id): baseline_ys})
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

			# Re-plot the baseline. This is important for the fill coloring of the fractions.
			# The filling function uses the previous trace so the baseline must be plotted again since the
			# peaks are plotted right after.
			plot_func(figure, graph_name + '_baseline', xsys_baseline_slice[x_id], xsys_baseline_slice[y_id], 'None', 'lines',
					  x_title, y_title, subplot_row, subplot_col, 'AKTA_baseline', sample_idx, param_dict,
					  color_list, color_count, user_input_dict, master_dict)
			plot_func(plot_figure, graph_name + '_baseline', xsys_baseline_slice[x_id], xsys_baseline_slice[y_id], 'None',
					  'lines',
					  x_title, y_title, subplot_row, subplot_col, 'AKTA_baseline', sample_idx, param_dict,
					  color_list, color_count, user_input_dict, master_dict)

			# Plotting the fraction
			plot_func(figure, graph_name + '_peak' + str(peak_count), xsys_slice[x_id], xsys_slice[y_id], 'N/A', 'lines', x_title,
					  y_title, subplot_row, subplot_col, 'AKTA_fraction', sample_idx, param_dict, color_list,
					  color_count, user_input_dict, master_dict)
			plot_func(plot_figure, graph_name + '_peak' + str(peak_count), xsys_slice[x_id], xsys_slice[y_id], 'N/A', 'lines', x_title,
					  y_title, subplot_row, subplot_col, 'AKTA_fraction', sample_idx, param_dict, color_list,
					  color_count, user_input_dict, master_dict)

			# Add the appropriate color to the table coloring list
			table_color_list_manager(user_input_dict['plot_color'], plot_dict['table_color_list'])
			#table_color_list_manager(color_list[color_count], plot_dict['table_color_list'])

		# If fit models have been supplied by user the script does peak fitting.
		# First calculate the fit through function and then plot the individual fits for thep peaks
		if user_input_dict['fit_models'][sample_idx] != 'None':
			fit_x, fit_y = peak_fit(sample_idx, xsys, x_id, y_id, peak_list, user_input_dict)
			for k, pk in enumerate(peak_list):
				plot_func(figure, graph_name + '_peak_fit' + str(k), fit_x[k], fit_y[k], 'None',
						  'lines', x_title, y_title, subplot_row, subplot_col, 'AKTA_peak_fit',
						  sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

		peak_frac_calculation(peak_AUC_list, AUC_sample_tot, AUC_baseline_tot, master_dict)

def akta_baseline(baseline, x_val, y_val, mode, param_dict):

	#baseline_values = [[], []]

	# Make sure the data is in arrays
	x_val = np.array(x_val).astype(float)
	y_val = np.array(y_val).astype(float)
	baseline = np.array(baseline).astype(float)

	# If only a single value is given, the baseline data is simply this value and then using the full x range
	if len(baseline) == 1:
		baseline_values = [x_val, [baseline[0]] * len(x_val)]
		baseline_xs = x_val
		baseline_ys = [baseline[0]] * len(x_val)

	elif len(baseline) == 2:

		base_low = baseline[0]
		base_high = baseline[1]

		# We then need to find the values in the data that is actually closest to the baseline bounds.
		# This is done by taking difference between baseline bounds and the x_vals and finding difference minimum
		diff_x1 = np.absolute(x_val - base_low)
		idx_x1 = diff_x1.argmin()
		baseline_x1 = x_val[idx_x1]
		baseline_y1 = y_val[idx_x1]

		diff_x2 = np.absolute(x_val - base_high)
		idx_x2 = diff_x2.argmin()
		baseline_x2 = x_val[idx_x2]
		baseline_y2 = y_val[idx_x2]

		# Calculate flat baseline similar to Unicorn software
		#baseline_x1 = float(x_val.iloc[(x_val - baseline[0]).abs().argsort().iloc[:1]])
		#baseline_y1 = float(y_val.iloc[(x_val - baseline[0]).abs().argsort().iloc[:1]])
		#baseline_x2 = float(x_val.iloc[(x_val - baseline[1]).abs().argsort().iloc[:1]])
		#baseline_y2 = float(y_val.iloc[(x_val - baseline[1]).abs().argsort().iloc[:1]])

		if mode == 'linear':
			a = (baseline_y2 - baseline_y1) / (baseline_x2 - baseline_x1)  # Calculate slope of linar baseline curve
			b = baseline_y1 - (a * baseline_x1)  # Calculate intercept of linear baseline curve

			baseline_xs = x_val
			baseline_ys = linear_model(x_val, a, b)

			baseline_values = [x_val, linear_model(x_val, a, b)]

		elif mode == 'peakutils':
			x_val_baseline = x_val[(x_val >= baseline_x1) & (x_val <= baseline_x2)]
			y_val_baseline = y_val[(x_val >= baseline_x1) & (
					x_val <= baseline_x2)]  # Getting all y values that are within the baseline interval

			baseline_values = [x_val_baseline, peakutils.baseline(y_val_baseline, deg=param_dict['AKTA baseline deg'])]

			baseline_xs = x_val_baseline
			baseline_ys = peakutils.baseline(y_val_baseline, deg=param_dict['AKTA baseline deg'])

		else:
			baseline_xs = 0
			baseline_ys = 0

	else:
		baseline_xs = 0
		baseline_ys = 0

	return baseline_xs, baseline_ys

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

def trace_stacking(i, ys, param_dict, master_dict):

	# The function moves traces vertially to allow stacking of data traces. It uses percent so if 100% is given
	# the function will place the lowest point of the top graph at the same y value as the highest point on the
	# graph below.

	# Setting the new ys value to the old ys value just as a default
	ys_new = ys

	# The default value is zero. The code should only change the traces if the default has been changed
	if param_dict['trace_stacking'] != 0:

		# Getting the factor that will be used to move traces
		factor = param_dict['trace_stacking'] / 100

		# The trace should only be stacked if i is greater than zero i.e. the first trace is not stacked.
		if i > 0:

			ys_new = ys + master_dict['trace_stack_jump']

		# Get the max value of y values and put it to master dict. Used to know how much next trace needs to be moved.
		master_dict['trace_stack_jump'] = master_dict['trace_stack_jump'] + max(ys)*factor

	return ys_new

def peak_fit(idx, xsys, x_id, y_id, peak_list, user_input_dict):

	# Get xs and ys as arrays and make sure they are floats
	xs = xsys[x_id].to_numpy().astype(float)
	ys = xsys[y_id].to_numpy().astype(float)

	# Defining lists for storing user specified peak limits
	peak_lim_low = []
	peak_lim_high = []

	peak_count = len(peak_list)

	# Getting the peak limits and putting them into the lists.
	for pk in peak_list:
		low = pk.split(';')[0]
		peak_lim_low.append(float(low))
		high = pk.split(';')[1]
		peak_lim_high.append(float(high))

	# Getting the initial guesses from the user.
	# If no guesses supplied will use center of the peak interval.
	fit_flags = user_input_dict['fit_misc'][idx].split('|')
	for j in fit_flags:
		if 'center_guess' in j:
			center_vals = j.split('=')[1]
			center_guesses = center_vals.split(';')
			center_guesses = np.array(center_guesses).astype(float)
		else:
			center_guesses = []
			for k, pk in enumerate(peak_list):
				cen_var = (peak_lim_high[k] - peak_lim_low[k])/2
				center_guesses.append(cen_var)

	# Create a dict for storing the models
	model_dict = {}

	# Getting the user specified models. Default is Gaussian
	models = user_input_dict['fit_models'][idx].lower().split(';')

	# Set the model for the first peak.
	# Note that the model is called "zero" simply to make more compatible with loop below
	if models[0] == 'gaussian':
		model_dict['mod0'] = GaussianModel(prefix='mod0_')
		pars = model_dict['mod0'].guess(ys, x=xs)
		#pars['mod0_center'].set(value=center_guesses[0], min=peak_lim_low[0], max=peak_lim_high[0])
		pars['mod0_center'].set(value=center_guesses[0], min=center_guesses[0]-0.1, max=center_guesses[0]+0.1)
		pars['mod0_sigma'].set(value=0.2, max=(peak_lim_high[0]-peak_lim_low[0]))
		pars['mod0_amplitude'].set(value=10, min=0.1)

	elif models[0] in ['skewedgaussian', 'skewed_gaussian']:
		model_dict['mod0'] = SkewedGaussianModel(prefix='mod0_')
		pars = model_dict['mod0'].guess(ys, x=xs)
		pars['mod0_center'].set(value=center_guesses[0], min=peak_lim_low[0], max=peak_lim_high[0])
		pars['mod0_sigma'].set(value=0.2, max=(peak_lim_high[0] - peak_lim_low[0]))
		pars['mod0_amplitude'].set(value=10, min=0.1)
		pars['mod0_gamma'].set(value=0.0)

	elif models[0] in ['expgaussian', 'exp_gaussian']:
		model_dict['mod0'] = ExponentialGaussianModel(prefix='mod0_')
		pars = model_dict['mod0'].guess(ys, x=xs)
		pars['mod0_center'].set(value=center_guesses[0], min=peak_lim_low[0], max=peak_lim_high[0])
		pars['mod0_sigma'].set(value=0.2, min=0.1, max=(peak_lim_high[0] - peak_lim_low[0]))
		pars['mod0_amplitude'].set(value=10, min=0.1, max=100)
		pars['mod0_gamma'].set(value=0.0, min=0.1, max=50)

	elif models[0] in ['lorentzian']:
		model_dict['mod0'] = LorentzianModel(prefix='mod0_')
		pars = model_dict['mod0'].guess(ys, x=xs)
		pars['mod0_center'].set(value=center_guesses[0], min=peak_lim_low[0], max=peak_lim_high[0])
		pars['mod0_sigma'].set(value=0.2, max=(peak_lim_high[0] - peak_lim_low[0]))
		pars['mod0_amplitude'].set(value=10, min=0.1)

	# Set the model for the other peaks
	for i in range(1, len(models)):
		if models[i] == 'gaussian':
			model_dict['mod' + str(i)] = GaussianModel(prefix=str('mod' + str(i) + '_'))
			pars.update(model_dict['mod' + str(i)].make_params())
			#pars['mod' + str(i) + '_center'].set(value=center_guesses[i], min=peak_lim_low[i], max=peak_lim_high[i])
			pars['mod' + str(i) + '_center'].set(value=center_guesses[i], min=center_guesses[i]-0.1, max=center_guesses[i]+0.1)
			pars['mod' + str(i) + '_sigma'].set(value=0.1, max=(peak_lim_high[i] - peak_lim_low[i])/4)
			pars['mod' + str(i) + '_amplitude'].set(value=10, min=0.1)

		elif models[i] in ['skewedgaussian', 'skewed_gaussian']:
			model_dict['mod' + str(i)] = SkewedGaussianModel(prefix=str('mod' + str(i) + '_'))
			pars.update(model_dict['mod' + str(i)].make_params())
			pars['mod' + str(i) + '_center'].set(value=center_guesses[i], min=center_guesses[i]-0.1, max=center_guesses[i]+0.1)
			pars['mod' + str(i) + '_sigma'].set(value=0.1, max=(peak_lim_high[i] - peak_lim_low[i])/4)
			pars['mod' + str(i) + '_gamma'].set(value=0.2, min=0.0, max=10)

		elif models[i] in ['expgaussian', 'exp_gaussian']:
			model_dict['mod' + str(i)] = ExponentialGaussianModel(prefix=str('mod' + str(i) + '_'))
			pars.update(model_dict['mod' + str(i)].make_params())
			pars['mod' + str(i) + '_center'].set(value=center_guesses[i], min=peak_lim_low[i], max=peak_lim_high[i])
			pars['mod' + str(i) + '_sigma'].set(value=1, min=0.1, max=(peak_lim_high[0] - peak_lim_low[0])/8)
			pars['mod' + str(i) + '_amplitude'].set(value=10, min=0.1, max=50)
			pars['mod' + str(i) + '_gamma'].set(value=1, min=0.1, max=30)

		elif models[i] in ['lorentzian']:
			model_dict['mod' + str(i)] = LorentzianModel(prefix=str('mod' + str(i) + '_'))
			pars.update(model_dict['mod' + str(i)].make_params())
			pars['mod' + str(i) + '_center'].set(value=center_guesses[i], min=peak_lim_low[i], max=peak_lim_high[i])
			pars['mod' + str(i) + '_sigma'].set(value=0.1, max=(peak_lim_high[i] - peak_lim_low[i]))

	# Define a final model based on the number of peaks
	if peak_count == 1:
		final_mod = model_dict['mod0']
	elif peak_count == 2:
		final_mod = model_dict['mod0'] + model_dict['mod1']
	elif peak_count == 3:
		final_mod = model_dict['mod0'] + model_dict['mod1'] + model_dict['mod2']
	elif peak_count == 4:
		final_mod = model_dict['mod0'] + model_dict['mod1'] + model_dict['mod2'] + model_dict['mod3']
	elif peak_count == 5:
		final_mod = model_dict['mod0'] + model_dict['mod1'] + model_dict['mod2'] + model_dict['mod3'] + model_dict['mod4']

	init = final_mod.eval(pars, x=xs)
	out = final_mod.fit(ys, pars, x=xs)

	print(out.fit_report())

	comps = out.eval_components(x=xs)

	xs_out = []
	ys_out = []
	for i in range(peak_count):
		xs_out.append(xs)
		ys_out.append(comps['mod' + str(i) + '_'])

	return xs_out, ys_out





