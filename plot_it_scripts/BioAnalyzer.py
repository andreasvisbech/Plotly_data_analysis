# import relevant modules
from scipy.signal import find_peaks, peak_widths
from lmfit import Model

# Import functions from other scripts
from plot_it_scripts.akta_module import *

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

	# Check if user has prompted for peak analysis
	if param_dict['peak_analysis'] == 'yes':

		# Calculate total AUC for sample
		AUC_total = auc(xs, ys)

		# Define a list for AUC values
		AUC_peak_list = []

		# If no peak limits have been specified the peak detection is automatic
		if user_input_dict['peak_limits'][sample_idx] == 'None':
			peak_list = automatic_peak_detection(xs, ys)

		# Else use the inputs from user.
		else:
			peak_list = peak_analysis_with_fixed_limits(sample_idx, xs, ys, user_input_dict)

		# Go over the individual peaks and perform analysis.
		for j, pk in enumerate(peak_list):

			# Make entries in master dict
			master_dict['ID_list_new'].append(graph_name)
			master_dict['notes_list'].append(user_input_dict['sample_notes'][sample_idx])
			master_dict['peak_id'].append('Peak' + str(j+1))

			peak = pk.split(';')

			xs_peak = xs[(xs >= float(peak[0])) & (xs <= float(peak[1]))]
			ys_peak = ys[(xs >= float(peak[0])) & (xs <= float(peak[1]))]

			# Calculate AUC values for the peak
			AUC_peak = auc(xs_peak, ys_peak)
			AUC_peak_list.append(AUC_peak)

			# Plotting the fraction
			plot_func(figure, graph_name + '_peak' + str(j+1), xs_peak, ys_peak, 'N/A', 'lines',
					  x_title, y_title, subplot_row, subplot_col, 'bioanalyzer_fraction', sample_idx, param_dict,
					  color_list, color_count, user_input_dict, master_dict)
			plot_func(plot_figure, graph_name + '_peak' + str(j + 1), xs_peak, ys_peak, 'N/A', 'lines',
					  x_title, y_title, subplot_row, subplot_col, 'bioanalyzer_fraction', sample_idx, param_dict,
					  color_list, color_count, user_input_dict, master_dict)

			# Add the appropriate color to the table coloring list
			table_color_list_manager(user_input_dict['plot_color'], plot_dict['table_color_list'])

		AUC_sum = np.array(AUC_peak_list).sum()

		for auc_val in AUC_peak_list:

			# Calculate the AUC relative to the sum of AUC peak values.
			AUC_calc1 = auc_val/AUC_sum*100
			master_dict['fraction_calculation'].append(round(AUC_calc1, 2))

			AUC_calc2 = auc_val/AUC_total*100
			master_dict['fraction_AUC_of_total'].append(round(AUC_calc2, 2))

	# Check if user has prompted for peak analysis
	if param_dict['peak_fitting'] == 'yes':

		#peaks, peak_widths = bioanalyzer_detect_peak(xs, ys)

		#peak_list = generate_peak_list(xs, ys, peaks, peak_widths)

		#user_input_dict['fit_models'][sample_idx] = 'skewed_gaussian'
		user_input_dict['fit_misc'][sample_idx] = user_input_dict['fit_misc'][sample_idx] + '|center_guesses=37.2'

		#plot_func(figure, 'test_'+str(sample_idx), xs[peaks], ys[peaks], 'None', 'dots', x_title, y_title, subplot_row, subplot_col,
		#		  'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

		#AUC_tot = auc(xs, ys)

		xsys = pd.DataFrame({x_id: xs, y_id:ys})
		fit_x, fit_y = bioanalyzer_peak_fit(sample_idx, xsys, x_id, y_id, peak_list, user_input_dict, param_dict)
		for k, pk in enumerate(peak_list):
			plot_func(figure, graph_name + '_peak_fit' + str(k+1), fit_x[k], fit_y[k], 'None',
					  'lines', x_title, y_title, subplot_row, subplot_col, 'AKTA_peak_fit',
					  sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

		plot_func(figure, graph_name + '_peak_fit_tot', xs, np.array(fit_y).sum(axis=0), 'None',
				  'lines', x_title, y_title, subplot_row, subplot_col, 'AKTA_peak_fit',
				  sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)
		#print(np.array(fit_y).sum(axis=0))

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

def automatic_peak_detection(xs, ys):

	from scipy.signal import find_peaks, peak_widths

	# Make list for storing peak values
	peak_list = []

	# Set prominence as 2% of max signal value
	ymax = max(ys)
	prom = (ymax/100) * 2

	peaks, _ = find_peaks(ys, height=2, prominence=prom)

	# Define peak widths using scipy module
	peak_width = peak_widths(ys, peaks, rel_height=0.07)

	for j, pk in enumerate(peaks):

		#print(pk)

		peak_low = xs[pk] - (peak_width[0][j]/2)
		peak_high = xs[pk] + (peak_width[0][j]/2)
		peak_str = str(peak_low) + ';' + str(peak_high)
		peak_list.append(peak_str)

	return peak_list

def generate_peak_list(xs, ys, peaks, peak_widths):

	peak_list = []

	peaks_low = xs[peaks] - (peak_widths/32)
	peaks_high = xs[peaks] + (peak_widths/32)

	for i, pk in enumerate(peaks):
		my_peak = str(peaks_low[i]) + ';' + str(peaks_high[i])
		peak_list.append(my_peak)

	return peak_list

def peak_analysis_with_fixed_limits(sample_idx, xs, ys, user_input_dict):

	peak_limits = str(user_input_dict['peak_limits'][sample_idx])

	if peak_limits.count('|') > 0:
		peak_list = list(peak_limits.replace(',', '.').split('|'))
	else:
		peak_list = [peak_limits.replace(',', '.')]

	return peak_list

def bioanalyzer_peak_fit(idx, xsys, x_id, y_id, peak_list, user_input_dict, param_dict):

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

	fit_param_dict = get_fit_params_init(idx, peak_list, user_input_dict)

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
		pars['mod0_sigma'].set(value=0.2, min=fit_param_dict['sigma_min'][0], max=fit_param_dict['sigma_max'][0])
		pars['mod0_amplitude'].set(value=10, min=0.1)
		pars['mod0_gamma'].set(value=0.0, min=fit_param_dict['gamma_min'][0], max=fit_param_dict['gamma_max'][0])

	elif models[0] in ['expgaussian', 'exp_gaussian']:
		model_dict['mod0'] = ExponentialGaussianModel(prefix='mod0_')
		pars = model_dict['mod0'].guess(ys, x=xs)
		pars['mod0_center'].set(value=center_guesses[0], min=peak_lim_low[0], max=peak_lim_high[0])
		pars['mod0_sigma'].set(value=0.2, min=fit_param_dict['sigma_min'][0], max=fit_param_dict['sigma_max'][0])
		pars['mod0_amplitude'].set(value=10, min=0.1, max=100)
		pars['mod0_gamma'].set(value=0.0, min=fit_param_dict['gamma_min'][0], max=fit_param_dict['gamma_max'][0])

	elif models[0] in ['lorentzian']:
		model_dict['mod0'] = LorentzianModel(prefix='mod0_')
		pars = model_dict['mod0'].guess(ys, x=xs)
		pars['mod0_center'].set(value=center_guesses[0], min=peak_lim_low[0], max=peak_lim_high[0])
		pars['mod0_sigma'].set(value=0.1, max=(peak_lim_high[0] - peak_lim_low[0]))
		pars['mod0_amplitude'].set(value=10, min=0.1)

	elif models[0] in ['pearson7model', 'pearson7']:
		model_dict['mod0'] = Pearson7Model(prefix='mod0_')
		pars = model_dict['mod0'].guess(ys, x=xs)
		pars['mod0_center'].set(value=center_guesses[0], min=peak_lim_low[0], max=peak_lim_high[0])
		pars['mod0_sigma'].set(value=0.2, max=10)
		pars['mod0_amplitude'].set(value=10, min=0.1)
		pars['mod0_expon'].set(value=2, min=2)



	# Set the model for the other peaks
	for i in range(1, len(models)):
		if models[i] == 'gaussian':
			model_dict['mod' + str(i)] = GaussianModel(prefix=str('mod' + str(i) + '_'))
			pars.update(model_dict['mod' + str(i)].make_params())
			#pars['mod' + str(i) + '_center'].set(value=center_guesses[i], min=peak_lim_low[i], max=peak_lim_high[i])
			pars['mod' + str(i) + '_center'].set(value=center_guesses[i], min=center_guesses[i]-0.1, max=center_guesses[i]+0.1)
			pars['mod' + str(i) + '_sigma'].set(value=0.1, min=fit_param_dict['sigma_min'][i], max=fit_param_dict['sigma_max'][i])
			pars['mod' + str(i) + '_amplitude'].set(value=10, min=0.1)

		elif models[i] in ['skewedgaussian', 'skewed_gaussian']:
			model_dict['mod' + str(i)] = SkewedGaussianModel(prefix=str('mod' + str(i) + '_'))
			pars.update(model_dict['mod' + str(i)].make_params())
			pars['mod' + str(i) + '_center'].set(value=center_guesses[i], min=center_guesses[i]-0.1, max=center_guesses[i]+0.1)
			pars['mod' + str(i) + '_sigma'].set(value=0.1, min=fit_param_dict['sigma_min'][i], max=fit_param_dict['sigma_max'][i])
			pars['mod' + str(i) + '_gamma'].set(value=0.2, min=fit_param_dict['gamma_min'][i], max=fit_param_dict['gamma_max'][i])

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

	if param_dict['fit_report'] == 'yes':
		print(out.fit_report())

	comps = out.eval_components(x=xs)

	xs_out = []
	ys_out = []
	for i in range(peak_count):
		xs_out.append(xs)
		ys_out.append(comps['mod' + str(i) + '_'])

	return xs_out, ys_out

def get_fit_params_init(idx, peak_list, user_input_dict):

	# Create empty dict for storing the values
	fit_param_dict = dict()

	fit_flags = user_input_dict['fit_misc'][idx].split('|')


	#for j in fit_flags:

		#if 'center_guess' in j:
		#	center_vals = j.split('=')[1]
		#	center_guesses = center_vals.split(';')
		#	center_guesses = np.array(center_guesses).astype(float)
		#else:
		#	center_guesses = []
		#	for k, pk in enumerate(peak_list):
		#		cen_var = (peak_lim_high[k] - peak_lim_low[k])/2
		#		center_guesses.append(cen_var)

	for j in fit_flags:
		if 'sigma_min' in j:
			sigma = j.split('=')[1]
			sigma_min = sigma.split(';')
			sigma_min = np.array(sigma_min).astype(float)
			fit_param_dict['sigma_min'] = sigma_min
		elif 'sigma_max' in j:
			sigma = j.split('=')[1]
			sigma_max = sigma.split(';')
			sigma_max = np.array(sigma_max).astype(float)
			fit_param_dict['sigma_max'] = sigma_max
		elif 'gamma_min' in j:
			gamma = j.split('=')[1]
			gamma_min = gamma.split(';')
			gamma_min = np.array(gamma_min).astype(float)
			fit_param_dict['gamma_min'] = gamma_min
		elif 'gamma_max' in j:
			gamma = j.split('=')[1]
			gamma_max = gamma.split(';')
			gamma_max = np.array(gamma_max).astype(float)
			fit_param_dict['gamma_max'] = gamma_max


	if not 'sigma_min' in fit_param_dict:
		sigma_min = []
		for k, pk in enumerate(peak_list):
			sigma_min.append(0.0)
	if not 'sigma_max' in fit_param_dict:
		sigma_max = []
		for k, pk in enumerate(peak_list):
			sigma_max.append(10)
	if not 'gamma_min' in fit_param_dict:
		gamma_min = []
		for k, pk in enumerate(peak_list):
			gamma_min.append(0.0)
	if not 'gamma_max' in fit_param_dict:
		gamma_max = []
		for k, pk in enumerate(peak_list):
			gamma_max.append(10)

	return fit_param_dict

def power_Gaussian(x, center, sigma, amplitude):
	function = amplitude * np.exp(-(x-center)**2/sigma)
	return function






