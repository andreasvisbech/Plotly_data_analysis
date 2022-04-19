"""
Author: Andreas
"""
# Load modules
import os
import argparse
import statistics
import pathlib
import tkinter as tk
from tkinter import ttk
import peakutils
import pandas as pd
import numpy as np

import plotly.graph_objects as go
from plotly.subplots import make_subplots

from sklearn.metrics import r2_score
from sklearn.metrics import auc

from scipy.stats import binned_statistic
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

__version__ = "3.0.0"

# globals
color_count_global = 0
python_misc_global = []
used_graph_names = []  # Create an empty list for storing used graph names
# Defining master dictionary
master_dict = {
	'notes_list': [], 'vertex_max': [], 'vertex_min': [], 'model_list': [], 'peak_onset': [], 'inflection_points': [],
	'fit_parameters': [], 'R_square': [], 'KD_fit': [], 'ID_list_new': [], 'table_coords': []
}
std_dev = []
error_values = []
subplot_col_count = 0


def parse_args():
	# initialize arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--plot_type", help="one of akta, fida, panta, bar, scatter or boxplot", required=True)
	parser.add_argument("-i", "--input_file", help="Excel file", type=pathlib.Path, required=True)
	parser.add_argument('-v', '--version', action='version', version=__version__)
	parser.add_argument("-fit", "--data_fit", help="y for fitting", action='store_true')
	parser.add_argument("-log", "--log_scale", help="make log scale", action='store_true')
	parser.add_argument("-out", "--output", type=pathlib.Path, default=pathlib.Path("./output"),
						help='the output path (default: %(default)s)')
	parser.add_argument("-plot", "--plotting", help="y/n to output", action='store_true')
	parser.add_argument("-advanced", "--advanced_option_box", help="y/n to advanced option box", action='store_true')
	return parser.parse_args()


def advanced_box():
	param_dict = {}

	def show_entry_fields():
		print('Values have been registered')
		param_dict['RI min'] = float(e1.get())
		param_dict['RI max'] = float(e2.get())
		param_dict['RIA min'] = float(e3.get())
		param_dict['RIA max'] = float(e4.get())
		param_dict['KD min'] = float(e5.get())
		param_dict['KD max'] = float(e6.get())
		param_dict['CI min'] = float(e7.get())
		param_dict['CI max'] = float(e8.get())
		param_dict['AKTA baseline mode'] = str(var.get())
		param_dict['AKTA baseline deg'] = int(var1.get())
		param_dict['AKTA_pathlength'] = float(e11.get())
		param_dict['AKTA neg val handle'] = str(var2.get())
		param_dict['graph width'] = float(e12.get())
		param_dict['xaxis_title_font_size'] = float(e14.get())
		param_dict['yaxis_title_font_size'] = float(e15.get())
		param_dict['xaxis_ticks_font_size'] = float(e16.get())
		param_dict['yaxis_ticks_font_size'] = float(e17.get())
		param_dict['plot_template'] = str(var3.get())
		param_dict['Baseline graph width'] = float(e19.get())
		param_dict['marker_size'] = float(e20.get())
		param_dict['vertex_point_window'] = int(e21.get())
		param_dict['peak_onset_var_deg'] = float(e22.get())
		param_dict['panta_intermediate_plot'] = str(var4.get())
		param_dict['Bmin_min'] = float(e24.get())
		param_dict['Bmin_max'] = float(e25.get())
		param_dict['Bmax_min'] = float(e26.get())
		param_dict['Bmax_max'] = float(e27.get())
		param_dict['KD_fit_min'] = float(e28.get())
		param_dict['KD_fit_max'] = float(e29.get())
		param_dict['k_coop_min'] = float(e30.get())
		param_dict['k_coop_max'] = float(e31.get())

	master = tk.Tk()
	ttk.Label(master, text="RI min [FIDA]").grid(row=0)
	ttk.Label(master, text='\t').grid(row=0, column=2)
	ttk.Label(master, text="RI max [FIDA]").grid(row=0, column=3)
	ttk.Label(master, text='\t').grid(row=0, column=5)
	ttk.Label(master, text="Graph width [Plot layot]").grid(row=0, column=6)
	ttk.Label(master, text="RIA min [FIDA]").grid(row=1)
	ttk.Label(master, text="RIA max [FIDA]").grid(row=1, column=3)
	ttk.Label(master, text="Baseline graph width [Plot layot]").grid(row=1, column=6)
	ttk.Label(master, text="Kd min [FIDA]").grid(row=2)
	ttk.Label(master, text="Kd max [FIDA]").grid(row=2, column=3)
	ttk.Label(master, text='Marker size').grid(row=2, column=6)
	ttk.Label(master, text="CI min [FIDA]").grid(row=3)
	ttk.Label(master, text="CI max [FIDA]").grid(row=3, column=3)
	ttk.Label(master, text='X axis title font size').grid(row=3, column=6)
	ttk.Label(master, text='Y axis title font size').grid(row=3, column=9)
	ttk.Label(master, text="Baseline mode [AKTA]").grid(row=4, column=0)
	ttk.Label(master, text="Plotting theme [layout]").grid(row=5, column=6)
	ttk.Label(master, text='X axis tick font size').grid(row=4, column=6)
	ttk.Label(master, text='Y axis tick font size').grid(row=4, column=9)
	ttk.Label(master, text="Baseline deg [AKTA]").grid(row=5, column=0)
	ttk.Label(master, text="Replace negative values with zero? [AKTA]").grid(row=6, column=0)
	ttk.Label(master, text="Path length in cm [AKTA]").grid(row=7, column=0)
	ttk.Label(master, text="Vertex point window [Panta]").grid(row=8, column=0)
	ttk.Label(master, text="Degree variation for onset (%) [Panta]").grid(row=9, column=0)
	ttk.Label(master, text="Plot intermediate data? [Panta]").grid(row=10, column=0)
	ttk.Label(master, text="Bmin min [scatter fit (Hill)]").grid(row=11, column=0)
	ttk.Label(master, text="Bmin max [scatter fit (Hill)]").grid(row=11, column=3)
	ttk.Label(master, text="Bmax min [scatter fit (Hill)]").grid(row=12, column=0)
	ttk.Label(master, text="Bmax max [scatter fit (Hill)]").grid(row=12, column=3)
	ttk.Label(master, text="KD min [scatter fit (Hill)]").grid(row=13, column=0)
	ttk.Label(master, text="KD max [scatter fit (Hill)]").grid(row=13, column=3)
	ttk.Label(master, text="k_coop min [scatter fit (Hill)]").grid(row=14, column=0)
	ttk.Label(master, text="k_coop max [scatter fit (Hill)]").grid(row=14, column=3)

	e1 = ttk.Entry(master)
	e2 = ttk.Entry(master)
	e3 = ttk.Entry(master)
	e4 = ttk.Entry(master)
	e5 = ttk.Entry(master)
	e6 = ttk.Entry(master)
	e7 = ttk.Entry(master)
	e8 = ttk.Entry(master)

	var = tk.StringVar(master)
	var.set('linear')
	e9 = ttk.OptionMenu(master, var, 'linear', 'peakutils')

	var1 = tk.StringVar(master)
	var1.set("1")
	e10 = ttk.OptionMenu(master, var1, "0", "1", "2", "3")

	var2 = tk.StringVar(master)
	var2.set('no')
	e13 = ttk.OptionMenu(master, var2, 'no', 'yes')

	e11 = ttk.Entry(master)
	e12 = ttk.Entry(master)
	e14 = ttk.Entry(master)
	e15 = ttk.Entry(master)
	e16 = ttk.Entry(master)
	e17 = ttk.Entry(master)

	var3 = tk.StringVar(master)
	var3.set('plotly')
	e18 = ttk.OptionMenu(master, var3, 'plotly', 'plotly_white', 'simple_white')

	e19 = ttk.Entry(master)
	e20 = ttk.Entry(master)
	e21 = ttk.Entry(master)
	e22 = ttk.Entry(master)

	var4 = tk.StringVar(master)
	var4.set('no')
	e23 = ttk.OptionMenu(master, var4, 'no', 'yes')

	e24 = ttk.Entry(master)
	e25 = ttk.Entry(master)
	e26 = ttk.Entry(master)
	e27 = ttk.Entry(master)
	e28 = ttk.Entry(master)
	e29 = ttk.Entry(master)
	e30 = ttk.Entry(master)
	e31 = ttk.Entry(master)

	e1.insert(10, "0")
	e2.insert(10, np.inf)
	e3.insert(10, "0")
	e4.insert(10, np.inf)
	e5.insert(10, "0")
	e6.insert(10, np.inf)
	e7.insert(10, "0")
	e8.insert(10, np.inf)
	e11.insert(10, "0.2")
	e12.insert(10, "2")
	e14.insert(10, "12")
	e15.insert(10, "12")
	e16.insert(10, "12")
	e17.insert(10, "12")
	e19.insert(10, "2")
	e20.insert(10, "10")
	e21.insert(10, "100")
	e22.insert(10, "0.5")
	e24.insert(10, "0")
	e25.insert(10, np.inf)
	e26.insert(10, "0")
	e27.insert(10, np.inf)
	e28.insert(10, "0")
	e29.insert(10, np.inf)
	e30.insert(10, -np.inf)
	e31.insert(10, np.inf)

	e1.grid(row=0, column=1)
	e2.grid(row=0, column=4)
	e3.grid(row=1, column=1)
	e4.grid(row=1, column=4)
	e5.grid(row=2, column=1)
	e6.grid(row=2, column=4)
	e7.grid(row=3, column=1)
	e8.grid(row=3, column=4)
	e9.grid(row=4, column=1)
	e10.grid(row=5, column=1)
	e13.grid(row=6, column=1)
	e11.grid(row=7, column=1)
	e12.grid(row=0, column=7)
	e14.grid(row=3, column=7)
	e15.grid(row=3, column=10)
	e16.grid(row=4, column=7)
	e17.grid(row=4, column=10)
	e18.grid(row=5, column=7)
	e19.grid(row=1, column=7)
	e20.grid(row=2, column=7)
	e21.grid(row=8, column=1)
	e22.grid(row=9, column=1)
	e23.grid(row=10, column=1)
	e24.grid(row=11, column=1)
	e25.grid(row=11, column=4)
	e26.grid(row=12, column=1)
	e27.grid(row=12, column=4)
	e28.grid(row=13, column=1)
	e29.grid(row=13, column=4)
	e30.grid(row=14, column=1)
	e31.grid(row=14, column=4)

	ttk.Button(master, text='Run', command=master.quit).grid(row=20, column=1, sticky=tk.W, pady=4)
	ttk.Button(master, text='Register', command=show_entry_fields).grid(row=20, column=0, sticky=tk.W, pady=4)

	master.mainloop()
	return param_dict


def set_param_dict():
	# Setting default values for the script
	return {
		'RI min': 0, 'RI max': np.inf, 'RIA min': 0, 'RIA max': np.inf, 'KD min': 0, 'KD max': np.inf, 'CI min': 0,
		'CI max': np.inf, 'AKTA baseline mode': 'linear', 'AKTA baseline deg': 2, 'AKTA_pathlength': 0.2,
		'AKTA neg val handle': 'no', 'graph width': 2, 'xaxis_title_font_size': 12, 'yaxis_title_font_size': 12,
		'xaxis_ticks_font_size': 12, 'yaxis_ticks_font_size': 12, 'plot_template': 'plotly', 'Baseline graph width': 2,
		'marker_size': 10, 'vertex_point_window': 100, 'peak_onset_var_deg': 0.5, 'panta_intermediate_plot': 'no',
		'Bmin_min': 0, 'Bmin_max': np.inf, 'Bmax_min': 0, 'Bmax_max': np.inf, 'KD_fit_min': 0, 'KD_fit_max': np.inf,
		'k_coop_min': -np.inf, 'k_coop_max': np.inf
	}


def linear_model(x, a, b):
	y = x * a + b
	return y


def model_fida_1to1(x, ri, ria, kd):
	y = (1 + (1 / kd) * x) / (((1 / ri) - (1 / ria)) + (1 + ((1 / kd) * x)) * (1 / ria))
	return y


def model_fida_excess(x, ri, ria, kd, ci):
	y = 1 / ((ci + x + kd - np.sqrt((ci + x + kd) ** 2 - 4 * ci * x)) / (2 * ria * ci) + (
			ci - x - kd + np.sqrt((ci + x + kd) ** 2 - 4 * ci * x)) / (2 * ri * ci))
	return y


def hill_equation(x, bmin, bmax, kd, k_coop):
	y = bmin + ((x ** k_coop) * (bmax - bmin)) / ((kd ** k_coop) + (x ** k_coop))
	return y


def hill_simple(x, bmin, bmax, kd):
	y = bmin + (x * (bmax - bmin)) / (kd + x)
	return y


def fit_4pl(x, bmin, bmax, kd, k_coop):
	y = bmax + (bmin - bmax) / (1 + (x / kd) ** k_coop)
	return y


def fida_text2floats(x_id, y_id, df):
	func_list = []
	func_list2 = []

	# Slice the dataframe to only include the rows that contain a concentration
	data_slice = df[df[x_id].str.contains('nM]', na=False)]
	xdata_list = data_slice[x_id].astype(str).str.replace(",", ".").tolist()
	ydata_list = pd.to_numeric(data_slice[y_id].astype(str).str.replace(",", ".")).tolist()
	for m in range(len(xdata_list)):
		conc = xdata_list[m][xdata_list[m].index('[') + 1: xdata_list[m].index(' nM]')]
		func_list.append(float(conc))

	func_list2.append(func_list)
	func_list2.append(ydata_list)
	return func_list2


def replicate_mean_error(unique_conc, redundant_conc, y_val):
	"""
	The function takes a list of unique concentrations, the redundant concentration list and a list of y values corresponding to the redundant concentration list.
	The function then calculates appropriate mean values and the standard deviation on the values used for the mean calculation.

	:param unique_conc:
	:param redundant_conc:
	:param y_val:
	:return:
	"""

	func_mean = []
	func_std = []

	y_val = list(y_val)

	for a4 in range(len(unique_conc)):
		# Getting indices of the specific concentration
		conc_indices = [k for k, x in enumerate(redundant_conc) if x == unique_conc[a4]]

		func_var = []
		for a5 in range(len(conc_indices)):
			func_var.append(y_val[conc_indices[a5]])
		func_mean.append(sum(func_var) / len(func_var))

		# If more than one replicate is present in the dataset calculate standard deviation.
		if len(conc_indices) == 1:
			func_std.append(0)
		elif len(conc_indices) > 1:
			func_std.append(statistics.stdev(func_var))

	return [func_mean, func_std]


# TODO: replace with numpy.linspace
def interval_generator(start, end):
	# The function is used to generate values of an interval with varying step size.
	# This is to have a nice distribution of data points throughout the interval without having
	# excessive number of points which make the html files large and laggy.

	counter = start
	interval = [counter]

	if counter == 0:
		counter = counter + 0.00001

	while counter < end:
		counter = counter + counter * 0.005  # Setting the step size increment
		interval.append(counter)

	return np.array(interval)


def akta_baseline(baseline, x_val, y_val, mode, param_dict):
	baseline_values = [[], []]

	for k in range(0, len(baseline)):
		baseline[k] = float(baseline[k])

	if len(baseline) == 1:
		baseline_values = [x_val, [baseline[0]] * len(x_val)]

	elif len(baseline) == 2:

		# Calculate flat baseline similar to Unicorn software
		baseline_x1 = float(x_val.iloc[(x_val - baseline[0]).abs().argsort()[:1]])
		baseline_y1 = float(y_val.iloc[(x_val - baseline[0]).abs().argsort()[:1]])
		baseline_x2 = float(x_val.iloc[(x_val - baseline[1]).abs().argsort()[:1]])
		baseline_y2 = float(y_val.iloc[(x_val - baseline[1]).abs().argsort()[:1]])

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


def color_selector(color_count, color_list):
	"""
	The function is used to set the color-determining index used in the plotting function.
	If more graphs are plotted than there are colors in the color_list then the counting will start over.

	:param color_count:
	:param color_list:
	:return:
	"""

	if color_count < len(color_list) - 1:
		color_count += 1
	else:
		color_count = 0

	return color_count


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
		# data[str(col_list[i])] = col_data

	# Check for duplicate columns
	for col in data.columns:
		col_count = data.columns.tolist().count(col)
		if col_count > 1:
			print('The column name ' + str(col) + ' appears more than once')

	return data


def plot_customize(figure, flags, ax_id):
	in_list = flags.split(';')

	for a1 in range(len(in_list)):
		if in_list[a1] == 'logx':
			figure['layout']['xaxis' + ax_id]['type'] = 'log'
			figure['layout']['xaxis' + ax_id]['dtick'] = 1

		elif in_list[a1] == 'logy':
			figure['layout']['yaxis' + ax_id]['type'] = 'log'
			figure['layout']['yaxis' + ax_id]['dtick'] = 1

	return None


def plot_func(figure, graph_name, x_val, y_val, marker, x_title, y_title, subplot_row, subplot_col, comment, i, param_dict, color_list):
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
	global used_graph_names, color_count_global

	used_graph_names.append(graph_name)
	if used_graph_names.count(graph_name) > 1 and figure != 'plot_fig':
		legend_show = False
	else:
		legend_show = True

	if comment == 'None':
		if marker.lower() in ['lines', 'line']:
			figure.add_trace(
				go.Scatter(
					name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show,
					line=dict(width=param_dict['graph width'], color=color_list[color_count_global])),
				row=subplot_row, col=subplot_col
			)

		elif marker.lower() in ['dot', 'dots', 'marker', 'markers']:
			figure.add_trace(
				go.Scatter(
					name=graph_name, x=x_val, y=y_val, mode='markers', legendgroup=graph_name, showlegend=legend_show,
					marker=dict(size=param_dict['marker_size'], color=color_list[color_count_global])),
				row=subplot_row, col=subplot_col
			)

	elif comment == 'Scatter_error':
		if marker.lower() in ['lines', 'line']:
			figure.add_trace(
				go.Scatter(
					name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show,
					error_y=dict(type='data', array=std_dev, visible=True),
					line=dict(width=param_dict['graph width'], color=color_list[color_count_global])),
				row=subplot_row, col=subplot_col
			)

		elif marker.lower() in ['dot', 'dots', 'marker', 'markers']:
			figure.add_trace(
				go.Scatter(
					name=graph_name, x=x_val, y=y_val, mode='markers', legendgroup=graph_name, showlegend=legend_show,
					error_y=dict(type='data', array=std_dev, visible=True),
					marker=dict(size=param_dict['marker_size'], color=color_list[color_count_global])),
				row=subplot_row, col=subplot_col
			)

	elif comment == 'Bar':
		figure.add_trace(
			go.Bar(name=graph_name, x=x_val, y=y_val, legendgroup=graph_name, showlegend=legend_show, marker=dict(color=color_list[color_count_global])),
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
				error_y=dict(type='data', array=error_values, visible=True),
				marker=dict(color=color_list[color_count_global])),
			row=subplot_row, col=subplot_col
		)

	elif comment == 'AKTA_baseline':
		figure.add_trace(
			go.Scatter(name=graph_name, x=x_val, y=y_val, mode='lines', legendgroup=graph_name, showlegend=legend_show, line=dict(width=param_dict['graph width'], color='rgb(255,0,0)')),
			row=subplot_row, col=subplot_col
		)

	elif comment == 'AKTA_fraction':
		figure.add_trace(
			go.Scatter(
				name=graph_name, x=x_val, y=y_val, mode='lines', fill='tozeroy', legendgroup=graph_name,
				showlegend=legend_show,
				line=dict(width=param_dict['graph width'], color=color_list[color_count_global])),
			row=subplot_row, col=subplot_col
		)

	# Determining the "number id" for the individual subplot. This can be used for targeting customization
	# such as axis labels to this subplot.
	subplot_id = ((subplot_row - 1) * subplot_col_count) + subplot_col
	if subplot_id == 1:
		ax_id = ''
	else:
		ax_id = str(subplot_id)

	figure['layout']['xaxis' + ax_id]['title'] = x_title
	figure['layout']['yaxis' + ax_id]['title'] = y_title

	figure.for_each_xaxis(lambda axis: axis.title.update(font=dict(size=param_dict['xaxis_title_font_size'])))
	figure.for_each_yaxis(lambda axis: axis.title.update(font=dict(size=param_dict['yaxis_title_font_size'])))

	figure.update_layout(
		template=param_dict['plot_template'],
		xaxis=dict(title_font=dict(size=param_dict['xaxis_title_font_size']), tickfont_size=param_dict['xaxis_ticks_font_size']),
		yaxis=dict(title_font=dict(size=param_dict['yaxis_title_font_size']), tickfont_size=param_dict['yaxis_ticks_font_size'])
	)

	if python_misc_global[i] != 'None':
		plot_customize(figure, python_misc_global[i], ax_id)


def vertex_detect(x_val, y_val, window_size):
	vertex_min = []
	vertex_max = []

	x_values = x_val.tolist()
	y_values = y_val.tolist()

	if len(x_values) != len(y_values):
		print('ERROR! Different length of arrays')
	else:
		for j in range(window_size, len(x_values) - window_size):
			window_values = y_values[j - window_size:j + window_size]
			if y_values[j] == max(window_values):
				vertex_max.append("{:.1f}".format(x_values[j]))
			elif y_values[j] == min(window_values):
				vertex_min.append("{:.1f}".format(x_values[j]))

	return [vertex_min, vertex_max]


def data_diff(x_val, y_val, savgol_window, savgol_pol, plot_intermediate, i, fig, graph_name, x_titles, y_titles, subplot_row, subplot_col, param_dict, color_list):
	"""

	:param x_val:
	:param y_val:
	:param savgol_window:
	:param savgol_pol:
	:param plot_intermediate:
	:param i:
	:param fig:
	:param graph_name:
	:param x_titles:
	:param y_titles:
	:param subplot_row:
	:param subplot_col:
	:param param_dict:
	:param color_list:
	:return:
	"""
	bin_size = int(len(y_val) / 10)
	bin_y = binned_statistic(list(x_val), list(y_val), statistic='mean', bins=bin_size)[0]  # Create bin for y values
	bin_x = binned_statistic(list(x_val), list(x_val), statistic='mean', bins=bin_size)[0]  # Create bin for x values

	# Update the bins to include first and last value of the original x_val and y_val.
	# This is to make sure the interpolation range can accomodate all values of the original data.
	bin_x = np.insert(bin_x, 0, list(x_val)[0])
	bin_x = np.insert(bin_x, len(bin_x), list(x_val)[-1])
	bin_y = np.insert(bin_y, 0, list(y_val)[0])
	bin_y = np.insert(bin_y, len(bin_y), list(y_val)[-1])

	# Create interpolation function based on the binned data
	f = interp1d(bin_x, bin_y)

	# Create new x_val and y_val using the interpolation function. This is basically to fill data between the bin
	# points i.e. making data big again.
	# I also create new x values that are evenly spaced in the original x value range
	x_val_new = np.linspace(list(x_val)[0], list(x_val)[-1], len(x_val))
	y_val_new = f(list(x_val_new))

	# Smooth the interpolated data using savgol filter and take first derivative.
	y_val_filter = savgol_filter(y_val_new, savgol_window, savgol_pol, deriv=1)

	# Smooth the first derivative data and take the second derivative.
	y_val_filter2 = savgol_filter(y_val_filter, savgol_window, savgol_pol, deriv=1)

	if plot_intermediate == 'yes':
		plot_func(fig, graph_name + '_bin', bin_x, bin_y, 'line', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)
		plot_func(fig, graph_name + '_interpolated', x_val_new, y_val_new, 'line', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)
		plot_func(fig, graph_name + '_1st deriv', x_val_new, y_val_filter, 'line', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)
		plot_func(fig, graph_name + '_2nd deriv', x_val_new, y_val_filter2, 'line', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)

	return [x_val_new, y_val_new, y_val_filter, y_val_filter2]


def auto_baseline(x_val, y_val, y_val_diff2, plot_baseline, fig, graph_name, x_titles, y_titles, subplot_row, subplot_col, i, param_dict, color_list):
	# Converting x values and 2nd derivative data to lists.
	x_val_list = list(x_val)
	y_val_list = list(y_val_diff2)

	# Determine baseline start as the point on the 2nd derivative data where the following
	# 10 data points are all within 1% of the max value of that derivative i.e. very close to zero
	for a2 in range(len(y_val_list)):
		data_window = y_val_list[a2:a2 + 10]
		counter = 0
		for a3 in range(len(data_window)):
			if abs(data_window[a3]) < 0.03 * max(y_val_diff2):
				counter = counter + 1
			else:
				# TODO: This does nothing...remove?
				pass

		if counter == 10:
			baseline_start_idx = a2
			break

		else:
			baseline_start_idx = 0

	for a4 in range(baseline_start_idx, len(y_val_list)):
		data_window = y_val_list[a4:a4 + 10]
		counter = 0
		for a5 in range(len(data_window)):
			if abs(data_window[a5]) > 0.1 * max(
					y_val_diff2) and a4 > baseline_start_idx:  # Make sure the index of baseline end is bigger than the start index
				counter = counter + 1
			else:
				# TODO: This does nothing...remove?
				pass
		if counter == 10:
			baseline_end_idx = a4
			break

		else:
			baseline_end_idx = len(y_val_list) - 1

	baseline_x = x_val_list[baseline_start_idx:baseline_end_idx]
	baseline_y = list(y_val)[baseline_start_idx:baseline_end_idx]

	parameters, covariance = curve_fit(linear_model, baseline_x, baseline_y, method='trf', maxfev=10000)

	# y_pred = linear_model(np.asarray(baseline_x), parameters[0], parameters[1])
	# r_square = r2_score(baseline_y, y_pred)
	# print('R square: ' + str(r_square))

	y_pred = linear_model(x_val, parameters[0], parameters[1])

	if plot_baseline == 'yes':
		plot_func(fig, graph_name + '_auto_baseline', x_val, y_pred, 'line', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)

	# Return coefficients for the linear fit
	return [parameters[0], parameters[1], baseline_x[-1]]


def onset_detect(x_val, y_val, slope, intercept, detect_limit, baseline_cutoff):
	x_val_list = list(x_val)
	y_val_list = list(y_val)

	for a6 in range(len(y_val_list)):
		data_obs = y_val_list[a6]
		data_calc = linear_model(x_val_list[a6], slope, intercept)

		# The onset is defined as the point where the observed data deviates more and 0.5% from the linear baseline.
		# It is further required that the x value is higher than the baseline cutoff to avoid predicting onsets on
		# noise in the beginning of the dataset.
		if data_obs > data_calc * (1 + (detect_limit / 100)) and x_val_list[a6] > baseline_cutoff:
			onset_val = x_val_list[a6]
			onset_val = "{:.1f}".format(onset_val)
			break

		else:
			onset_val = 0

	return onset_val


def ip_detect(x_val, y_val, window_size, onset_cutoff):
	func_list = []

	func_vertex = vertex_detect(x_val, y_val, window_size)

	func_ip = func_vertex[1]

	# Iterate through the calculated inflection points and only include those that are above the onset
	for a7 in range(len(func_ip)):
		if float(func_ip[a7]) > float(onset_cutoff):
			func_list.append(func_ip[a7])

	return func_list


def table_plot(figure, col_names_list, col_values_list, ID_list):
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

	table_pos_list = []

	# Defining a dataframe with the first column
	data_func = {col_names_list[0]: col_values_list[0]}
	df_func = pd.DataFrame(data=data_func)

	# Go over the lists loaded into the function and append them to the dataframe.
	for b1 in range(1, len(col_names_list)):
		df_func[col_names_list[b1]] = col_values_list[b1]

	# Check if the number of samples given by user is the same as given to the function. Deviations can occur e.g. in
	# global fittings where the code is generating extra traces.
	if len(ID_list) == len(col_values_list[0]):
		# Code for extracting the user defined table coordinates i.e. what line in the table the data should be inputted on.
		# First loop iterates through ID_list i.e. the samples in the data sheet and extracts the flags separated by semi colon.
		# The code goes over each flag and if table_pos= flag is present it extracts the integer following the equal sign.

		for b2 in range(len(ID_list)):
			func_var = python_misc_global[b2].split(';')

			# If user has specified table_pos in the python_misc column the entry will be extracted into dummy list.
			func_list = []
			for b3 in range(len(func_var)):
				if 'table_pos=' in func_var[b3]:
					func_list.append(func_var[b3])

			# If no table position has been specified the code will just use the row index from excel sheet.
			if len(func_list) == 0:
				table_pos_list.append(b2 + 1)
			else:
				func_coord = func_var[b3].split('=')[1]
				func_coord = int(func_coord)
				table_pos_list.append(func_coord)

		# Add the table position list to the dataframe
		df_func['table_coords'] = table_pos_list

		# Sort the dataframe by the table coordiates to get the desired order
		df_func2 = df_func.sort_values(by='table_coords', ascending=True)
		del df_func2['table_coords']

		figure.add_trace(
			go.Table(
				header=dict(values=list(df_func2.columns), align='left'),
				cells=dict(values=df_func2.transpose().values.tolist(), align='left', height=50))
		)

	else:
		figure.add_trace(
			go.Table(
				header=dict(values=list(df_func.columns), align='left'),
				cells=dict(values=df_func.transpose().values.tolist(), align='left', height=50))
		)


def scatter_fit(model_function, model_name, x_val, y_val, fitting_mode, fitting_min, fitting_max, i, ID_list,
				sample_notes, fit_models, fig, graph_name, x_titles, y_titles, subplot_row, subplot_col, param_dict,
				color_list):
	out_dict = {}

	x_val_list = []
	y_val_list = []

	unique_x = list(dict.fromkeys(x_val))

	if model_name == 'Hill':
		bound_param = [
			(param_dict['Bmin_min'], param_dict['Bmax_min'], param_dict['KD_fit_min'], param_dict['k_coop_min']),
			(param_dict['Bmin_max'], param_dict['Bmax_max'], param_dict['KD_fit_max'], param_dict['k_coop_max'])]

	elif model_name == 'Hill_simple':
		bound_param = [
			(param_dict['Bmin_min'], param_dict['Bmax_min'], param_dict['KD_fit_min']),
			(param_dict['Bmin_max'], param_dict['Bmax_max'], param_dict['KD_fit_max'])
		]

	elif model_name == '4PL':
		bound_param = [
			(param_dict['Bmin_min'], param_dict['Bmax_min'], param_dict['KD_fit_min'], param_dict['k_coop_min']),
			(param_dict['Bmin_max'], param_dict['Bmax_max'], param_dict['KD_fit_max'], param_dict['k_coop_max'])
		]

	elif model_name == 'FIDA_1to1':
		bound_param = [
			(param_dict['RI min'], param_dict['RIA min'], param_dict['KD min']),
			(param_dict['RI max'], param_dict['RIA max'], param_dict['KD max'])
		]

	elif model_name == 'FIDA_excess':
		bound_param = (
			(param_dict['RI min'], param_dict['RIA min'], param_dict['KD min'], param_dict['CI min']),
			(param_dict['RI max'], param_dict['RIA max'], param_dict['KD max'], param_dict['CI max'])
		)

	if fitting_mode == 'Local':

		# In local fitting mode the list of unique x values is the only x list needed.
		x_val_list.append(unique_x)

		# Calculate mean and std error on replicate values in the data
		y_mean = replicate_mean_error(unique_x, x_val, y_val)[0]
		y_val_list.append(y_mean)
		y_mean_err = replicate_mean_error(unique_x, x_val, y_val)[1]
		out_dict['y_mean'] = y_mean

	elif fitting_mode == 'Global':

		# Get the max number of times a concentration value occurs in raw data. This is the number of graphs to calculate.
		func_conc_count = []
		for m in unique_x:
			func_conc_count.append(list(x_val).count(m))
		graph_number = max(func_conc_count)

		# Create list of empty lists for filling with the individual graphs.
		x_val_list = [[] for _ in range(graph_number)]
		y_val_list = [[] for _ in range(graph_number)]

		# Sort the data values into the appropriate lists
		for j in unique_x:
			conc_indices = [k for k, x in enumerate(list(x_val)) if
							x == j]  # Getting indices of the specific concentration
			for k in range(len(conc_indices)):
				x_val_list[k].append(list(x_val)[conc_indices[k]])
				y_val_list[k].append(list(y_val)[conc_indices[k]])

	func_KD_list = []
	func_R2_list = []

	for c in range(len(x_val_list)):
		parameters, covariance = curve_fit(model_function, list(x_val_list[c]), list(y_val_list[c]), bounds=bound_param, method='trf', maxfev=10000)

		if len(x_val_list) == 1:
			master_dict['ID_list_new'].append(str(ID_list[i]))
			master_dict['notes_list'].append(sample_notes[i])
			master_dict['model_list'].append(fit_models[i] + ', ' + fitting_mode)

		elif len(x_val_list) > 1:
			master_dict['ID_list_new'].append(str(ID_list[i]) + '_fit' + str(c + 1))
			master_dict['notes_list'].append(' ')
			master_dict['model_list'].append(' ')

		if model_name == 'Hill':

			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = hill_equation(x_fit, parameters[0], parameters[1], parameters[2], parameters[3])
			plot_func(fig, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'line', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)

			# Calcularing R squared value
			y_fit_small = hill_equation(x_val_list[c], parameters[0], parameters[1], parameters[2], parameters[3])
			method_name(c, func_KD_list, func_R2_list, parameters, y_fit_small, y_val_list)

		elif model_name == 'Hill_simple':
			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = hill_simple(x_fit, parameters[0], parameters[1], parameters[2])
			plot_func(fig, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'line', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)

			# Calcularing R squared value
			y_fit_small = hill_simple(np.asarray(x_val_list[c]), parameters[0], parameters[1], parameters[2])
			r_square = r2_score(y_val_list[c], y_fit_small)
			master_dict['R_square'].append("{:.3f}".format(r_square))
			func_R2_list.append(r_square)

			func_Bmin = "{:.3f}".format(parameters[0])
			func_Bmax = "{:.3f}".format(parameters[1])
			func_KD = "{:.3f}".format(parameters[2])
			func_KD_list.append(parameters[2])

			master_dict['KD_fit'].append(func_KD)
			master_dict['fit_parameters'].append('Bmin=' + str(func_Bmin) + ', Bmax=' + str(func_Bmax))

		elif model_name == '4PL':
			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = fit_4pl(x_fit, parameters[0], parameters[1], parameters[2], parameters[3])
			plot_func(fig, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'line', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)

			# Calcularing R squared value
			y_fit_small = fit_4pl(x_val_list[c], parameters[0], parameters[1], parameters[2], parameters[3])
			method_name(c, func_KD_list, func_R2_list, parameters, y_fit_small, y_val_list)

		elif model_name == 'FIDA_1to1':
			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = model_fida_1to1(x_fit, parameters[0], parameters[1], parameters[2])
			plot_func(fig, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'line', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)

			# Calcularing R squared value
			y_fit_small = model_fida_1to1(np.asarray(x_val_list[c]), parameters[0], parameters[1], parameters[2])
			r_square = r2_score(y_val_list[c], y_fit_small)
			master_dict['R_square'].append("{:.3f}".format(r_square))
			func_R2_list.append(r_square)

			func_RI = "{:.3f}".format(parameters[0])
			func_RIA = "{:.3f}".format(parameters[1])
			func_KD = "{:.3f}".format(parameters[2])
			func_KD_list.append(parameters[2])

			master_dict['KD_fit'].append(func_KD)
			master_dict['fit_parameters'].append('RI=' + str(func_RI) + ', RIA=' + str(func_RIA))

		elif model_name == 'FIDA_excess':

			# Generating a fitting curve with many points for the plot
			x_fit = interval_generator(fitting_min, fitting_max)
			y_fit = model_fida_excess(x_fit, parameters[0], parameters[1], parameters[2], parameters[3])
			plot_func(fig, graph_name + '_fit' + str(c + 1), x_fit, y_fit, 'line', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)

			# Calcularing R squared value
			y_fit_small = model_fida_excess(np.asarray(x_val_list[c]), parameters[0], parameters[1], parameters[2], parameters[3])
			r_square = r2_score(y_val_list[c], y_fit_small)
			master_dict['R_square'].append("{:.3f}".format(r_square))
			func_R2_list.append(r_square)

			func_RI = "{:.3f}".format(parameters[0])
			func_RIA = "{:.3f}".format(parameters[1])
			func_KD = "{:.3f}".format(parameters[2])
			func_KD_list.append(parameters[2])
			func_CI = "{:.3f}".format(parameters[3])

			master_dict['KD_fit'].append(func_KD)
			master_dict['fit_parameters'].append(
				'RI=' + str(func_RI) + ', RIA=' + str(func_RIA) + ', CI=' + str(func_CI))

	# If global fitting (i.e. more than one KD and R^2 has been calculated for the data) add final result to table
	if fitting_mode == 'Global':
		master_dict['ID_list_new'].append(str(ID_list[i]))
		master_dict['notes_list'].append(sample_notes[i])
		master_dict['model_list'].append(fit_models[i] + ', ' + fitting_mode)
		master_dict['fit_parameters'].append(' ')

		global_KD = statistics.mean(func_KD_list)
		global_R2 = statistics.stdev(func_KD_list)
		master_dict['KD_fit'].append("{:.3f}".format(global_KD) + ' ' + u"\u00B1" + ' ' + str("{:.3f}".format(global_R2)))
		master_dict['R_square'].append(' ')

	# TODO: out_dict is never filled with anything...
	return out_dict


def method_name(c, func_KD_list, func_R2_list, parameters, y_fit_small, y_val_list):
	r_square = r2_score(y_val_list[c], y_fit_small)
	master_dict['R_square'].append("{:.3f}".format(r_square))
	func_R2_list.append(r_square)
	func_Bmin = "{:.3f}".format(parameters[0])
	func_Bmax = "{:.3f}".format(parameters[1])
	func_KD = "{:.3f}".format(parameters[2])
	func_KD_list.append(parameters[2])
	func_k_coop = "{:.3f}".format(parameters[3])
	master_dict['KD_fit'].append(func_KD)
	master_dict['fit_parameters'].append('Bmin=' + str(func_Bmin) + ', Bmax=' + str(func_Bmax) + ', k_coop=' + str(func_k_coop))


def main_func(args, param_dict):
	global python_misc_global, color_count_global, master_dict, std_dev, error_values, subplot_col_count

	# Load data from excel sheet and prep it for analysis
	print('Loading and cleaning your data')
	df = pd.concat(pd.read_excel(args.input_file, sheet_name=None), ignore_index=True)
	df = data_clean(df)

	# Load user input from the excel sheet
	ID_list = df['IDs'].tolist()
	ID_list = [x for x in ID_list if str(x) != 'nan']

	# Loading a list with colors for plotting. The colors come from https://plotly.com/python/discrete-color/.
	color_list = ['#1F77B4', '#FF7F0E', '#2CA02C', '#9467BD', '#FECB52', '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF']
	# If the number of samples equals the length of the color list the script will append one extra color to avoid the same
	# color comparisons in subplots.
	if len(ID_list) % len(color_list) == 0:
		color_list.append('#B6E880')

	# Get the data intervals for slicing the data to only analyse sections of full dataset.
	data_interval = df['Data_interval'].fillna(0).replace(',', '.').tolist()
	data_interval = data_interval[0:len(ID_list)]
	plot_markers = df['Plot_markers'].fillna('line').tolist()
	sample_notes = df['Notes'].fillna('None').tolist()
	fit_models = df['Fit_model'].fillna('None').tolist()
	fit_modes = df['Fit_approach'].fillna('Local').tolist()
	fit_intervals = df['Fitting_interval'].fillna('0;0').tolist()

	python_misc_global = df['Python_misc'].fillna('None').tolist()

	axis_title = df['Axis_titles'].fillna(';').tolist()
	x_titles = []
	y_titles = []
	for id_b in range(len(axis_title)):
		x_titles.append(axis_title[id_b].split(';')[0])
		y_titles.append(axis_title[id_b].split(';')[1])

	# Loading a dataframe for outputting
	pd_out = pd.DataFrame({'IDs': ID_list})

	# Setting up the figure for plotting.
	subplot_row = []
	subplot_col = []
	subplot_coord = df['Sub_plot'].fillna('1;1').tolist()
	for id_a in range(len(subplot_coord)):
		# The first subplot coordinate is appended to subplot row list as integer
		subplot_row.append(int(subplot_coord[id_a].split(';')[0]))
		# The first subplot coordinate is appended to subplot row list as integer
		subplot_col.append(int(subplot_coord[id_a].split(';')[1]))
	subplot_row_count = max(subplot_row)
	subplot_col_count = max(subplot_col)
	fig = make_subplots(rows=subplot_row_count, cols=subplot_col_count, vertical_spacing=0.07, horizontal_spacing=0.07)
	plot_fig = make_subplots(rows=subplot_row_count, cols=subplot_col_count, vertical_spacing=0.07, horizontal_spacing=0.07)

	if args.plot_type in ['AKTA', 'akta', 'Akta']:
		sample_areas_tot = []
		baseline_area_tot = []
		fraction_retentions = []
		fraction_areas = []
		fraction_baseline = []
		fraction_calculation = []
		fraction_concentrations = []
		fraction_yield = []
		culture_yield = []

		# Get listed extinction coefficient from excel sheet and replace empty values with zero.
		# Also make sure the the right decimal seperator is used in case people use e.g. Danish excel
		ext_coeff = df['AKTA_extinc_coeff'].fillna(0).replace(',', '.').tolist()
		ext_coeff = ext_coeff[0:len(ID_list)]  # Restrict the list so it has the same length as ID list

		# Get the volumes loaded. This is the total volume of supernatant loaded during experiments.
		# Replace empty values with zero
		volume_load = df['AKTA_volume_load'].fillna(0).replace(',', '.').tolist()
		volume_load = volume_load[0:len(ID_list)]  # Restrict the list so it has the same length as ID list

		for i in range(len(ID_list)):
			print('Running data from: ' + str(ID_list[i]) + ' (' + str(sample_notes[i]) + ')')

			x_id = 'x' + str(i + 1)
			y_id = 'y' + str(i + 1)

			graph_name = ID_list[i]
			master_dict['notes_list'].append(sample_notes[i])

			# Extracting the raw x and y values from the excel sheet
			xs = df[x_id][pd.to_numeric(df[x_id], errors='coerce').notnull()]
			ys = df[y_id][pd.to_numeric(df[y_id], errors='coerce').notnull()]

			# Slicing the data so only data within the specified data interval is included.
			if data_interval[i] != 0:
				# interval_var = list(data_interval[i].split(';'))
				interval_var = data_interval[i].split(';')
				xsys_interval = pd.concat([xs, ys], axis=1)
				xsys_interval_slice = xsys_interval[
					(xsys_interval[x_id] >= float(interval_var[0])) & (xsys_interval[x_id] <= float(interval_var[1]))]
				xs = xsys_interval_slice[x_id]
				ys = xsys_interval_slice[y_id]

			# If specified by user all the negative y values will be replaced by zeros. Can be relevant since the negative
			# values negatively impact the area calculations
			if param_dict['AKTA neg val handle'] == 'yes':
				ys[ys < 0] = 0

			if str(df['AKTA_baseline'][i]).count(';') == 0:
				baseline_coord = [float(str(df['AKTA_baseline'][i]).replace(',', '.'))]
			elif df['AKTA_baseline'][i].count(';') > 0:
				baseline_coord = list(df['AKTA_baseline'][i].replace(',', '.').split(';'))

			# NOTE THIS IS THE PLACE TO CHANGE THE METHOD FOR BASELINE CALLING. THE MODE CAN BE 'peakutils' OR 'linear'.
			baseline_values = akta_baseline(baseline_coord, xs, ys, param_dict['AKTA baseline mode'], param_dict)

			AUC_sample_tot = auc(xs, ys)  # Calculate the total AUC from entire sample
			sample_areas_tot.append("{:.3f}".format(AUC_sample_tot))

			AUC_baseline_tot = auc(baseline_values[0], baseline_values[1])
			baseline_area_tot.append("{:.3f}".format(AUC_baseline_tot))

			plot_func(fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)
			plot_func(fig, graph_name + '_baseline', baseline_values[0], baseline_values[1], 'lines', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'AKTA_baseline', i, param_dict, color_list)
			plot_func(plot_fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)
			plot_func(plot_fig, graph_name + '_baseline', baseline_values[0], baseline_values[1], plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'AKTA_baseline', i, param_dict, color_list)

			# Make section in case fraction is not included in analysis
			if str(df['AKTA_fraction'][i]).count(';') == 0:
				fraction_areas.append('N/A')
				fraction_baseline.append('N/A')
				fraction_concentrations.append('N/A')
				fraction_calculation.append('N/A')
				fraction_yield.append('N/A')
				culture_yield.append('N/A')
				fraction_retentions.append('N/A')

			# Make section for data where sample has been fractionated e.g. for elution
			elif str(df['AKTA_fraction'][i]).count(';') == 1:

				# Get fraction boundaries from excel and turn them to floats
				fraction = list(df['AKTA_fraction'][i].replace(',', '.').split(';'))
				for j in range(0, len(fraction)):
					fraction[j] = float(fraction[j])
				frac_volume = abs(fraction[1] - fraction[0])

				# Slicing AKTA data to only get data points that are part of the fraction interval
				xsys = pd.concat([xs, ys], axis=1)
				xsys_slice = xsys[(xsys[x_id] >= fraction[0]) & (xsys[x_id] <= fraction[1])]

				# Slicing baseline data to only get data points that are part of fraction interval
				# Create a dataframe from the lists containing the baseline values
				xsys_baseline = pd.DataFrame({str(x_id): baseline_values[0], str(y_id): baseline_values[1]})
				xsys_baseline_slice = xsys_baseline[(xsys_baseline[x_id] >= fraction[0]) & (xsys_baseline[x_id] <= fraction[1])]

				# The retention time/volume is defined as the x value where the y value hits max
				Retention_frac = float(xsys_slice[x_id][xsys_slice[y_id] == max(xsys[y_id])])
				fraction_retentions.append("{:.3f}".format(Retention_frac))

				AUC_frac = auc(xsys_slice[x_id], xsys_slice[y_id])
				fraction_areas.append("{:.3f}".format(AUC_frac))

				AUC_frac_baseline = auc(xsys_baseline_slice[x_id], xsys_baseline_slice[y_id])
				fraction_baseline.append("{:.3f}".format(AUC_frac_baseline))

				AUC_calculation = AUC_frac - AUC_frac_baseline
				fraction_calculation.append("{:.3f}".format(AUC_calculation))

				frac_yield = AUC_calculation / (float(ext_coeff[i]) * param_dict['AKTA_pathlength'])
				frac_yield = frac_yield / 1000  # Divide by 1000 to get the yield in mg
				fraction_yield.append("{:.3f}".format(frac_yield))

				fraction_conc = frac_yield / frac_volume
				fraction_concentrations.append("{:.3f}".format(fraction_conc))

				cul_yield = (frac_yield / volume_load[i] * 1000)
				culture_yield.append("{:.3f}".format(cul_yield))

				plot_func(fig, graph_name + '_fraction', xsys_slice[x_id], xsys_slice[y_id], 'lines', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'AKTA_fraction', i, param_dict, color_list)
				plot_func(plot_fig, graph_name + '_fraction', xsys_slice[x_id], xsys_slice[y_id], 'lines', x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'AKTA_fraction', i, param_dict, color_list)

				# Adding the raw data generated in code above to a dataframe that can be outputted to look at the raw data.
				# In the concattenated dataframes the index is reset because the dataframes do not have the same length
				pd_out = pd.concat([
					pd_out,
					pd.DataFrame({x_id: xs, y_id: ys}).reset_index(drop=True),
					pd.DataFrame(
						{x_id + '_baseline': baseline_values[0], y_id + '_baseline': baseline_values[1]}).reset_index(
						drop=True),
					pd.DataFrame(
						{x_id + '_fraction': xsys_slice[x_id], y_id + '_fraction': xsys_slice[y_id]}).reset_index(
						drop=True)],
					ignore_index=False, axis=1
				)
			# Update color value using color_selector function
			color_count_global = color_selector(color_count_global, color_list)

		fig.add_trace(
			go.Table(
				header=dict(values=[
					'Sample ID', 'Sample notes', 'Total sample area', 'Total baseline area',
					'Retention time/volume (beta)', 'Area of fraction', 'Baseline area (fraction)',
					'Area used for calculation', 'Concentration in fraction [mg/mL]', 'Fraction yield [mg]',
					'Culture yield [ug/mL]'],
					align='left'),
				cells=dict(values=[
					ID_list, master_dict['notes_list'], sample_areas_tot, baseline_area_tot, fraction_retentions,
					fraction_areas, fraction_baseline, fraction_calculation, fraction_concentrations, fraction_yield,
					culture_yield],
					align='left', height=50))
		)

		# Add the values from the plotly table to the output table
		pd_out = pd.concat([
			pd_out,
			pd.DataFrame({
				'Sample ID': ID_list, 'Sample notes': master_dict['notes_list'],
				'Total sample area': sample_areas_tot, 'Total baseline area': baseline_area_tot,
				'Area of fraction': fraction_areas, 'Baseline area (fraction)': fraction_baseline,
				'Area used for calculation': fraction_calculation,
				'Concentration in fraction [mg/mL]': fraction_concentrations,
				'Fraction yield [mg]': fraction_yield, 'Culture yield [ug/mL]': culture_yield}).reset_index(drop=True)
		], ignore_index=False, axis=1)

		# Add dropdown
		fig.update_layout(
			updatemenus=[
				dict(
					type="buttons",
					direction="left",
					buttons=list([dict(args=["type", "scatter"], label="Graphs", method="restyle"), dict(args=["type", "table"], label="Stats", method="restyle")]),
					pad={"r": 10, "t": 10}, showactive=True, x=0.11, xanchor="left", y=1.1, yanchor="top"), ])

	elif args.plot_type in ['FIDA', 'Fida', 'fida']:

		for i in range(len(ID_list)):

			print('Analysing data: ' + str(ID_list[i]))

			x_id = 'x' + str(i + 1)
			y_id = 'y' + str(i + 1)

			graph_name = ID_list[i]

			xs = fida_text2floats(x_id, y_id, df)[0]  # Extract concentration float values from the text FIDA output
			ys = fida_text2floats(x_id, y_id, df)[1]  # Extract appropriate y values from the FIDA data
			unique_conc = list(dict.fromkeys( xs))  # Create a list with unique concentration values. These should include ALL relevant x values.

			# Slicing the data so only data within the specified data interval is included.
			if data_interval[i] != 0:
				interval_var = list(data_interval[i].split(';'))
				xsys_interval = pd.concat([xs, ys], axis=1)
				xsys_interval_slice = xsys_interval[(xsys_interval[x_id] >= float(interval_var[0])) & (xsys_interval[x_id] <= float(interval_var[1]))]
				xs = xsys_interval_slice[x_id]
				ys = xsys_interval_slice[y_id]

			# Check if the data should just be plotted or if user is trying to fit
			if fit_intervals[i] != '0;0' and fit_models[i] != 'None':

				# If local fitting approach the data will be averaged and plotted with error bars
				if fit_modes[i] == 'Local':
					y_mean = replicate_mean_error(unique_conc, xs, ys)[0]
					std_dev = replicate_mean_error(unique_conc, xs, ys)[1]
					plot_func(fig, graph_name, unique_conc, y_mean, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'Scatter_error', i, param_dict, color_list)
					plot_func(plot_fig, graph_name, unique_conc, y_mean, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'Scatter_error', i, param_dict, color_list)
				else:
					plot_func(fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)
					plot_func(plot_fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)

				# Extract boundaries for the fit from the excel sheet and turn to list of floats. If no interval is supplied code will automaticallt assign 0;0.
				fitting_interval = fit_intervals[i].split(';')
				for j in range(0, len(fitting_interval)):
					fitting_interval[j] = float(
						fitting_interval[j].replace(",", "."))  # Replace commas with proper dots to get proper numbers.

				if fit_models[i].lower() == 'fida_1to1' or fit_models[i].lower() == '1to1':
					# Running fitting function
					scatter_fit(
						model_fida_1to1, 'FIDA_1to1', xs, ys, fit_modes[i], fitting_interval[0], fitting_interval[1], i,
						ID_list, sample_notes, fit_models, fig, graph_name, x_titles, y_titles, subplot_row,
						subplot_col, param_dict, color_list)
				elif fit_models[i].lower() == 'fida_excess' or fit_models[i].lower() == 'excess':
					# Running fitting function
					scatter_fit(
						model_fida_excess, 'FIDA_excess', xs, ys, fit_modes[i], fitting_interval[0],
						fitting_interval[1], i, ID_list, sample_notes, fit_models, fig, graph_name, x_titles,
						y_titles, subplot_row, subplot_col, param_dict, color_list)

				color_count_global = color_selector(color_count_global, color_list)  # Update color value using color_selector function

			else:
				master_dict['ID_list_new'].append(ID_list[i])
				master_dict['notes_list'].append(sample_notes[i])
				master_dict['model_list'].append(' ')
				master_dict['fit_parameters'].append(' ')
				master_dict['R_square'].append(' ')
				master_dict['KD_fit'].append(' ')

				plot_func(fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)
				plot_func(plot_fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)
				color_count_global = color_selector(color_count_global, color_list)  # Update color value using color_selector function

		table_plot(
			fig,
			['Samples', 'Sample notes', 'Fitting models', 'Fitted KD/EC50', 'R^2', 'Fit parameters'],
			[master_dict['ID_list_new'], master_dict['notes_list'], master_dict['model_list'], master_dict['KD_fit'],
			 master_dict['R_square'], master_dict['fit_parameters']], ID_list
		)

		# Add dropdown
		fig.update_layout(
			updatemenus=[
				dict(
					type="buttons",
					direction="left",
					buttons=list([dict(args=["type", "scatter"], label="Graphs", method="restyle"), dict(args=["type", "table"], label="Stats", method="restyle")]),
					pad={"r": 10, "t": 10}, showactive=True, x=0.11, xanchor="left", y=1.1, yanchor="top")])

	elif args.plot_type in ['Panta', 'panta', 'PANTA']:

		for i in range(len(ID_list)):

			master_dict['notes_list'].append(sample_notes[i])

			print('Running data from: ' + str(ID_list[i]) + ' (' + str(sample_notes[i]) + ')')

			x_id = 'x' + str(i + 1)
			y_id = 'y' + str(i + 1)
			graph_name = ID_list[i]

			# Extracting the raw x and y values from the excel sheet
			xs = df[x_id][pd.to_numeric(df[x_id], errors='coerce').notnull()]
			ys = df[y_id][pd.to_numeric(df[y_id], errors='coerce').notnull()]

			# Slicing the data so only data within the specified data interval is included.
			if data_interval[i] != 0:
				interval_var = list(data_interval[i].split(';'))
				xsys_interval = pd.concat([xs, ys], axis=1)
				xsys_interval_slice = xsys_interval[
					(xsys_interval[x_id] >= float(interval_var[0])) & (xsys_interval[x_id] <= float(interval_var[1]))]
				xs = xsys_interval_slice[x_id]
				ys = xsys_interval_slice[y_id]

			plot_func(fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)
			plot_func(plot_fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'None', i, param_dict, color_list)

			# Smoothing data and getting first and second derivative of the data. The last argument enables plotting of bin data,
			# interpolated data, first derivative and second derivative ('yes')
			data_diff_out = data_diff(xs, ys, 81, 2, param_dict['panta_intermediate_plot'], i, fig, graph_name, x_titles, y_titles, subplot_row, subplot_col, param_dict, color_list)
			xs_new = data_diff_out[0]
			ys_interpol = data_diff_out[1]
			ys_diff1 = data_diff_out[2]
			ys_diff2 = data_diff_out[3]

			# Creating a linear baseline based on the second derivative. The last argument enables plotting of the baseline ('yes').
			auto_baseline_coeff = auto_baseline(xs_new, ys_interpol, ys_diff2, param_dict['panta_intermediate_plot'],
												fig, graph_name, x_titles, y_titles, subplot_row, subplot_col, i,
												param_dict, color_list)

			# Determine the onset of the peak
			onset = onset_detect(xs_new, ys_interpol, auto_baseline_coeff[0], auto_baseline_coeff[1], param_dict['peak_onset_var_deg'], auto_baseline_coeff[2])
			master_dict['peak_onset'].append(onset)

			# Determine inflection points as vertex points on the first derivative data.
			ip_points = ip_detect(xs_new, ys_diff1, param_dict['vertex_point_window'], onset)
			master_dict['inflection_points'].append(ip_points)

			# The color counter is placed in the end of Panta code snippet since some of the functions above also potentially include plotting
			color_count_global = color_selector(color_count_global, color_list)  # Update color value using color_selector function

			# Determine vertex points
			vertex_out = vertex_detect(xs, ys, param_dict['vertex_point_window'])
			master_dict['vertex_min'].append(vertex_out[0])
			master_dict['vertex_max'].append(vertex_out[1])

		table_plot(
			fig,
			['Sample ID', 'Sample notes', 'Vertex points min (beta)', 'Vertex points max (beta)', 'Peak onset (beta)',
			 'Inflection points (beta)'],
			[
				ID_list, master_dict['notes_list'], master_dict['vertex_min'], master_dict['vertex_max'],
				master_dict['peak_onset'], master_dict['inflection_points']
			],
			ID_list
		)

		# Add dropdown
		fig.update_layout(
			updatemenus=[
				dict(
					type="buttons",
					direction="left",
					buttons=list([dict(args=["type", "scatter"], label="Graphs", method="restyle"),
								  dict(args=["type", "table"], label="Stats", method="restyle")]),
					pad={"r": 10, "t": 10}, showactive=True, x=0.11, xanchor="left", y=1.1, yanchor="top"), ])

	elif args.plot_type in ['Bar', 'bar', 'Bar_group', 'bar_group']:

		x_labels = df['x1']

		for i in range(len(ID_list)):

			print('Analysing data: ' + str(ID_list[i]))

			x_id = 'x' + str(i + 1)
			y_id = 'y' + str(i + 1)
			error_id = 'error' + str(i + 1)

			Group_name = ID_list[i]

			y_values = df[y_id]
			y_values = y_values.replace(np.nan, 0)

			if error_id in df.columns:  # Check if error column exists
				error_values = df[error_id]
				error_values = error_values.replace(np.nan, 0)
				plot_func(fig, Group_name, x_labels, y_values, 'N/A', x_titles[i], y_titles[i], subplot_row[i],
						  subplot_col[i], 'Bar_error', i, param_dict, color_list)
				plot_func(plot_fig, Group_name, x_labels, y_values, 'N/A', x_titles[i], y_titles[i], subplot_row[i],
						  subplot_col[i], 'Bar_error', i, param_dict, color_list)

			else:
				plot_func(fig, Group_name, x_labels, y_values, 'N/A', x_titles[i], y_titles[i], subplot_row[i],
						  subplot_col[i], 'Bar', i, param_dict, color_list)
				plot_func(plot_fig, Group_name, x_labels, y_values, 'N/A', x_titles[i], y_titles[i], subplot_row[i],
						  subplot_col[i], 'Bar', i, param_dict, color_list)

			color_count_global = color_selector(color_count_global,
												color_list)  # Update color value using color_selector function

		fig.update_layout(bargap=0.3, bargroupgap=0.05)
		plot_fig.update_layout(bargap=0.3, bargroupgap=0.05)

	elif args.plot_type in ['Scatter', 'scatter']:

		for i in range(len(ID_list)):

			print('Analysing data: ' + str(ID_list[i]))

			x_id = 'x' + str(i + 1)
			y_id = 'y' + str(i + 1)

			graph_name = ID_list[i]

			# Extracting the raw x and y values from the excel sheet.
			# Replace commas with dots in case user has e.g. Danish Excel.
			xs = df[x_id][pd.to_numeric(df[x_id], errors='coerce').notnull()]
			ys = df[y_id][pd.to_numeric(df[y_id], errors='coerce').notnull()]
			xs = pd.to_numeric(xs.astype(str).str.replace(",", "."))
			ys = pd.to_numeric(ys.astype(str).str.replace(",", "."))
			unique_conc = list(dict.fromkeys(
				xs))  # Create a list with unique concentration values. These should include ALL relevant x values.

			# Slicing the data so only data within the specified data interval is included.
			if data_interval[i] != 0:
				interval_var = list(data_interval[i].split(';'))
				xsys_interval = pd.concat([xs, ys], axis=1)
				xsys_interval_slice = xsys_interval[
					(xsys_interval[x_id] >= float(interval_var[0])) & (xsys_interval[x_id] <= float(interval_var[1]))]
				xs = xsys_interval_slice[x_id]
				ys = xsys_interval_slice[y_id]

			# Check if the data should just be plotted or if user is trying to fit
			if fit_intervals[i] != '0;0' and fit_models[i] != 'None':

				# If local fitting approach the data will be averaged and plotted with error bars. Else the full data set will be plotted.
				if fit_modes[i] == 'Local':
					y_mean = replicate_mean_error(unique_conc, xs, ys)[0]
					std_dev = replicate_mean_error(unique_conc, xs, ys)[1]
					plot_func(fig, graph_name, unique_conc, y_mean, plot_markers[i], x_titles[i], y_titles[i],
							  subplot_row[i], subplot_col[i], 'Scatter_error', i, param_dict, color_list)
					plot_func(plot_fig, graph_name, unique_conc, y_mean, plot_markers[i], x_titles[i], y_titles[i],
							  subplot_row[i], subplot_col[i], 'Scatter_error', i, param_dict, color_list)
				else:
					plot_func(fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i],
							  subplot_col[i], 'None', i, param_dict, color_list)
					plot_func(plot_fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i],
							  subplot_col[i], 'None', i, param_dict, color_list)

				# Extract boundaries for the fit from the excel sheet and turn to list of floats. If no interval is supplied code will automaticallt assign 0;0.
				fitting_interval = fit_intervals[i].split(';')
				for j in range(0, len(fitting_interval)):
					fitting_interval[j] = float(
						fitting_interval[j].replace(",", "."))  # Replace commas with proper dots to get proper numbers.

				# Running fitting function and extracting parameters.
				if fit_models[i] == 'Hill':
					scatter_fit(hill_equation, 'Hill', xs, ys, fit_modes[i], fitting_interval[0], fitting_interval[1],
								i, ID_list, sample_notes, fit_models, fig, graph_name, x_titles, y_titles, subplot_row,
								subplot_col, param_dict, color_list)

				elif fit_models[i].lower() == 'hill_simple':
					scatter_fit(hill_simple, 'Hill_simple', xs, ys, fit_modes[i], fitting_interval[0],
								fitting_interval[1], i, ID_list, sample_notes, fit_models, fig, graph_name, x_titles,
								y_titles, subplot_row, subplot_col, param_dict, color_list)

				elif fit_models[i].lower() == '4pl':
					scatter_fit(fit_4pl, '4PL', xs, ys, fit_modes[i], fitting_interval[0], fitting_interval[1], i,
								ID_list, sample_notes, fit_models, fig, graph_name, x_titles, y_titles, subplot_row,
								subplot_col, param_dict, color_list)

				# Update color value using color_selector function
				color_count_global = color_selector(color_count_global,	color_list)

			else:
				master_dict['ID_list_new'].append(ID_list[i])
				master_dict['notes_list'].append(sample_notes[i])
				master_dict['model_list'].append(' ')
				master_dict['fit_parameters'].append(' ')
				master_dict['R_square'].append(' ')
				master_dict['KD_fit'].append(' ')

				plot_func(fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i],
						  subplot_col[i], 'None', i, param_dict, color_list)
				plot_func(plot_fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i],
						  subplot_col[i], 'None', i, param_dict, color_list)
				color_count_global = color_selector(color_count_global,
													color_list)  # Update color value using color_selector function

		table_plot(fig, ['Samples', 'Sample notes', 'Fitting models', 'Fitted KD/EC50', 'R^2', 'Fit parameters'],
				   [master_dict['ID_list_new'], master_dict['notes_list'], master_dict['model_list'],
					master_dict['KD_fit'],
					master_dict['R_square'], master_dict['fit_parameters']], ID_list)

		# Add dropdown
		fig.update_layout(
			updatemenus=[
				dict(
					type="buttons",
					direction="left",
					buttons=list([dict(args=["type", "scatter"], label="Graphs", method="restyle"),
								  dict(args=["type", "table"], label="Stats", method="restyle")]),
					pad={"r": 10, "t": 10}, showactive=True, x=0.11, xanchor="left", y=1.1, yanchor="top")])

	elif args.plot_type in ['boxplot', 'Boxplot']:

		for i in range(len(ID_list)):

			print('Analysing data: ' + str(ID_list[i]))

			x_id = 'x' + str(i + 1)
			y_id = 'y' + str(i + 1)

			graph_name = ID_list[i]

			# Extracting the raw x and y values from the excel sheet.
			# Replace commas with dots in case user has e.g. Danish Excel.
			xs = df[x_id][pd.to_numeric(df[x_id], errors='coerce').notnull()]
			ys = df[y_id][pd.to_numeric(df[y_id], errors='coerce').notnull()]
			xs = pd.to_numeric(xs.astype(str).str.replace(",", "."))
			ys = pd.to_numeric(ys.astype(str).str.replace(",", "."))

			# Extract values from python_misc list
			python_misc_list = python_misc_global[i].split(';')

			if python_misc_list.count('Box_all') > 0:
				plot_func(fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'Boxplot_all', i, param_dict, color_list)
			else:
				plot_func(fig, graph_name, xs, ys, plot_markers[i], x_titles[i], y_titles[i], subplot_row[i], subplot_col[i], 'Boxplot', i, param_dict, color_list)

	os.makedirs(args.output, exist_ok=True)
	pd_out.to_csv(args.output.joinpath(args.input_file.stem + '.tsv'), sep='\t')

	# Write output file
	output_file_name = args.output.joinpath(args.input_file.stem + '.html')
	if os.path.isfile(output_file_name):
		print('\n')
		new_file = input('A file with that name already exist. Do you wish to overwrite (Y/N)? ')
		if new_file.lower() == 'y' or new_file.lower() == 'yes':
			fig.write_html(output_file_name)
		else:
			print('I did not make any output file')
	else:
		fig.write_html(output_file_name)

	if args.plotting:
		plot_fig.write_image(args.output.joinpath(args.input_file.stem + '.svg'), width=1000, height=650)


if __name__ == "__main__":
	args = parse_args()
	if args.advanced_option_box:
		pdict = advanced_box()
	else:
		pdict = set_param_dict()
	main_func(args, pdict)
