# Import modules
import tkinter as tk
from tkinter import ttk
import numpy as np


def create_advanced_box():
	# Define an empty dict for storing run-specific data in
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
		param_dict['error_marker'] = str(var5.get())
		param_dict['Baseline graph width'] = float(e19.get())
		param_dict['marker_size'] = float(e20.get())
		param_dict['min_peak_prominence'] = float(e21.get())
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
		param_dict['scatter_residuals'] = str(var6.get())
		param_dict['alpha_level'] = float(e34.get())

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
	ttk.Label(master, text="Error markers [layout]").grid(row=6, column=6)
	ttk.Label(master, text='X axis tick font size').grid(row=4, column=6)
	ttk.Label(master, text='Y axis tick font size').grid(row=4, column=9)
	ttk.Label(master, text="Baseline deg [AKTA]").grid(row=5, column=0)
	ttk.Label(master, text="Replace negative values with zero? [AKTA]").grid(row=6, column=0)
	ttk.Label(master, text="Path length in cm [AKTA]").grid(row=7, column=0)
	ttk.Label(master, text="Vertex point minimum prominence [Panta]").grid(row=8, column=0)
	ttk.Label(master, text="Degree variation for onset (%) [Panta]").grid(row=9, column=0)
	ttk.Label(master, text="Plot intermediate data? [Panta]").grid(row=10, column=0)
	ttk.Label(master, text="Bmin min [scatter fit]").grid(row=11, column=0)
	ttk.Label(master, text="Bmin max [scatter fit]").grid(row=11, column=3)
	ttk.Label(master, text="Bmax min [scatter fit]").grid(row=12, column=0)
	ttk.Label(master, text="Bmax max [scatter fit]").grid(row=12, column=3)
	ttk.Label(master, text="KD min [scatter fit]").grid(row=13, column=0)
	ttk.Label(master, text="KD max [scatter fit]").grid(row=13, column=3)
	ttk.Label(master, text="k_coop min [scatter fit]").grid(row=14, column=0)
	ttk.Label(master, text="k_coop max [scatter fit]").grid(row=14, column=3)
	ttk.Label(master, text="Plot residuals [scatter fit/FIDA]").grid(row=15, column=0)
	ttk.Label(master, text="Alpha level [statistics]").grid(row=16, column=0)

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
	e18 = ttk.OptionMenu(master, var3, 'plotly', 'plotly', 'plotly_white', 'simple_white')

	e19 = ttk.Entry(master)
	e20 = ttk.Entry(master)
	e21 = ttk.Entry(master)
	e22 = ttk.Entry(master)

	var4 = tk.StringVar(master)
	var4.set('no')
	e23 = ttk.OptionMenu(master, var4, 'no', 'no', 'yes')

	e24 = ttk.Entry(master)
	e25 = ttk.Entry(master)
	e26 = ttk.Entry(master)
	e27 = ttk.Entry(master)
	e28 = ttk.Entry(master)
	e29 = ttk.Entry(master)
	e30 = ttk.Entry(master)
	e31 = ttk.Entry(master)

	var5 = tk.StringVar(master)
	var5.set('Error bands')
	e32 = ttk.OptionMenu(master, var5, 'Error bands', 'Error bands', 'Error bars')

	var6 = tk.StringVar(master)
	var6.set('no')
	e33 = ttk.OptionMenu(master, var6, 'no', 'no', 'yes')

	e34 = ttk.Entry(master)

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
	e21.insert(10, "0")
	e22.insert(10, "0.5")
	e24.insert(10, "0")
	e25.insert(10, np.inf)
	e26.insert(10, "0")
	e27.insert(10, np.inf)
	e28.insert(10, "0")
	e29.insert(10, np.inf)
	e30.insert(10, -np.inf)
	e31.insert(10, np.inf)
	e34.insert(10, 0.05)

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
	e32.grid(row=6, column=7)
	e33.grid(row=15, column=1)
	e34.grid(row=16, column=1)

	ttk.Button(master, text='Run', command=master.quit).grid(row=20, column=1, sticky=tk.W, pady=4)
	ttk.Button(master, text='Register', command=show_entry_fields).grid(row=20, column=0, sticky=tk.W, pady=4)

	master.mainloop()
	return param_dict


def default_param_dict():
	return {
		'RI min': 0,
		'RI max': np.inf,
		'RIA min': 0,
		'RIA max': np.inf,
		'KD min': 0,
		'KD max': np.inf,
		'CI min': 0,
		'CI max': np.inf,
		'AKTA baseline mode': 'linear',
		'AKTA baseline deg': 2,
		'AKTA_pathlength': 0.2,
		'AKTA neg val handle': 'no',
		'graph width': 2,
		'xaxis_title_font_size': 12,
		'yaxis_title_font_size': 12,
		'xaxis_ticks_font_size': 12,
		'yaxis_ticks_font_size': 12,
		'plot_template': 'plotly',
		'error_marker': 'Error bands',
		'Baseline graph width': 2,
		'marker_size': 10,
		'min_peak_prominence': 0,
		'peak_onset_var_deg': 0.5,
		'panta_intermediate_plot': 'no',
		'Bmin_min': 0,
		'Bmin_max': np.inf,
		'Bmax_min': 0,
		'Bmax_max': np.inf,
		'KD_fit_min': 0,
		'KD_fit_max': np.inf,
		'k_coop_min': -np.inf,
		'k_coop_max': np.inf,
		'scatter_residuals': 'no',
		'alpha_level': 0.05
	}
