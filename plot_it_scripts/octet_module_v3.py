
# Import modules
import tkinter as tk
from tkinter import ttk

from scipy.signal import savgol_filter
import lmfit
from scipy.optimize import curve_fit
from decimal import Decimal
from sklearn.metrics import r2_score

from plot_it_scripts.plotting_script import *
from plot_it_scripts.main_functions import *

def octet_main(sample_idx, user_input_dict, plot_dict, param_dict, master_dict, df):

    # Loading local variables for plotting
    figure = plot_dict['figure']
    plot_figure = plot_dict['plot_figure']
    graph_name = plot_dict['graph_names'][sample_idx]
    plot_marker = user_input_dict['plot_markers'][sample_idx]
    x_title = user_input_dict['x_titles'][sample_idx]
    y_title = user_input_dict['y_titles'][sample_idx]
    color_list = plot_dict['color_list']
    color_count = plot_dict['color_count']

    fit_dict = {}
    fit_dict['y_anal_global'] = []
    fit_dict['y_fit_global'] = []
    fit_dict['Req'] = []

    fit_mode = user_input_dict['fit_modes'][sample_idx]

    # Define a conversion scheme based on the concentration unit.
    if param_dict['octet_conc_unit'] == 'nM':
        fit_dict['unit_conversion'] = 10**9
    elif param_dict['octet_conc_unit'] == 'uM':
        fit_dict['unit_conversion'] = 10**6

    # Set the data IDs and get the appropriate data from dataframe.
    x_id = 'x' + str(sample_idx + 1)
    y_id = 'y' + str(sample_idx + 1)
    z_id = 'z' + str(sample_idx + 1)

    # Getting the right data from the dataframe. The data contains 2nd baseline, association and dissociation
    data = collect_data(df, x_id, y_id, z_id)

    # Get the analyte concentrations for each sensor.
    sensor_conc_list = df['conc' + str(sample_idx + 1)].dropna().tolist()

    # Prepare the data by:
    # 1=normalizing to reference, 2=align to y axis, 3=smooth with Savitz-Golay.
    data, sensor_ids, sensor_conc_list = prepare_data(data, x_id, y_id, z_id, sample_idx, sensor_conc_list, param_dict)
    fit_dict['concentrations'] = sensor_conc_list

    [master_dict['octet_sensor_conc'].append(x) for x in sensor_conc_list]

    for sensor in sensor_ids:

        # Add the sample name to the master dict
        master_dict['ID_list_new'].append(graph_name)

        # Add the sensor to the master dict
        master_dict['octet_sensors'].append('Sensor ' + sensor.upper())

        # Define the t0 values for both the association and dissociation
        fit_dict['t0_ass'], fit_dict['t0_dis'] = get_t0(data[x_id + sensor], data[z_id])

        # Getting simple data for plotting. This should be the full association/dissociation sensorgrams
        x_plot, y_plot = simple_plot_data(data[x_id + sensor], data[y_id + sensor], data[z_id])

        # Get the part of the data that is intended for use in the fitting procedure
        x_anal, y_anal = trim_data(data[x_id + sensor], data[y_id + sensor], data[z_id], param_dict, fit_dict)
        fit_dict['y_anal_global'].append(y_anal)

        # Define the time cutoff between association and dissociation
        fit_dict['cutoff'], fit_dict['cutoff_base'] = define_cutoff(data[x_id + sensor], data[z_id])

        # Plot the raw sensorgrams for the assocition/dissociation
        plot_func(figure, graph_name + '_' + sensor.upper(), x_plot, y_plot, 'None', plot_marker, x_title, y_title, 1,
                  1, 'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

        # If specified by user we also plot the discrete points on the association/dissociation
        if param_dict['octet_show_points'] == 'Yes':
            plot_func(figure, graph_name + '_' + sensor.upper(), x_plot, y_plot, 'None', 'Dots', x_title, y_title,
                      1, 1, 'None', sample_idx, param_dict, color_list, color_count, user_input_dict, master_dict)

        # Plot the full completely raw sensorgrams that show all phases of experiments
        plot_func(figure, graph_name + '_' + sensor.upper(), df[x_id + sensor], df[y_id + sensor], 'None',
                  plot_marker, x_title, y_title, 4, 1, 'None', sample_idx, param_dict, color_list, color_count,
                  user_input_dict, master_dict)


    if fit_mode.lower() == 'global':
        fit_result = global_fitting_procedure(np.array(x_anal), np.array(fit_dict['y_anal_global']), fit_dict)
        #print(fit_result.params.pretty_print())

        for s, sensor in enumerate(sensor_ids):

            # Calculate the fit and the residuals from the fit
            my_fit, residuals = calculate_fit(np.array(x_anal), np.array(fit_dict['y_anal_global'])[s], fit_result, fit_dict, s)

            # Save the fitted values.
            fit_dict['y_fit_global'].append(np.array(my_fit))

            # Plot the fit
            plot_func(figure, graph_name + '_' + sensor.upper() + '_fit',
                      np.array(x_anal), my_fit, 'None', plot_marker, x_title, y_title, 1,
                      1, 'AKTA_baseline', sample_idx, param_dict, color_list, color_count,
                      user_input_dict, master_dict)

            # Plot the residuals
            plot_func(figure, graph_name + '_' + sensor.upper(), np.array(x_anal), residuals, 'None',
                      plot_marker, x_title, y_title, 3, 1, 'None', sample_idx, param_dict, color_list, color_count,
                      user_input_dict, master_dict)

            # Calculate the equilibrium
            Rmax = fit_result.params['Rmax_'+str(s+1)].value
            kd = fit_result.params['kd_' + str(s + 1)].value
            ka = fit_result.params['ka_' + str(s + 1)].value
            Req = association_func_full(np.inf, sensor_conc_list[s], Rmax, kd, ka)
            fit_dict['Req'].append(Req)

            # Add the fitted parameters to the master dict for plotting
            fit_params_to_master_dict(s, fit_result, fit_dict['unit_conversion'], master_dict)

            # Plot the steady state
            plot_func(figure, graph_name + '_' + sensor.upper(), np.array(sensor_conc_list[s]), np.array(Req), 'None',
                      'dots', x_title, y_title, 1, 4, 'None', sample_idx, param_dict, color_list, color_count,
                      user_input_dict, master_dict)

            # Add the appropriate color to the table coloring list
            table_color_list_manager(color_list[color_count], plot_dict['table_color_list'])

        # Calculate full R2 value and add to master dict using function
        full_R2(fit_dict['y_anal_global'] , fit_dict['y_fit_global'], master_dict)

        # Fit the steady state
        ss_par, ss_cov = curve_fit(steady_state_model, np.array(sensor_conc_list), np.array(fit_dict['Req']))
        ss_Req = ss_par[0]
        ss_KD = ss_par[1]

        # Add the fitted steady-state KD to the output table
        ss_KD_table = '%.2E' % Decimal(ss_KD)
        [master_dict['octet_KD_SS'].append(ss_KD_table) for x in range(len(sensor_ids))]

        # Plot the steady state
        ss_pred = steady_state_model(np.array(sensor_conc_list), ss_Req, ss_KD)
        plot_func(figure, graph_name + '(steady state)', np.array(sensor_conc_list), ss_pred, 'None',
                  plot_marker, x_title, y_title, 1, 4, 'None', sample_idx, param_dict, color_list, color_count,
                  user_input_dict, master_dict)

        ss_R2 = r2_score(np.array(fit_dict['Req']), ss_pred)
        ss_R2 = round(ss_R2, 3)
        [master_dict['octet_R2_SS'].append(ss_R2) for x in range(len(sensor_ids))]


def octet_box(param_dict):
    def show_entry_fields():
        print('Input for Octet analysis have been registered')
        param_dict['octet_conc_unit'] = str(var1.get())
        param_dict['octet_ref_sensors'] = str(e1.get())
        param_dict['octet_savgol_window'] = str(e2.get())
        param_dict['octet_savgol_pol'] = str(e3.get())
        param_dict['octet_savgol_pol'] = str(e3.get())
        param_dict['octet_Yaxis_align'] = str(e4.get())
        param_dict['octet_interstep_correct'] = str(var2.get())
        param_dict['octet_assoc_start'] = str(e5.get())
        param_dict['octet_assoc_end'] = str(e6.get())
        param_dict['octet_dissoc_start'] = str(e7.get())
        param_dict['octet_dissoc_end'] = str(e8.get())
        param_dict['octet_show_points'] = str(var3.get())


    master = tk.Tk()

    ttk.Label(master, text="Concentration unit").grid(row=0, column=0)
    var1 = tk.StringVar(master)
    var1.set('nM')
    e9 = ttk.OptionMenu(master, var1, 'nM', 'nM', 'uM', 'mM', 'pM')
    e9.grid(row=0, column=1)

    ttk.Label(master, text="Reference sensor (letter)").grid(row=1, column=0)
    e1 = ttk.Entry(master)
    e1.insert(10, "0")
    e1.grid(row=1, column=1)

    ttk.Label(master, text="Savgol window").grid(row=2, column=0)
    e2 = ttk.Entry(master)
    e2.insert(10, "10")
    e2.grid(row=2, column=1)

    ttk.Label(master, text="Savgol polynomium").grid(row=3, column=0)
    e3 = ttk.Entry(master)
    e3.insert(10, "3")
    e3.grid(row=3, column=1)

    ttk.Label(master, text="Align y axis (None, Baseline 2, Association or Dissociation).").grid(row=4, column=0)
    e4 = ttk.Entry(master)
    e4.insert(10, "None")
    e4.grid(row=4, column=1)

    ttk.Label(master, text="Interstep correction").grid(row=5, column=0)
    var2 = tk.StringVar(master)
    var2.set('None')
    e10 = ttk.OptionMenu(master, var2, 'None', 'None', 'To association', 'To dissociation', 'To baseline')
    e10.grid(row=5, column=1)

    ttk.Label(master, text="Association start").grid(row=6, column=0)
    e5 = ttk.Entry(master)
    e5.insert(10, "Start")
    e5.grid(row=6, column=1)

    ttk.Label(master, text="Association end").grid(row=6, column=3)
    e6 = ttk.Entry(master)
    e6.insert(10, "End")
    e6.grid(row=6, column=4)

    ttk.Label(master, text="Dissociation start").grid(row=7, column=0)
    e7 = ttk.Entry(master)
    e7.insert(10, "Start")
    e7.grid(row=7, column=1)

    ttk.Label(master, text="Dissociation end").grid(row=7, column=3)
    e8 = ttk.Entry(master)
    e8.insert(10, "End")
    e8.grid(row=7, column=4)

    ttk.Label(master, text="Show points?").grid(row=8, column=0)
    var3 = tk.StringVar(master)
    var3.set('No')
    e11 = ttk.OptionMenu(master, var3, 'No', 'No', 'Yes')
    e11.grid(row=8, column=1)

    ttk.Button(master, text='Run', command=master.quit).grid(row=20, column=1, sticky=tk.W, pady=4)
    ttk.Button(master, text='Register', command=show_entry_fields).grid(row=20, column=0, sticky=tk.W, pady=4)

    master.mainloop()
    return param_dict

def collect_data(df, x_id, y_id, z_id):

    # Set column names lower case
    df.columns = df.columns.str.lower()

    # Only include columns that start with the right x or y id.
    filter_col = [col for col in df if col.lower().startswith((x_id, y_id, z_id))]
    data_full = df[filter_col]

    # Slice data to contain baseline 2, association and dissociation
    data_interaction = data_full[data_full[z_id].str.lower().isin(['baseline 2','association', 'dissociation'])]

    return data_interaction

def prepare_data(data, x_id, y_id, z_id, sample_idx, conc_list, param_dict):

    # Generating a letter list for the sensors.
    sensor_ids = np.unique([x.lower().replace(x_id, '').replace(y_id, '').replace(z_id, '') for x in data.columns]).tolist()
    sensor_ids.remove("")

    for sensor in sensor_ids:

        # Subtract reference data
        data = subtract_ref(data, x_id, y_id, sample_idx, sensor, param_dict)

        # Align the data to y axis
        data = y_axis_align(data, x_id, y_id, z_id, sensor, sample_idx, param_dict)

        # Perform interstep correction
        data = interstep_correct(data, x_id, y_id, z_id, sensor, sample_idx, param_dict)

        # Smooth data
        data = data_smooth(data, x_id, y_id, sensor, param_dict)

    # Remove the reference sensor trace from the data to avoid fitting reference.
    data = remove_ref_sensor_trace(data, x_id, y_id, sample_idx, param_dict)

    # Generating a letter list for the sensors.
    sensor_ids = np.unique([x.lower().replace(x_id, '').replace(y_id, '').replace(z_id, '') for x in data.columns]).tolist()
    sensor_ids.remove("")

    sensor_ids_new = []
    conc_list_new = []
    for a, conc in enumerate(conc_list):
        if conc != 0:
            conc_list_new.append(conc)
            sensor_ids_new.append(sensor_ids[a])

    return data, sensor_ids_new, conc_list_new

def data_smooth(data, x_id, y_id, sensor_id, param_dict):

    # Define the Savgol filter parameters
    if param_dict['octet_savgol_window'] == 'Default':
        savgol_window = 5
    else:
        savgol_window = int(param_dict['octet_savgol_window'])

    if param_dict['octet_savgol_pol'] == 'Default':
        savgol_pol = 3
    else:
        savgol_pol = int(param_dict['octet_savgol_pol'])

    # Smooth the data
    new_col = savgol_filter(data[y_id + sensor_id], savgol_window, savgol_pol)
    data = data.drop([y_id + sensor_id], axis='columns')
    data[y_id + sensor_id] = new_col

    return data

def subtract_ref(data, x_id, y_id, sample_idx, sensor_id, param_dict):

    # Getting sensor-specific data into arrays
    xs = data[x_id + sensor_id].to_numpy()
    ys = data[y_id + sensor_id].to_numpy()
    #zs = data[z_id].to_numpy()

    if param_dict['octet_ref_sensors'] == '0':
        data = data

    else:

        # Get the letter of the reference sensor that is supplied by user.
        ref_sensor = param_dict['octet_ref_sensors'].lower()

        # y values for the reference sensor
        ys_ref = data[y_id + ref_sensor].to_numpy()

        # Make new y values by subtracting reference sensor from the sensor data.
        ys_new = ys - ys_ref

        # Drop the old y value column and insert the new
        data = data.drop([y_id + sensor_id], axis='columns')
        data[y_id + sensor_id] = ys_new

    return data

def y_axis_align(data, x_id, y_id, z_id, sensor_id, sample_idx, param_dict):

    # Getting sensor-specific data into arrays
    xs = data[x_id + sensor_id].to_numpy()
    ys = data[y_id + sensor_id].to_numpy()
    zs = data[z_id].to_numpy()

    if param_dict['octet_Yaxis_align'] == 'None':
        data = data.copy()

    else:

        align_id = param_dict['octet_Yaxis_align']

        # The align data is all the y values where z value equals the user input.
        # The align value is the 25 points corresponding to the last 5 seconds.
        align_data = ys[zs == align_id]
        align_value = align_data[-25:].mean()

        # Drop the old y value column and insert the new
        data = data.drop([y_id + sensor_id], axis='columns')
        data[y_id + sensor_id] = ys - align_value

    return data

def interstep_correct(data, x_id, y_id, z_id, sensor_id, sample_idx, param_dict):

    # Getting sensor-specific data into arrays
    xs = data[x_id + sensor_id].to_numpy()
    ys = data[y_id + sensor_id].to_numpy()
    zs = data[z_id].to_numpy()

    # Getting y values for the baseline and the association
    ys_baseline = ys[zs == 'Baseline 2']
    ys_association = ys[zs == 'Association']
    ys_dissociation = ys[zs == 'Dissociation']

    # If interstep correcting to assocation we are moving the baseline on y axis to align to association.
    # This should be the same as moving the association and dissociation step together.
    if param_dict['octet_interstep_correct'] == 'To association':

        # Getting the time point where the baseline before assocation ends.
        baseline_end = ys_baseline[-1]

        # Getting the time point where the assocation after the baseline starts
        association_start = ys_association[0]

        # Calculate interstep difference
        step_diff = association_start - baseline_end

        ys[zs == 'Association'] = ys[zs == 'Association'] - step_diff
        ys[zs == 'Dissociation'] = ys[zs == 'Dissociation'] - step_diff

    # Move assocation step so last point of assocation aligns with dissociation
    elif param_dict['octet_interstep_correct'] == 'To dissociation':

        # Getting the time point where the assocation ends
        association_end = ys_association[-1]

        # Getting time point where dissociation starts
        dissociation_start = ys_dissociation[0]

        step_diff = association_end - dissociation_start

        ys[zs == 'Association'] = ys[zs == 'Association'] - step_diff


    # If interstep correcting to the baseline we are moving the association step so the start of association aligns
    # to the end of the baseline
    elif param_dict['octet_interstep_correct'] == 'To baseline':

        # Getting the time point where the baseline before assocation ends.
        baseline_end = ys_baseline[-1]

        # Getting the time point where the assocation after the baseline starts
        association_start = ys_association[0]

        # Calculate interstep difference
        step_diff = association_start - baseline_end

        # Move the assocation data so the first association point equals the last baseline point
        ys[zs == 'Association'] = ys[zs == 'Association']-step_diff

        # Insert the modified values back into the dataframe
        data[y_id + sensor_id] = ys

    return data

def remove_ref_sensor_trace(data, x_id, y_id, sample_idx, param_dict):

    if param_dict['octet_ref_sensors'] == '0':
        data = data

    else:

        ref_sensor = param_dict['octet_ref_sensors'].lower()

        #ref_sensor_list = param_dict['octet_ref_sensors'].split(',')

        # The variable stores the letter corresponding to the ref sensor from user
        #ref_sensor = ref_sensor_list[sample_idx].lower()

        data = data.drop(columns=[x_id + ref_sensor , y_id + ref_sensor])

    return data

def simple_plot_data(x_data, y_data, z_data):
    # Convert to arrays
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    z_data = np.array(z_data)

    # Create arrays specific for association and dissociation
    x_baseline = x_data[z_data == 'Baseline 2']
    x_baseline = x_baseline[len(x_baseline)-100:]
    y_baseline = y_data[z_data == 'Baseline 2']
    y_baseline = y_baseline[len(y_baseline) - 100:]
    x_ass = x_data[z_data == 'Association']
    y_ass = y_data[z_data == 'Association']
    x_dis = x_data[z_data == 'Dissociation']
    y_dis = y_data[z_data == 'Dissociation']

    x_out = np.concatenate((x_baseline, x_ass, x_dis), axis=None)
    y_out = np.concatenate((y_baseline, y_ass, y_dis), axis=None)

    return x_out, y_out

def trim_data(x_data, y_data, z_data, param_dict, fit_dict):

    # Convert to arrays
    x_data = np.array(x_data)
    y_data = np.array(y_data)
    z_data = np.array(z_data)

    # Create arrays specific for association and dissociation
    x_ass = x_data[z_data=='Association']
    y_ass = y_data[z_data == 'Association']
    x_dis = x_data[z_data == 'Dissociation']
    y_dis = y_data[z_data == 'Dissociation']

    # Create arrays with x values that are normalized to start in zero.
    # This is because each phase in principle starts in zero.
    x_ass_base = x_ass - fit_dict['t0_ass']
    x_dis_base = x_dis - fit_dict['t0_dis']

    # Set starting point and end point for the association
    if param_dict['octet_assoc_start'] != 'Start':
        ass_start = float(param_dict['octet_assoc_start'])
    else:
        ass_start = x_ass_base[0]

    if param_dict['octet_assoc_end'] != 'End':
        ass_end = float(param_dict['octet_assoc_end'])
    else:
        ass_end = x_ass_base[-1]

    # Set starting point and end point for the dissociation
    if param_dict['octet_dissoc_start'] != 'Start':
        dis_start = float(param_dict['octet_dissoc_start'])
    else:
        dis_start = x_dis_base[0]

    if param_dict['octet_dissoc_end'] != 'End':
        dis_end = float(param_dict['octet_dissoc_end'])
    else:
        dis_end = x_dis_base[-1]

    # Create arrays with x and y data that falls within the specified intervals
    x_ass_out = x_ass[(x_ass_base >= ass_start) & (x_ass_base <= ass_end)]
    y_ass_out = y_ass[(x_ass_base >= ass_start) & (x_ass_base <= ass_end)]

    x_dis_out = x_dis[(x_dis_base >= dis_start) & (x_dis_base <= dis_end)]
    y_dis_out = y_dis[(x_dis_base >= dis_start) & (x_dis_base <= dis_end)]

    x_out = np.concatenate((x_ass_out, x_dis_out), axis=None)
    y_out = np.concatenate((y_ass_out, y_dis_out), axis=None)

    return x_out, y_out

def define_cutoff(x_data, z_data):

    # Convert to arrays
    x_data = np.array(x_data)
    z_data = np.array(z_data)

    # Create arrays specific for association and dissociation
    x_ass = x_data[z_data == 'Association']
    x_dis = x_data[z_data == 'Dissociation']

    cutoff = (x_ass[-1] + x_dis[0])/2

    cutoff_base = x_dis[0] - x_ass[0]

    return cutoff, cutoff_base

def global_fitting_procedure(t, data, fit_dict):

    # Defining the list of concentrations.
    conc_list = fit_dict['concentrations']

    # Defining the fitting parameters.
    fitting_params = lmfit.Parameters()


    fitting_params.add('cutoff', value=fit_dict['cutoff'], vary=False)

    # Defining the conversion factor depending on which unit the user has inputted.
    fitting_params.add('conversion_factor', value=fit_dict['unit_conversion'], vary=False)

    for iy, y in enumerate(data):
        conc = conc_list[iy]

        # Loading in parameters for the fit
        #fitting_params.add(f'analyte_{iy + 1}', value=conc*(10**(-9)), vary=False)
        fitting_params.add(f'analyte_{iy + 1}', value=conc, vary=False)
        fitting_params.add(f'Rmax_{iy + 1}', value=0.1, min=0)
        fitting_params.add(f'ka_{iy + 1}', value=0.1, min=0)
        fitting_params.add(f'kd_{iy + 1}', value=0.1, min=0)

        # Define the x0 as the time point between association and dissociation.
        fitting_params.add('x0', value=fit_dict['cutoff_base'], vary=False)

        fitting_params.add(f'y0_{iy + 1}', expr=f'Rmax_{iy + 1}' + '*(1/(1+(' + f'kd_{iy + 1}' + '/(' + f'ka_{iy + 1}' + '*' + f'analyte_{iy + 1}' + '))))*(1-exp(-(' + f'ka_{iy + 1}' + '*' + f'analyte_{iy + 1}' + '+' + f'kd_{iy + 1}' + ')*x0))')

    # Fixing ka and kd across across all samples since it is a global fit
    for iy in range(1, len(data)):
        fitting_params[f'kd_{iy + 1}'].expr = 'kd_1'
        fitting_params[f'ka_{iy + 1}'].expr = 'ka_1'

    # Calculate the KD based on kinetic parameters. We just use the first value since these should all be the same
    # in the global fit.
    fitting_params.add('KD_kin', expr='kd_1/ka_1')

    # Running the minimization.
    out = lmfit.minimize(global_resid_minimize, fitting_params, args=(t, data), method='least_squares')

    return out

def global_resid_minimize(params, t, data):

    ndata, _ = data.shape
    resid = 0.0 * data[:]

    cutoff = params['cutoff'].value

    for i in range(ndata):

        data_new = data[i,:]

        analyte = params['analyte_%i' % (i + 1)].value
        Rmax = params['Rmax_%i' % (i + 1)].value
        kd = params['kd_%i' % (i + 1)].value
        ka = params['ka_%i' % (i + 1)].value
        y0 = params['y0_%i' % (i + 1)].value

        # If association data is present we calculate
        t_ass = t[t <= cutoff].astype(float)
        if len(t_ass) > 0:
            t_ass_base = t_ass - t_ass[0]
            data_ass = data_new[t <= cutoff]
            resid_ass = data_ass - association_func_full(t_ass_base, analyte, Rmax, kd, ka)

        t_dis = t[t > cutoff].astype(float)
        if len(t_dis) > 0:
            t_dis_base = t_dis - t_dis[0]
            data_dis = data_new[t > cutoff]
            resid_dis = data_dis - dissociation_func_full(t_dis_base, y0, kd)

        if len(t_ass) > 0 and len(t_dis) > 0:
            resid[i,:] = np.concatenate((resid_ass, resid_dis), axis=None)
        elif len(t_ass) > 0 and len(t_dis) == 0:
            resid[i,:] = resid_ass
        elif len(t_ass) == 0 and len(t_dis) > 0:
            resid[i,:] = resid_dis

    return resid.flatten()

def association_func_full(t, analyte, Rmax, kd, ka):

    A = 1 + (kd/(ka*analyte))
    B = np.exp(-(ka*analyte+kd)*t)

    return Rmax*(1/A)*(1-B)

def dissociation_func_full(t, y0, kd):

    return y0*np.exp(-kd*t)

def get_t0(x_data, z_data):

    # Convert to arrays
    x_data = np.array(x_data)
    z_data = np.array(z_data)

    # Create arrays specific for association and dissociation
    x_ass = x_data[z_data == 'Association']
    x_dis = x_data[z_data == 'Dissociation']

    t0_ass = x_ass[0]
    t0_dis = x_dis[0]

    return t0_ass, t0_dis

def calculate_fit(xs, ys, fit_result, fit_dict, s):

    cutoff = fit_dict['cutoff']
    t0_ass = fit_dict['t0_ass']
    t0_dis = fit_dict['t0_dis']

    analyte = fit_result.params['analyte_'+str(s+1)].value
    Rmax = fit_result.params['Rmax_' + str(s + 1)].value
    kd = fit_result.params['kd_' + str(s + 1)].value
    ka = fit_result.params['ka_' + str(s + 1)].value
    y0 = fit_result.params['y0_' + str(s + 1)].value

    # If association data is present we calculate
    xs_ass = xs[xs <= cutoff].astype(float)
    if len(xs_ass) > 0:
        xs_ass_base = xs_ass - t0_ass
        ys_ass = ys[xs <= cutoff]
        ass_pred = association_func_full(xs_ass_base, analyte, Rmax, kd, ka)
        ass_resid = ys_ass - ass_pred

    xs_dis = xs[xs > cutoff].astype(float)
    if len(xs_dis) > 0:
        xs_dis_base = xs_dis - t0_dis
        ys_dis = ys[xs > cutoff]
        dis_pred = dissociation_func_full(xs_dis_base, y0, kd)
        dis_resid = ys_dis - dis_pred

    if len(xs_ass) > 0 and len(xs_dis) > 0:
        pred_out = np.concatenate((ass_pred, dis_pred), axis=None)
        resid_out = np.concatenate((ass_resid, dis_resid), axis=None)
    elif len(xs_ass) > 0 and len(xs_dis) == 0:
        pred_out = ass_pred
        resid_out = ass_resid
    elif len(xs_ass) == 0 and len(xs_dis) > 0:
        pred_out = dis_pred
        resid_out = dis_resid

    return pred_out, resid_out

def steady_state_model(conc, Req, KD):

    return (Req*conc)/(KD + conc)

def fit_params_to_master_dict(sensor_idx, fit_result, conversion_factor, master_dict):

    Rmax = fit_result.params['Rmax_' + str(sensor_idx + 1)].value
    Rmax_err = fit_result.params['Rmax_' + str(sensor_idx + 1)].stderr

    ka = fit_result.params['ka_' + str(sensor_idx + 1)].value
    ka = ka * conversion_factor
    ka_round = '%.2E' % Decimal(ka)

    ka_err = fit_result.params['ka_' + str(sensor_idx + 1)].stderr
    ka_err = ka_err * conversion_factor
    ka_err_round = '%.2E' % Decimal(ka_err)

    kd = fit_result.params['kd_' + str(sensor_idx + 1)].value
    kd_round = '%.2E' % Decimal(kd)

    kd_err = fit_result.params['kd_' + str(sensor_idx + 1)].stderr
    kd_err_round = '%.2E' % Decimal(kd_err)

    kinetic_KD = fit_result.params['KD_kin'].value
    kinetic_KD_round = '%.2E' % Decimal(kinetic_KD)

    kinetic_KD_err = fit_result.params['KD_kin'].stderr
    kinetic_KD_err_round = '%.2E' % Decimal(kinetic_KD_err)

    master_dict['octet_ka'].append(ka_round)
    master_dict['octet_ka_err'].append(ka_err_round)
    master_dict['octet_kd'].append(kd_round)
    master_dict['octet_kd_err'].append(kd_err_round)
    master_dict['octet_kinetic_KD'].append(kinetic_KD_round)
    master_dict['octet_kinetic_KD_err'].append(kinetic_KD_err_round)

def full_R2(y_true, y_fit, master_dict):

    # Make sure both y arrays are actually arrays
    y_true = np.array(y_true)
    y_fit = np.array(y_fit)

    # Get the length of the y values. We basically count how many fits have been done.
    y_len = len(y_true)

    # Flatten both arrays
    y_true_flat = y_true.flatten()
    y_fit_flat = y_fit.flatten()

    # Calculate the R2
    r2_full = r2_score(y_true_flat , y_fit_flat)
    r2_full = round(r2_full, 3)

    # Add the value to the master dict. We need to add the value to the master dict the same number of times
    # as a fit has been made in order for the output table to match.
    for a in range(y_len):
        master_dict['octet_R2_full'].append(r2_full)
