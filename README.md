# Plotly_data_analysis

A one-stop solution analysis scripts for Andreas' PhD. Script is continuously developed as I learn how to do more
experiments :-)

The script is developed for routine analysis of data that are deposited in Excel sheets. As input the scripts takes an
Excel file that contains certain columns specifying how the analysis should be performed. The data to be analyzed within
the file is specified by putting x1, x2, x3 etc. and y1, y2, y3 etc. above the data columns.

Running the script:  
python $SCRIPT.py -s $ANALYSIS_TYPE (AKTA, FIDA, Bar, Scatter, Panta, Octet, Boxplot) -i $EXCEL_FILE.xlsx -advanced (optional)

Columns for specifying how the analysis should be done:

- "IDs": list of data samples to be analyzed. First entry in the column corresponds to x1/y1 and so on. If different
  data samples are given the same name these will be grouped together in the plot legend.
- "Notes": notes will be included in the analysis table.
- "Axis_titles": two strings separated by ";" denote title of x axis and y axis.
- "Plot_markers": specify the type of plotting marker to use. Can be "Line" or "Dots".
- "Sub_plot": allows user to create subplots in the plot. Takes two integers separated by ";". Denotes the row number
  and column number, of the subplot to put the data into.
- "Data_interval": takes to numbers separated by ";". Allows user to slice data on x value to only include part of the
  data in the analysis without manually changing the data in the excel file.
  
AKTA MODULE:
- "AKTA_fraction": Specific to AKTA analysis. Two numbers separated by ";". Specifies which part of the data is considered to be the peak i.e. which part of the data to use for yield calculations. Multiple peaks can be specified using "|". 
- "AKTA_baseline": Specific to AKTA analysis. Two numbers separated by ";". Specifies the outer bounds for creating a
  linear baseline. The baseline can also be set in a different automatic manner under advanced analysis option.
- "AKTA_extinc_coeff": Specific to AKTA analysis. Takes a protein-specific extinction coefficient for yield calculation.
  The extinction coefficient should be given in Abs 0.1% i.e. the M-1 cm-1 extinction coeff. divided by MW.
- "AKTA_volume_load": Specific to AKTA analysis. Takes the amount of supernatant loaded during chromatography. Used for
  calculating the culture yield.
- "Fit_model": Specify the models used for peak fitting. Separate with semicolon (;). Can be: "gaussian", "skewed_gaussian", "exp_gaussian", "lorentzian". 
- "Fit_misc": 
    * "center_guess": used as initial values for guessing peak centers. Specify a value for each peak separated by semicolon (;). 
  
FIDA MODULE: 
- "Fitting_interval": Used for fitting data. Two numbers separated by ";". Denotes the x value interval that should be
  passed to the fitting function.
- "Fit_model": Used for fitting data. Can be "1to1 or "Excess" for FIDA fitting. Can be "Hill", "Hill_simple" or "4PL"
  for general scatter fitting.
- "Fit_approach": Used for fitting data. Can be "Local" or "Global". Local fitting will calculate mean y values and
  standard errors for identical x values and then do fitting on the mean values. "Global fitting" will treat each data
  set individually.
- It is possible to create an "ignore$ID" column next to the data that can be used for excluding certain data points from the fitting. 
  The $ID refers to the ID of the data. 
  
SCATTER MODULE: 
- "Fitting_interval": Used for fitting data. Two numbers separated by ";". Denotes the x value interval that should be
  passed to the fitting function.
- "Fit_model": Used for fitting data. Can be "1to1 or "Excess" for FIDA fitting. Can be "Hill", "Hill_simple" or "4PL"
  for general scatter fitting.
- "Fit_approach": Used for fitting data. Can be "Local" or "Global". Local fitting will calculate mean y values and
  standard errors for identical x values and then do fitting on the mean values. "Global fitting" will treat each data
  set individually.
  
PANTA MODULE:

OCTET MODULE:
- Currently only support global fitting. 


- "Python_misc": Enables user to customize the analysis. Takes several arguments separated by ";".
    1. "logx" forces x axis logarithmic.
    2. "logy" forces y axis logarithmic.
    3. "secondary_y" specifies the data to be plotted on a secondary y axis.
    4. "baseline_x=$VALUE". Allows user to specify an x value that will be used as a zero point i.e. data will be normalized so y is zero in this point. 
    5. "normalize_to_max". Will normalize the data to the max value of the data so the data is really in % of max. 
    6. "table_pos=$VALUE". Allows user to specify positions of data samples in the analysis table. Can make comparison
       easier in case of large datasets.
    5. "Box_all" will include all the datapoints next to the boxplot in the figure.  
    6. "ignore_data" can be used for excluding data from the analysis while allowing it to stay in the excel sheet. 
