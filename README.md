# SeoaneYurdagul2019
Seoane, Hernan; Yurdagul, Emircan (2019), Supporting Materials for "Trend shocks and Sudden Stops" by Seoane and Yurdagul (2019)




This file contains the codes and the data used to produce the results in "Trend Shocks and Sudden Stops" (August, 2018) by Hernan Seoane and Emircan Yurdagul.
First, we go over the files included for the model solution, and then through the data source and the figures in the paper.

README: SOLUTION

The 'solution_policies' file is the code we run to solve the policy functions. 
Load variable should be set to 0 to run without any guesses. Though the 'params' file should be copied in an 'outputs/' subfolder.
The 'simulations' file reads the policy function and gives the simulated sample.
To run the planner's problem, we just have to set the 'planner' variable to 1 in the 'solution_policies' file.
The codes are initially written to allow for shocks to kappa and omega. The current version still includes these in a degenerate way.
The random variables used for simulations are read from the files titled "randomXXX" included in this folder. They should also be copied in an 'outputs/' subfolder.

README: DATA & MATLAB
A file with the smoothed estimates of unobservables is included.

't_nt_crosscountry' includes data for figure 1, this data comes directly from the WDI and has multiple countries.

File 'ss_data_arg.mat' includes the average path for observables and smoothed estimates around sudden stops. This data is based on the raw data from Ferreres (2010) as cited in the paper.

Matlab files: 'figure_1', 'figures_2_3' plot the figures in the paper. For 'figures_2_3' you need to specify first the path to the simulations of the model with only transitory shocks and with trend and transitory shocks.
