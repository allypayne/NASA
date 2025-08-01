import scipy as sc
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
 

coulter_barnes_fortney= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/coulterbarnesfortney_fig4_jup.txt')
    
phase_zero_wavelength=[]
phase_zero_albedo=[]

for row in coulter_barnes_fortney:
    if row[0]<2:
        #print(row[1])
        #row[1], row[3])
        phase_zero_wavelength.append(row[1])
        phase_zero_albedo.append(row[3]) 
jupiter_ideal_x= phase_zero_wavelength
jupiter_ideal_y= phase_zero_albedo

# CRISM CO2 ABSORPTION DATA VALUES
# --------------------------------
psg_NEW_data= np.genfromtxt('vary_params/output_psg_spectra/Jupiter_WaterIce_min0_max12_scale5_OUT.txt')
#output_file_psg='/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/output_psg_spectra/'+'Jupiter_'+str(absorber)+'_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'_OUT.txt'
#psg_NEW_data= np.genfromtxt(output_file_psg)
new_fit_psg_x= psg_NEW_data[:,0]
new_fit_psg_y= psg_NEW_data[:,1]

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#######################__________CUBIC SPLINE FUNCTION__________###################
# This will create a best fit to the data and interpolate it to ensure that each line has the same # of points/ spatial resolution amongst the data to do the spectral analysis


f_idealspectrum = CubicSpline(jupiter_ideal_x,jupiter_ideal_y, bc_type='natural')
# sets up a local cubicspline function for the co2 data (fco2)
# The return statement of this function doesn't really work for whatever reason... (we move forward with just needing the fco2 portion to interpolate the data)
def fit_new_function(jupiter_ideal_x, jupiter_ideal_y, new_fit_psg_x, new_fit_psg_y):
    # set up cubicspline for the function that we want to match
    #f_idealspectrum = CubicSpline(jupiter_ideal_x,jupiter_ideal_y, bc_type='natural')
    #f_new_psg = CubicSpline(new_fit_psg_x,new_fit_psg_y, bc_type='natural')

    psg_default_data_jupiter= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/psg_rad_5micron.txt')
    x_vals_default= psg_default_data_jupiter[:,0]
    y_vals_default= psg_default_data_jupiter[:,1]
    new_psg_vs_cbf= np.sqrt(np.sum((y_vals_default-f_idealspectrum(new_fit_psg_x))**2)/(len(x_vals_default)-2))
    old_psg_vs_cbf= np.sqrt(np.sum((y_vals_default-f_idealspectrum(x_vals_default))**2)/(len(x_vals_default)-2))
    print(new_psg_vs_cbf)
    print(old_psg_vs_cbf)
    plt.plot(jupiter_ideal_x,jupiter_ideal_y, label='CBF (Ideal) Spectrum')
    plt.plot(new_fit_psg_x,new_fit_psg_y, label='New fit')    
    plt.plot(psg_default_data_jupiter[:,0],psg_default_data_jupiter[:,1], label='PSG DEFAULT', linestyle='--', color='lightgreen')
    plt.grid()
    plt.legend()
    plt.show()
    
    # params= a,b
    # switched to real mars y vals
    return print("This is the RSME value between the NEW PSG data and the CBF spectrum for Jupiter.", new_psg_vs_cbf), print("This is the RSME value between the defalt (old) PSG data and the CBF spectrum for Jupiter.", old_psg_vs_cbf)
fit_new_function(jupiter_ideal_x, jupiter_ideal_y, new_fit_psg_x, new_fit_psg_y)