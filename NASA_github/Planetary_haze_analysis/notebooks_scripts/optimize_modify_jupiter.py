import os
import numpy as np
import matplotlib.pyplot as plt
# Path to folder containing the 5 atmosphere files from GV
# Files E-0 to E-5


# -------------------------------------------------------------------------------------------------------------

#                       IF YOU JUST WANT THE INPUT/ CONFIG FILES RUN THIS FUNCTION

# -------------------------------------------------------------------------------------------------------------

folder='/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/'
# direction to the default jupiter configuration file from PSG
# the jupiter compatible file also includes the type of output data that we will want (i/f albedo) and spectral range 0.2-3 micron (among other variations)
jupiter_base_cfg= open(folder+"psg_cfg_5micron_compat_test.txt", "r");

# reads/ stores each line of file in basecfg (as a list)
basecfg= jupiter_base_cfg.readlines();
#close the file everytime
jupiter_base_cfg.close();

## (Here we find "loc_molecules"- line where the base_config lists the molecules that) are considered in jupiters atmosphere (use this index to print the molecules)
# BASE CONFIG: "<ATMOSPHERE-LAYERS-MOLECULES>H2O,CO2,O3..."
loc_molecules=0
for line in basecfg:
    if "ATMOSPHERE-LAYERS-MOLECULE" in line:  
        break
    loc_molecules= loc_molecules+1

molecules=[]
for line in basecfg:
    if '<ATMOSPHERE-LAYERS-MOLECUL' in line:
        molecule= line.split(',')
        molecule[0]= molecule[0][-4:]
        for value in molecule:
            molecules.append(value)

print("MOLECULES in this template file:")
for value in molecules:
    print(value)

# keep track of the # of layers here
for line in basecfg:
    if "<ATMOSPHERE-LAYERS>" in line:
        num_layers= line.split('>')[1]
        num_layers_int= int(num_layers)

loc_layer1=0
for line in basecfg:
    if "<ATMOSPHERE-LAYER-1" in line:  
        break
    loc_layer1= loc_layer1+1
    
coulter_barnes_fortney= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/coulterbarnesfortney_fig4_jup.txt')
    
phase_zero_wavelength=[]
phase_zero_albedo=[]
for row in coulter_barnes_fortney:
    if row[0]<2:
        #print(row[1])
        #row[1], row[3])
        phase_zero_wavelength.append(row[1])
        phase_zero_albedo.append(row[3])
    
    
# Molecules listed: 'C2H2', 'C2H4', 'C2H6', 'CH4', 'H2O', 'NH3', 'WaterIce', 'NH4SH', 'Ammonia'
         
# Create code to adjust the cloud height of a given absorber in the list
# enter absorber as a string
        
def change_cloud_height_getfile(absorber, scale_abun, min_level, max_level):   #min_level, max_value, inc_abun (the amount we are increasing the abundance by)
    # consider this a scale factor for the cloud abundance, it will be automatically set to 1 unless specified in function
    # this creates range(0,9), have to add 2 bc the first two values are pressue and temp
    for value in range(len(molecules)):
        molecule_index=0
        if molecules[value][0:2]== 'C2':
            if molecules[value][0:4]== absorber[0:4]:
                molecule_index= value
                break
        if molecules[value][0:3]== absorber[0:3]:
            molecule_index= value
            #stop the loop once this requirement is satisfied
            break

    # This is where the true data from the ~original~ config file is stored 
    store_data_molecules_original = np.zeros(shape=(num_layers_int, len(molecules)+2))
    loc_layers=0
    for line in basecfg:
        if "<ATMOSPHERE-LAYER-" in line:
            abun_1= line.split(',')
            # adjust the formatting from the original config file
            abun_1[0]= abun_1[0].split('>')[1]
            #store_data_molecules[:,]
            loc_layers= loc_layers+1
            # this line stores all of the data into an array that can be adjusted to change the height of the cloud deck
            store_data_molecules_original[loc_layers-1]= abun_1
    # Will take the maximum abundance value from all of the rows and then multiply it by a given scale factor to increase the abundance of 
            #the clouds of a given absorber 
            # need to add 2 to the indexing bc P&T are the first 2 values
    
    #rounded_data = np.round(store_data_molecules_original, decimals=5)
    #np.savetxt('original_molecules.csv', rounded_data, delimiter=',')

    max_value_multiply= np.max(store_data_molecules_original[:,molecule_index+2])

    # copy the data frame so this version can be modified and compared to the original one
    store_data_molecules_newvals= np.zeros(shape=(num_layers_int, len(molecules)+2))
    for value in np.arange(0,len(store_data_molecules_original),1):
        store_data_molecules_newvals[value]= store_data_molecules_original[value]
   

    # Make the config file:
    file_name= '/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/new_cfg_files/Jupiter_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'.txt'
    f = open(file_name, 'w')
    # copy the original jupiter file to vary the values before saving
    new_file= basecfg.copy() 
    
    for line_val in np.arange(min_level, max_level):
        # Will take the maximum abundance value from all of the rows and then multiply it by a given scale factor to increase the abundance of 
        # the clouds of a given absorber
        new_val= (store_data_molecules_original[line_val,(molecule_index+2)]+ scale_abun*max_value_multiply)
        store_data_molecules_newvals[line_val,(molecule_index+2)]= new_val

    # New vals has length 30 (This corresponds to the 30 changed lines; does not have data for every layer)
    # for line in store_data_molecules_newvals:
    #     print("NEWVALS LINES:", line)
    
    for value in np.arange(min_level, max_level,1): 
        new_file[value+loc_layer1-1]= "<ATMOSPHERE-LAYER-" + str(value) + ">"
        line_values = ','.join(map(str, store_data_molecules_newvals[value]))
        new_file[value+loc_layer1-1] += line_values + '\n'

    for line in new_file:
        f.write(str(line))
    f.close() 

    return print(r"The atmospheric abundance values for {} have been changed...".format(molecules[molecule_index])), print(r"The old values for {}:".format(molecules[molecule_index]), store_data_molecules_original[:,molecule_index+2]), print(r"The new values for {}:".format(molecules[molecule_index]), store_data_molecules_newvals[:,molecule_index+2]), print("The file is saved as", file_name)

#return print(r"The atmospheric abundance values for {} have been changed...".format(molecules[molecule_index])), print(r"The old values for {}:".format(molecules[molecule_index]), store_data_molecules_original[:,molecule_index+2]), print(r"The new values for {}:".format(molecules[molecule_index]), store_data_molecules_newvals[:,molecule_index+2]), print("The file is saved as", file_name)

#print("Running cloud function...")
#change_cloud_height_getfile(absorber='Ammonia', scale_abun=20, min_level=30, max_level=45)

# -------------------------------------------------------------------------------------------------------------

#                       TO RUN THE FUNCTION ON PSG AND GET A RETURN PLOT

# -------------------------------------------------------------------------------------------------------------

def change_cloud_height_psgspectrum(absorber, scale_abun, min_level, max_level):   #min_level, max_value, inc_abun (the amount we are increasing the abundance by)
    # consider this a scale factor for the cloud abundance, it will be automatically set to 1 unless specified in function
    # this creates range(0,9), have to add 2 bc the first two values are pressue and temp
    print("Running cloud function...")
    for value in range(len(molecules)):
        molecule_index=0
        if molecules[value][0:1]== 'C2':
            if molecules[value][0:4]== absorber[0:4]:
                molecule_index= value
                break    
        if molecules[value][0:2]== absorber[0:2]:
            molecule_index= value
            #stop the loop once this requirement is satisfied
            break

    # This is where the true data from the ~original~ config file is stored 
    store_data_molecules_original = np.zeros(shape=(num_layers_int, len(molecules)+2))
    loc_layers=0
    for line in basecfg:
        if "<ATMOSPHERE-LAYER-" in line:
            abun_1= line.split(',')
            # adjust the formatting from the original config file
            abun_1[0]= abun_1[0].split('>')[1]
            #store_data_molecules[:,]
            loc_layers= loc_layers+1
            # this line stores all of the data into an array that can be adjusted to change the height of the cloud deck
            store_data_molecules_original[loc_layers-1]= abun_1
    # Will take the maximum abundance value from all of the rows and then multiply it by a given scale factor to increase the abundance of 
            #the clouds of a given absorber 
            # need to add 2 to the indexing bc P&T are the first 2 values
    max_value_multiply= np.max(store_data_molecules_original[:,molecule_index+2])

    # copy the data frame so this version can be modified and compared to the original one
    store_data_molecules_newvals= np.zeros(shape=(num_layers_int, len(molecules)+2))
    for value in np.arange(0,len(store_data_molecules_original),1):
        store_data_molecules_newvals[value]= store_data_molecules_original[value]
   

    # Make the config file:
    file_name= '/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/new_cfg_files/Jupiter_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'.txt'
    f = open(file_name, 'w')
    # copy the original jupiter file to vary the values before saving
    new_file= basecfg.copy() 
    
    for line_val in np.arange(min_level, max_level):
        # Will take the maximum abundance value from all of the rows and then multiply it by a given scale factor to increase the abundance of 
        # the clouds of a given absorber
        new_val= (store_data_molecules_original[line_val,(molecule_index+2)]+ scale_abun*max_value_multiply)
        store_data_molecules_newvals[line_val,(molecule_index+2)]= new_val
    
    for value in np.arange(min_level, max_level,1): 
        new_file[value+loc_layer1-1]= "<ATMOSPHERE-LAYER-" + str(value) + ">"
        line_values = ','.join(map(str, store_data_molecules_newvals[value]))
        new_file[value+loc_layer1-1] += line_values + '\n'

    for line in new_file:
        f.write(str(line))
    f.close() 

    print(r"The old values for {}:".format(molecules[molecule_index]), store_data_molecules_original[:,molecule_index+2])
    print(r"The new values for {}:".format(molecules[molecule_index]), store_data_molecules_newvals[:,molecule_index+2])

    # RUN THROUGH PSG
    print("Sending amazing science to PSG...")
    print("  ")
    #curl --data-urlencode file@config.txt https://psg.gsfc.nasa.gov/api.php
    run_file= 'curl --data-urlencode file@'+file_name+' https://psg.gsfc.nasa.gov/api.php > ./vary_params/output_psg_spectra/'+'Jupiter_'+str(absorber)+'_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'_OUT.txt'
    os.system(run_file)
    print("Good news, it worked!")

    # Create the plot here:

    #Project_Part2_Colors/jupiter_spectra/psg_rad_wide_default.txt
    # Read in data files
    jupiter_published= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/madden_spectra/CatalogofSolarSystemObjects/Albedos/Jupiter_Lundock080507_Albedo.txt')
    psg_default_data_jupiter= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/psg_rad_wide_default.txt')
    output_file_psg='/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/output_psg_spectra/'+'Jupiter_'+str(absorber)+'_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'_OUT.txt'
    psg_NEW_data= np.genfromtxt(output_file_psg)

    #figsize=(7,3)
    plt.figure(dpi=250)
    plt.plot(jupiter_published[:,0],jupiter_published[:,1], label='Published vals', color='forestgreen')
    plt.plot(psg_default_data_jupiter[:,0],psg_default_data_jupiter[:,1], label='PSG DEFAULT', linestyle='--', color='lightgreen')
    plt.plot(psg_NEW_data[:,0],psg_NEW_data[:,1], label='Jupiter Adjusted Clouds', color= 'blue')
    plt.plot(phase_zero_wavelength,phase_zero_albedo, label='Disk Integrated Jupiter: CBF 2022', color= 'aquamarine')
    plt.grid()
    plt.legend()
    plt.ylabel('Albedo')
    plt.xlabel('Wavelength ($\mu$m)')
    plt.title('Adjusted Jupiter Params: {} clouds at layers {}-{}, scaled by {}x'.format(molecules[molecule_index],str(min_level),str(max_level),str(scale_abun)))
    plt.ylim(0,1)
    plt.xlim(0.5,3)
    plt.show()
    #plt.savefig('/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/output_plots/'+'Jupiter_min'+str(min_level)+'_max'+str(max_level)+'_scale'+str(scale_abun)+'_PLOT.png');
    
    return print(r"Jupiter's {} cloud position/ abundance values have been changed...".format(molecules[molecule_index])), print("The input configuration file is saved as", file_name), print("The output is saved to /Users/apayne3/Desktop/Project_Part2_Colors/vary_params/output_psg_spectra")

#change_cloud_height_psgspectrum(absorber='Ammonia', scale_abun=20, min_level=30, max_level=45)












#############-----------------------##################
############# OPTIMIZE THE PLOTS HERE ################

import scipy as sc
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

# AFTER DATA IS CREATED FROM PSG IT WILL BE IN THIS ARRAY: psg_NEW_data

# Define the x and y vals that we want it to be close to (from coulter_barnes_fortney (line 57))
jupiter_ideal_x= phase_zero_wavelength
jupiter_ideal_y= phase_zero_albedo

# # implement the find_nearest function to get the index values that define the boundary of our lower "resolution" data set (less data points over this range)
# xmin_realdata,xmax_realdata = find_nearest_l(real_mars_xvalues,um_min),find_nearest_l(real_mars_xvalues,um_max)

# CRISM CO2 ABSORPTION DATA VALUES
# --------------------------------
psg_NEW_data= np.genfromtxt('./output_psg_spectra/Jupiter_WaterIce_min0_max12_scale5_OUT.txt')
new_fit_psg_x= psg_NEW_data[:,0]
new_fit_psg_y= psg_NEW_data[:,1]
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#######################__________CUBIC SPLINE FUNCTION__________###################
# This will create a best fit to the data and interpolate it to ensure that each line has the same # of points/ spatial resolution amongst the data to do the spectral analysis


#f_idealspectrum = CubicSpline(jupiter_ideal_x,jupiter_ideal_y, bc_type='natural')
# sets up a local cubicspline function for the co2 data (fco2)
# The return statement of this function doesn't really work for whatever reason... (we move forward with just needing the fco2 portion to interpolate the data)
def fit_new_function(jupiter_ideal_x, jupiter_ideal_y, new_fit_psg_x, new_fit_psg_y):
    # set up cubicspline for the function that we want to match
    #f_idealspectrum = CubicSpline(jupiter_ideal_x,jupiter_ideal_y, bc_type='natural')
    #f_new_psg = CubicSpline(new_fit_psg_x,new_fit_psg_y, bc_type='natural')

    jupiter_ideal_x= phase_zero_wavelength
    jupiter_ideal_y= phase_zero_albedo
    psg_default_data_jupiter= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/jupiter_spectra/psg_rad_wide_default.txt')
    x_vals_default= psg_default_data_jupiter[:,0]
    y_vals_default= psg_default_data_jupiter[:,1]

    f_cbf_initial_spectrum = CubicSpline(x_vals_default,y_vals_default, bc_type='natural')
    f_cbf_new_data_spectrum = CubicSpline(new_fit_psg_x,new_fit_psg_y, bc_type='natural')
    print("length when attempting to interpolate to ideal spectrum length:", len(f_cbf_initial_spectrum(jupiter_ideal_x)))
    print("length of ideal jupiter x",len(jupiter_ideal_x))
    print("length of x_default:", len(x_vals_default))

    #plt.plot(jupiter_ideal_x, f_cbf_initial_spectrum(jupiter_ideal_x), label= 'test plot')
    new_psg_vs_cbf= np.sqrt(np.sum((jupiter_ideal_y-f_cbf_new_data_spectrum(jupiter_ideal_x))**2)/(len(x_vals_default)-2))
    old_psg_vs_cbf= np.sqrt(np.sum((jupiter_ideal_y-f_cbf_initial_spectrum(jupiter_ideal_x))**2)/(len(x_vals_default)-2))
    print("New data vs ideal:",new_psg_vs_cbf)
    print("Default data vs ideal:", old_psg_vs_cbf)
    plt.plot(jupiter_ideal_x,jupiter_ideal_y, label='CBF (Ideal) Spectrum')
    plt.plot(new_fit_psg_x,new_fit_psg_y, label='New fit')    
    plt.plot(psg_default_data_jupiter[:,0],psg_default_data_jupiter[:,1], label='PSG DEFAULT', linestyle='--', color='lightgreen')
    plt.grid()
    plt.legend()
    plt.show()
    
    # params= a,b
    # switched to real mars y vals
    #return print("This is the RSME value between the NEW PSG data and the CBF spectrum for Jupiter.", new_psg_vs_cbf), print("This is the RSME value between the defalt (old) PSG data and the CBF spectrum for Jupiter.", old_psg_vs_cbf)

fit_new_function(jupiter_ideal_x, jupiter_ideal_y, new_fit_psg_x, new_fit_psg_y)
### This function collects the offset values between each of the plots and it is stored in the "new_func" term
# the val/ diff_vals term stores the value associated with the original fit_co2_ice_abundance (this value seems to have an error as every answer appears exactly the same even though it should vary with the offsets)
'''
# FUNCTION 2
diff_vals=[]
a_vals=[]
b_vals=[]
graph_diff=[]
for value in range(len(a_list)):
    for value_2 in range(len(b_list)):
        val= fit_co2_ice_abundance([a_list[value],b_list[value_2]], real_mars_xvalues[xmin_realdata:xmax_realdata])
        # the real difference in values below (as seen on the plot between the model (including offsets) and the true data
        new_func= np.sum(  np.abs( real_mars_yvalues[xmin_realdata:xmax_realdata]  -(a_list[value]*fCO2(real_mars_xvalues[xmin_realdata:xmax_realdata])-b_list[value_2]))  )
        graph_diff.append(new_func)
        diff_vals.append(val)
        a_vals.append(a_list[value])
        b_vals.append(b_list[value_2])
        print("a value:", a_list[value])
        print("b value:", b_list[value_2])
        print(val)
        
        # If you want to see all of the plots (not just the best fit below) then you can uncomment these lines (these are mainly just useful for constraining the boundaries of a and b lists

        #plt.plot(real_mars_xvalues[xmin_realdata:xmax_realdata],a_list[value]*fCO2(real_mars_xvalues[xmin_realdata:xmax_realdata])-b_list[value_2],label=f'CO2 MODEL, a:{a_list[value]}')
        #plt.plot(co2_abs_x,-a_list[value]*co2_abs_y-b_list[value_2], label='CO2 CRISM')
        # The line below is unneccessary but I used to it to constrain the bounds that I choose for a_list and b_list
        #if b_list[value_2] <0:
            #plt.plot(co2_abs_x,-a_list[value]*co2_abs_y-b_list[value_2], label='CO2 CRISM', color='green', linewidth=3)
#plt.grid()
#plt.plot(real_mars_xvalues[xmin_realdata:xmax_realdata],real_mars_yvalues[xmin_realdata:xmax_realdata],label='real mars', linewidth=4, color='red')
#plt.legend()
#plt.show()


# Get the min value here:
values_frame = pd.DataFrame({'difference': diff_vals, 'a_values': a_vals, 'b_values': b_vals, 'graph_diff': graph_diff})

# We want the a and b values that allow are model to align the best with the original data set (has the smallest graph_diff, or the lowest spatial separation between points in the y direction of the plots)
# This will tell you the values in the terminal for easier reference

min_val= values_frame['graph_diff'].min()
print("MIN VALUES")
min_a=[]
min_b=[]
graph_min_diff=[]
for index,row in values_frame.iterrows():
    if row['graph_diff']== min_val:
        min_a= row['a_values']
        min_b= row['b_values'] 
        graph_min_diff= row['graph_diff']
        print(row)
print("MIN A:,", min_a)
print("MIN B:,", min_b)

# FINAL PLOT: This will show you the best fit function where you can plot the original data set, the original data that you are modeling, and the interpolated model over the specified x value range
plt.plot(real_mars_xvalues[xmin_realdata:xmax_realdata],min_a*fCO2(real_mars_xvalues[xmin_realdata:xmax_realdata])-min_b,label=f'CO2 MODEL, a:{min_a}, b:{min_b}')
plt.plot(co2_abs_x,min_a*co2_abs_y-min_b, label='CO2 CRISM')
plt.grid()
plt.plot(real_mars_xvalues[xmin_realdata:xmax_realdata],real_mars_yvalues[xmin_realdata:xmax_realdata],label='real mars', linewidth=4, color='red')
plt.legend()
plt.show()
'''