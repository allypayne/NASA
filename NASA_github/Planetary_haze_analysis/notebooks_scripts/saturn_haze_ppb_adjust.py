import os
import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

# -------------------------------------------------------------------------------------------------------------

#                              READING IN DATA FOR ALL FUNCTIONS/ SYTNAX

# -------------------------------------------------------------------------------------------------------------

folder='/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/'
# direction to the default jupiter configuration file from PSG
# the saturn compatible file also includes the type of output data that we will want (i/f albedo) and spectral range 0.2-3 micron (among other variations)
saturn_base_cfg= open(folder+"psg_cfg_SATURN_best_haze.txt", "r");

# reads/ stores each line of file in basecfg (as a list)
basecfg= saturn_base_cfg.readlines();
#close the file everytime
saturn_base_cfg.close();

# Read in Karkoschka data
karkoschka_fromgv= np.genfromtxt('../jupiter_spectra/KarkoschkaAlbedoJupiter.txt')


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

print("Welcome to the Hazes Function!")

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

# compare to CBF data also (will compare to 2 of the different phase angles of saturn)
coulter_barnes_fortney_saturn_fig6= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/saturn_spectra/dbf6.txt')
coulter_barnes_fortney_saturn_fig7= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/saturn_spectra/dbf7.txt')

# since the data is organized by phase angle, separate the data into lists for the plots below
cbf_fig6= np.genfromtxt('./saturn_spectra/dbf6.txt') # read in figure 6 from cbf paper

phase_40_saturn=[]
phase_40_wavelength_saturn=[]
phase_40_albedo_saturn=[]
phase_60_saturn=[]
phase_60_wavelength_saturn=[]
phase_60_albedo_saturn=[]
phase_90_saturn=[]
phase_90_wavelength_saturn=[]
phase_90_albedo_saturn=[]
phase_110_saturn=[]
phase_110_wavelength_saturn=[]
phase_110_albedo_saturn=[]
for line in cbf_fig6:
    if line[0]<40:
        phase_40_saturn.append(line)
        phase_40_wavelength_saturn.append(line[1])
        phase_40_albedo_saturn.append(line[3])
    if line[0]>40 and line[0]<60:
        phase_60_saturn.append(line)
        phase_60_wavelength_saturn.append(line[1])
        phase_60_albedo_saturn.append(line[3])
    if line[0]>90 and line[0]<100:
        phase_90_saturn.append(line)
        phase_90_wavelength_saturn.append(line[1])
        phase_90_albedo_saturn.append(line[3])
    if line[0]>110:
        phase_110_saturn.append(line)
        phase_110_wavelength_saturn.append(line[1])
        phase_110_albedo_saturn.append(line[3])

albedo_cloudy_fig7=[]
wavelength_cloudy_fig7=[]
for line in coulter_barnes_fortney_saturn_fig7:
    if line[0]== 'a':
        albedo_cloudy_fig7.append(line[1])
        wavelength_cloudy_fig7.append(line[3])


# Since we don't have data for phase angle=0 use the phase integral (PSG eqn 77) to try to simulate by dividing the values for 40 by the coefficient for 40 degree phase
def get_scale_factor(phase_angle): # phase angle in degrees
    scale_val= (1/np.pi)*(np.sin(np.deg2rad(phase_angle))+(np.pi-np.deg2rad(phase_angle))*np.cos(np.deg2rad(phase_angle)))
    return scale_val 
    
# Molecules included: 'C2H2', 'C2H4', 'C2H6', 'CH4', 'H2O', 'NH3', 'WaterIce', 'NH4SH', 'Ammonia'

# -------------------------------------------------------------------------------------------------------------

#            TO RUN A RANGE OF FUNCTIONS ON PSG AND FIND THE PARAMETERS THAT IMPROVE THE FIT THE MOST

# -------------------------------------------------------------------------------------------------------------

diff_vals=[]
def vary_particle_params(abund_min, abund_max, min_size, max_size, abundance_interval, particle_size_interval):   #min_level, max_value, inc_abun (the amount we are increasing the abundance by)
    # this is the range of abundance scale factors we want to consider (units of ppm)
    for value_1 in np.arange(abund_min,abund_max,abundance_interval):
        # this is the range of particle sizes we want to consider (in units of micron)
        for value_2 in np.arange(min_size,max_size,particle_size_interval):

            # Make the config file:
            file_name= '/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/new_cfg_files/Saturn_HAZE_PPBabundance_'+str(value_1)+"_micronsize_"+str(value_2)+'.txt'
            f = open(file_name, 'w')
            # copy the original jupiter file to vary the values before saving
            new_file= basecfg.copy() 

            # Define the locations of the lines that I want to change
            loc_aabun=0
            for line in basecfg:
                if "<ATMOSPHERE-AABUN" in line:  
                    break
                loc_aabun= loc_aabun+1
            #print(loc_aabun)

            loc_asize=0
            for line in basecfg:
                if "<ATMOSPHERE-ASIZE" in line:  
                    break
                loc_asize= loc_asize+1
            #print(loc_asize)

            for line in basecfg:
                if "<ATMOSPHERE-AABUN" in line:
                    line_new_1= '<ATMOSPHERE-AABUN>1,1,'+str(value_1)+'\n'
                    #print(line_new_1)

                if "<ATMOSPHERE-ASIZE" in line:
                    line_new_2= "<ATMOSPHERE-ASIZE>10,0.01,"+str(value_2)+'\n'
                    #print(line_new_2)


            new_file[loc_aabun]= line_new_1
            new_file[loc_asize]= line_new_2

            for line in new_file:
                f.write(str(line))
            f.close() 

            # RUN THROUGH PSG
            print("Sending amazing science to PSG...")
            print("  ")

            run_file= 'curl --data-urlencode file@'+file_name+' https://psg.gsfc.nasa.gov/api.php > ./output_psg_spectra/'+'Saturn_HAZE_PPBabundance_'+str(value_1)+"_micronsize_"+str(value_2)+'OUT_.txt'
            os.system(run_file)
            print("PSG Status: Complete")
            
        #  ~~~~~~~~~~~~~~ Plotting + Optimize
            # Read in data files
            output_file_NEW_psg=np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/vary_params/output_psg_spectra/'+'Saturn_HAZE_PPBabundance_'+str(value_1)+"_micronsize_"+str(value_2)+'OUT_.txt')

            # the spectrum that we want to match is Karkoschka
            saturn_ideal_x= karkoschka_fromgv[:,0]/1000 # to make it in microns (the psg rad file is already in microns)
            saturn_ideal_y= karkoschka_fromgv[:,4] # saturn index

            # use the n and l max adjusted parameters to compare to (close to the true data)
            psg_default_data_saturn= np.genfromtxt('/Users/apayne3/Desktop/Project_Part2_Colors/psg_rad_saturn_n1_L80.txt')
            x_vals_default= psg_default_data_saturn[:,0] # wavelength in microns
            y_vals_default= psg_default_data_saturn[:,1]

            # Cubic Spline functions
            f_karkoschka_initial_spectrum = CubicSpline(x_vals_default,y_vals_default, bc_type='natural')
            f_karkoschka_new_data_spectrum = CubicSpline(output_file_NEW_psg[:,0],output_file_NEW_psg[:,1], bc_type='natural')

            new_psg_vs_kark= np.sqrt(np.sum((saturn_ideal_y-f_karkoschka_new_data_spectrum(saturn_ideal_x))**2)/(len(saturn_ideal_y)-2))
            old_psg_vs_kark= np.sqrt(np.sum((saturn_ideal_y-f_karkoschka_initial_spectrum(saturn_ideal_x))**2)/(len(saturn_ideal_y)-2))
            delta_rmse= old_psg_vs_kark-new_psg_vs_kark

            diff_vals.append(delta_rmse)

            if new_psg_vs_kark<old_psg_vs_kark:
                print("This change improved the fit!")
                print("  ")
                title='BETTER'
                # COULD ADD LINE HERE TO MAKE THE BASE CFG FOR THE NEXT RUN=OUTPUT OF THIS ONE (IF IMPROVED)
            else:
                print("The new fit is worse than the default data.")
                title='WORSE'
            print("New data vs ideal:",new_psg_vs_kark)
            print("Default data vs ideal:", old_psg_vs_kark)

            ###
            # ~~~~~~~~~~~~~~~~~ The small/ zoomed plot is here ~~~~~~~~~~~~~~~~
            ###

            # plt.plot((karkoschka_fromgv[:,0])/1000, karkoschka_fromgv[:,3], label='Karkoschka', linewidth=1, color='royalblue')
            plt.figure(figsize=(7,4), dpi=150)
            plt.plot((saturn_ideal_x), saturn_ideal_y, label='Karkoschka (IDEAL spectrum)', linewidth=0.7, color='dodgerblue')
            plt.plot(phase_40_wavelength_saturn,phase_40_albedo_saturn, label='CBF Saturn Phase angle= 40 degrees', linewidth=0.7, color='royalblue')
            plt.plot(output_file_NEW_psg[:,0], output_file_NEW_psg[:,1], label='PSG NEW Model; RMSE:{}'.format(new_psg_vs_kark), color= 'mediumseagreen')
            plt.plot(psg_default_data_saturn[:,0], psg_default_data_saturn[:,1], label='PSG Default; RMSE:{}'.format(old_psg_vs_kark), linewidth=0.7, linestyle='--', color='violet')
            plt.grid()
            plt.ylim(0,1)
            plt.xlim(0.2,1)
            plt.title('Zoomed on Karkoschka Data; {} Fit; PSG Haze Abundance {}, Particle Size: {}'.format(title, value_1, value_2))
            plt.xlabel(r'Wavelength [$\mu$m]')
            plt.ylabel('I/F')
            plt.xticks(np.arange(0,4.8, 0.2))
            plt.grid(True, linewidth=0.4, linestyle='--') 
            plt.legend(fontsize=7)
            plt.xlim(0.3,1)
            plt.show()
            
            ###
            # The larger plot is here 
            ###
            plt.figure(figsize=(9,5), dpi=200)
            plt.plot((saturn_ideal_x), saturn_ideal_y, label='Karkoschka (IDEAL spectrum)', linewidth=0.7, color='dodgerblue')
            plt.plot(phase_40_wavelength_saturn,phase_40_albedo_saturn, label='CBF Saturn Phase angle= 40 degrees', linewidth=0.7, color='royalblue')
            plt.plot(phase_40_wavelength_saturn,phase_40_albedo_saturn/get_scale_factor(40), label='Simulated Phase Angle=0', linewidth=0.7, color='teal', linestyle='--')
            plt.plot(output_file_NEW_psg[:,0], output_file_NEW_psg[:,1], label='PSG NEW Model; RMSE:{}'.format(new_psg_vs_kark), color= 'mediumseagreen')
            plt.plot(psg_default_data_saturn[:,0], psg_default_data_saturn[:,1], label='PSG Default; RMSE:{}'.format(old_psg_vs_kark), linewidth=0.7, linestyle='--', color='violet')
            plt.plot(wavelength_cloudy_fig7,albedo_cloudy_fig7, label='CBF Dayside, cloudy Saturn', linewidth=0.7, color='blueviolet', linestyle='--')
            plt.grid()
            plt.ylim(0,1)
            plt.xlim(0,4.8)
            plt.title('Saturn Spectra; {}'.format(title))
            plt.xlabel(r'Wavelength [$\mu$m]')
            plt.ylabel('I/F')
            plt.xticks(np.arange(0,4.8, 0.5))
            plt.grid(True, linewidth=0.4, linestyle='--') 
            plt.legend(fontsize=7)
            # plt.text(2, 0.33, 'Features here are more impacted by larger particle sizes', fontsize=7, color='blue')
            # plt.text(0.1, 0.57, 'This side is more impacted by', fontsize=7, color='blue')
            # plt.text(0.1, 0.55, 'smaller particles', fontsize=7, color='blue')
            plt.show();

    return diff_vals


#vary_particle_params(abund_min=1, abund_max=4, min_size=0.01, max_size=.11)