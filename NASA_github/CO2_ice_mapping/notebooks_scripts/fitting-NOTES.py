import scipy as sc
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt



# CREATING THIS ANNOTATED COPY FOR REFERENCE
# THIS IS CO2 ABSORPTION DATA?

# f =  open('./CO2-CO2_2018.cia', "r"); CIAf = f.readlines(); f.close()
# cias280 = np.genfromtxt(CIAf[9917:10718])

# READING IN CO2_ICE DATA FROM PSG/CRISM
co2_crism_data= np.genfromtxt('crism_co2_ice.txt')
                             
# THIS IS THE SPECTRAL RANGE, FOR THE CO2 ICE FEATURE I AM ADJUSTING THIS AND THE UNITS ARE IN MICRONS NOW
um_min = 1.82 #microns
um_max = 2.2 #microns


# THIS WILL FIND THE DATA POINTS/ INDEX VALUES IN THE ARRAY OF REFLECTANCE DATA THAT ARE CLOSEST TO THE SPECIFIED UM_MIN/UM_MAX VALUES FOR THE CO2 SPECTRAL FEATURE
def find_nearest_l(array, value):
    n = (np.abs(np.array(array)-value)).argmin()
    return n

# READ IN TRUE DATA FOR MARS ALBEDO/ REFLECTANCE VALUES
real_mars= pd.read_csv('organized_mars_albedo_data.csv')
del real_mars['Unnamed: 0']

#define the list of wavelength values (doing this through pandas for convenience of the csv file format)
wavelength_truedata= real_mars['Wavelength']
print(wavelength_truedata)


# USE THE NEAREST FUNCTION TO FIND THE RANGE OF DATA WE WANT TO PULL FROM THE REAL_MARS DATA
# INITIALLY SETTING THIS UP FOR JUST REGION 1 --> WILL CREATE FOR LOOP NEXT

# THE CODE THAT WAS HERE BEFORE:
#h2obg = np.genfromtxt('./psg_trn_h2o_293.txt')
# xmax_r,xmin_r = find_nearest_l(h2obg[:,0],cm_min),find_nearest_l(h2obg[:,0],cm_max)
# fh2o = CubicSpline(h2obg[xmin_r:xmax_r,0][::-1], 1-h2obg[xmin_r:xmax_r,2][::-1], bc_type='natural')

# ADJUSTING IT TO INCLUDE THE ARRAYS THAT I AM WORKING WITH NOW
xmin_r,xmax_r = find_nearest_l(wavelength_truedata,um_min),find_nearest_l(wavelength_truedata,um_max)
# XMIN_R AND XMAX_R CORRESPOND TO THE -INDEX VALUES- IN THE LIST "WAVELENGTH_TRUEDATA" THAT ARE WITHIN THE SPECIFIED RANGE
print(xmax_r,xmin_r)
#print(wavelength_truedata[xmin_r:xmax_r])


# Cupic Spline Function between CO2 and true data

# KEEP THIS CONSISTENT ACROSS CO2 ABS DATA
# fco2 = CubicSpline(wavelength_truedata[xmin_r:xmax_r], 1-(real_mars['0'][xmin_r:xmax_r]), bc_type='natural')
#print(fco2)

# Dont think I need this portion (only considering 1 abs feature)
                   
# testfile = './20231102_CO2_temp52_pre230_res1_scan64_apt3_0.txt'
# df = np.genfromtxt(testfile,delimiter=',')[::-1]

# xmin,xmax = find_nearest_l(df[:,0],cm_min),find_nearest_l(df[:,0],cm_max)


                   
# USE FIND NEAREST TO FIND THE INDEX VALUES IN "co2_crism_data" THAT CORRESPOND TO THE WAVELENGTH RANGE
# xmin,xmax == 1.9,2.1, this find nearest code finds the index of data point that best corresponds to the wavelength that I want
                   
# y_values= co2_crism_data[:,1] 
xmin_co2crism,xmax_co2crism = find_nearest_l(co2_crism_data[:,0],um_min),find_nearest_l(co2_crism_data[:,0],um_max)




# a,b = .3,0
# plt.plot(co2_crism_data[xmin_co2crism:xmax_co2crism,0],a*co2_crism_data[xmin_co2crism:xmax_co2crism,1]-b, label= 'crism')
# plt.plot(wavelength_truedata[xmin_r:xmax_r],real_mars[xmin_r:xmax_r]['0'], label='Mars Reflectance')
# plt.legend()
# plt.show()



print(xmin_co2crism, xmax_co2crism)



####


# qr = np.arange(0,len(df[xmin:xmax,1]))

a, b = 0.3 , .5

co2x= co2_crism_data[xmin_co2crism:xmax_co2crism,0]
co2y= co2_crism_data[xmin_co2crism:xmax_co2crism,1]
fCO2 = CubicSpline(co2x,co2y, bc_type='natural')


a,b = .3,.1
plt.plot(wavelength_truedata[xmin_r:xmax_r],fCO2(wavelength_truedata[xmin_r:xmax_r]), label= 'crism_CO2')
plt.plot(wavelength_truedata[xmin_r:xmax_r],real_mars[xmin_r:xmax_r]['0'], label='Mars Reflectance')
plt.legend()
#plt.show()


# FCO2 HERE
# fco2 = CubicSpline(co2_crism_data[xmin_co2crism:xmax_co2crism,0], co2_crism_data[xmin_co2crism:xmax_co2crism,1], bc_type='natural')

# params= 0.3, 5
# co2x= co2_crism_data[xmin_co2crism:xmax_co2crism,0]
# co2y= co2_crism_data[xmin_co2crism:xmax_co2crism,1]


def fit_co2_ice_abundance(params,real_mars_obs):
    fco2 = CubicSpline(co2x,co2y, bc_type='natural')
    params= a,b
    return np.sum(  np.abs( real_mars_obs  -a*fco2(wavelength_truedata[xmin_r:xmax_r]) -  b)  )

#print(fit_co2_ice_abundance(params,co2y))


# fitting the crism data to the real data
# fitting the lower resolution fit from cubic spline to the true data
print(fit_co2_ice_abundance([a,b],real_mars[xmin_r:xmax_r]['0']))

plt.plot(wavelength_truedata[xmin_r:xmax_r], real_mars[xmin_r:xmax_r]['0'], label=' original data')
plt.plot(wavelength_truedata[xmin_r:xmax_r],        a*fCO2(wavelength_truedata[xmin_r:xmax_r]) -  b)
plt.legend()
#plt.show()


# first do function to minimize
a_list= np.linspace(0,2,20)
b_list= np.linspace(0,2,20)

# for value in range(len(a_list)):
#     print(sc.optimize.minimize(fit_co2_ice_abundance ,[a_list[value],b_list[value]], args = real_mars[xmin_r:xmax_r]['0']))
 
    
# The output of this code shows that the values that represent the difference of the 2 functions are not varying as the a and b parameters are adjusted (as they should)
# Next step is to figure out why they don't vary
# Resolve this and then use the optimize function to pick the best parameters
print("values here")
for value in range(len(a_list)):
    print(fit_co2_ice_abundance([a_list[value],b_list[:-value]],real_mars[xmin_r:xmax_r]['Wavelength']))
    print("a value:", a_list[value])
    print("b value:", b_list[value])

    
'''
print(sc.optimize.minimize(fit_co2_ice_abundance ,[a,b], args = real_mars[xmin_r:xmax_r]['0']))



# THIS IS THE TRUE DATA SPECTRUM




def correct_baseline(params,yvalues):
    fHT = CubicSpline(wavelength_truedata,np.exp(-3*real_mars['0']), bc_type='natural')
    a, b  = params
    return np.sum(  np.abs( fHT(real_mars[xmin_r:xmax_r,0]) - yvalues  -  b * range(len(real_mars[xmin:xmax,0])))  )
'''

'''
yvalues= real_mars['0']
# VINCENT CODE
def correct_baseline(params,yvalues):
    fco2 = CubicSpline(co2x,co2y, bc_type='natural')
    a, b  = params
    return np.sum(  np.abs( fco2(real_mars[xmin_r:xmax_r,0]) - real_mars[xmin_r:xmax_r,1]  -  b * range(len(real_mars[xmin_r:xmax_r,0])))  )

#sc.optimize.minimize(fit_co2_ice_abundance,[0.01,0.01],args=(real_mars[xmin_r:xmax_r]['0']))


## Works up to here 
a,b =  sc.optimize.minimize(fit_co2_ice_abundance,[0.01, 0.01], args=(real_mars[xmin_r:xmax_r]['0']))
# correct water is the function to minimize, which is done by scaling the values a b and c to get as close to 0 as possible
# args=(df[xmin:xmax,1],N) -- the y values, N= density

print(a,b)
#plt.plot(df[xmin:xmax,0],correct_water_test([0.3462850314894456,1.0917760712327427,0.6576826280286024],df[xmin:xmax,1],N)+1)
plt.plot(real_mars[xmin_r:xmax_r,0],real_mars[xmin_r:xmax_r,1])

#plt.plot(df[xmin:xmax,0],fh2o(df[xmin:xmax,0])+1)
y_corrected = sc.signal.medfilt((real_mars[xmin_r:xmax_r,1]-b*fco2(co2_crism_data[xmin_co2crism:xmax_co2_crism,0])),7)
plt.plot(real_mars[xmin_r:xmax_r,0],y_corrected  )
'''