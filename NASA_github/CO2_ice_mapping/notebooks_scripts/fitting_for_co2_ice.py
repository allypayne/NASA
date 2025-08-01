import scipy as sc
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt


##f =  open('crism_, "r"); CIAf = f.readlines(); f.close()
##cias280 = np.genfromtxt(CIAf[9917:10718])

# Read in reflectance data from CRISM for CO2 ice
co2_reflectance_data= np.genfromtxt('crism_co2_ice.txt')


# set parameters that correspond to the width of the CO2 feature
um_min = 1.82 
um_max = 2.2


def find_nearest_l(array, value):
    n = (np.abs(np.array(array)-value)).argmin()
    return n


# Read in the true data for Mars albedo/ reflectance
##h2obg = np.genfromtxt('./psg_trn_h2o_293.txt')

real_mars= np.genfromtxt('organized_mars_albedo_data.csv')
xmax_r,xmin_r = find_nearest_l(real_mars[:,0],um_min),find_nearest_l(real_mars[:,1],um_max)

fmars = CubicSpline(real_mars[xmin_r:xmax_r,0][::-1], 1-real_mars[xmin_r:xmax_r,2][::-1], bc_type='natural')


testfile = 'organized_mars_albedo_data.csv'
df = np.genfromtxt(testfile,delimiter=',')[::-1]

xmin,xmax = find_nearest_l(df[:,0],um_min),find_nearest_l(df[:,0],um_max)




# co2 spectrum
# cut co2 spectrum to relevant range

#cubic spline of Co2 spectrum

# use co2 ice data from crism


# xmin,xmax == 1.9,2.1
xmin,xmax = find_nearest_l(df[:,0],um_min),find_nearest_l(df[:,0],um_max)
qr = np.arange(0,len(df[xmin:xmax,1]))

fCO2 = CubicSpline(co2x,co2y, bc_type='natural')
a, b = 0.3 , +.5

def fit_co2_ice_abundance(params,yvalues,N ):
    fCO2 = CubicSpline(co2x,co2y, bc_type='natural')
    a, b = params
    return np.sum(  np.abs( yvalues  -a*fCO2(df[xmin:xmax,0]) -  b)  )

def correct_baseline(params,yvalues,N ):
    fHT = CubicSpline(cias280[:,0],np.exp(-3*cias280[:,1]*N**2), bc_type='natural')
    a, b, c  = params
    return np.sum(  np.abs( fHT(df[xmin:xmax,0]) - yvalues  -  b * range(len(df[xmin:xmax,0])) - c)  )


#a,b,c =  sc.optimize.minimize(correct_water,[a=0.01,b=0.01,c=0.01],args=(df[xmin:xmax,1],N)  )['x']
a,b,c =  sc.optimize.minimize(correct_water,[a,b,c],args=(df[xmin:xmax,1],N)  )['x']
# correct water is the function to minimize, which is done by scaling the values a b and c to get as close to 0 as possible
# args=(df[xmin:xmax,1],N) -- the y values, N= density
print(a,b,c)
#plt.plot(df[xmin:xmax,0],correct_water_test([0.3462850314894456,1.0917760712327427,0.6576826280286024],df[xmin:xmax,1],N)+1)
plt.plot(df[xmin:xmax,0],df[xmin:xmax,1])

#plt.plot(df[xmin:xmax,0],fh2o(df[xmin:xmax,0])+1)
y_corrected = sc.signal.medfilt((df[xmin:xmax,1]-b*fh2o(df[xmin:xmax,0])),7)
plt.plot(df[xmin:xmax,0],y_corrected  )