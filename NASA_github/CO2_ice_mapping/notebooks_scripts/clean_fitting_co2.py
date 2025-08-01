import scipy as sc
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt


# READING IN CO2_ICE ABSORPTION DATA FROM PSG/CRISM
co2_crism_data= np.genfromtxt('crism_co2_ice.txt')
                             
# SPECTRAL RANGE FOR THE CO2 ICE FEATURE
um_min = 1.82 #microns
um_max = 2.2  #microns

# FUNCTION FIND THE DATA POINTS/ INDEX VALUES IN THE ARRAY OF REFLECTANCE DATA THAT ARE CLOSEST TO THE SPECIFIED UM_MIN/UM_MAX VALUES
def find_nearest_l(array, value):
    n = (np.abs(np.array(array)-value)).argmin()
    return n

# TRUE DATA FOR MARS ALBEDO/ REFLECTANCE VALUES
real_mars= pd.read_csv('organized_mars_albedo_data.csv')
del real_mars['Unnamed: 0']

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#######################__________DEFINING THE DATA POINTS TO BE CONSIDERED IN ANALYSIS __________###################
a, b = 0.3 , .5


# REAL_MARS DATA
# --------------

# list of WAVELENGTH/x values 
real_mars_xvalues= real_mars['Wavelength']
real_mars_yvalues= real_mars['0']

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#######################_______ FIND INDEX VALUES FOR DATA__________###################


# REAL_MARS DATA
# --------------
xmin_realdata,xmax_realdata = find_nearest_l(real_mars_xvalues,um_min),find_nearest_l(real_mars_xvalues,um_max)

# print(xmax_r,xmin_r) # or #print(wavelength_truedata[xmin_r:xmax_r]) # can check index values here if I want to

# CRISM CO2 ABSORPTION DATA VALUES
# --------------------------------
# y_values= co2_crism_data[:,1] 
xmin_co2crism,xmax_co2crism = find_nearest_l(co2_crism_data[:,0],um_min),find_nearest_l(co2_crism_data[:,0],um_max)


# CRISM CO2 ABSORPTION DATA VALUES
# --------------------------------
co2_abs_x= co2_crism_data[xmin_co2crism:xmax_co2crism,0]
co2_abs_y= co2_crism_data[xmin_co2crism:xmax_co2crism,1]


# 😃😃😃😃😃😃😃😃




#######################__________CUBIC SPLINE FUNCTIONS__________###################
# This will create a best fit to the data and interpolate it to ensure that each line has the same # of points/ spatial resolution amongst the data to do the spectral analysis

# REAL_MARS DATA
# --------------

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# CRISM CO2 ABSORPTION DATA
# --------------------------------

# This fits a function for the co2 crism data (this will allow us to interpolate the data for the points that we need (the location of the real_mars x values))
fCO2 = CubicSpline(co2_abs_x,co2_abs_y, bc_type='natural')


#     print(sc.optimize.minimize(fit_co2_ice_abundance ,[a_list[value],b_list[value]], args = real_mars[xmin_r:xmax_r]['0']))


# sets up a local cubicspline function for the co2 data

def fit_co2_ice_abundance(params,real_mars_xvals):
    fco2 = CubicSpline(co2_abs_x,co2_abs_y, bc_type='natural')
    params= a,b
    # switched to real mars y vals
    return np.sum(  np.abs( real_mars_yvalues[xmin_realdata:xmax_realdata]  -a*fco2(real_mars_xvalues[xmin_realdata:xmax_realdata]) -  b)  )



# def fit_co2_ice_abundance_2(params,real_mars_xvals):
#     fco2 = CubicSpline(co2_abs_x,co2_abs_y, bc_type='natural')
#     params= a,b
#     return fco2(real_mars_xvalues[xmin_realdata:xmax_realdata])
    #return np.sum(  np.abs( co2_abs_y  -a*fco2(real_mars_xvalues[xmin_realdata:xmax_realdata]) -  b)  )
    


a_list= np.linspace(0,0.5,10)
b_list= np.linspace(-1,2,20)

#print(np.sum(np.abs(real_mars_xvalues  -a*fCO2(real_mars_xvalues[xmin_realdata:xmax_realdata]) -  b) ))

#a,b =  sc.optimize.minimize(fit_co2_ice_abundance,[0.01,0.01],args=(real_mars_xvalues[xmin_realdata:xmax_realdata])  )
#print(a,b)`
   
#for value in range(len(a_list)):
 #   print(np.sum(  np.abs( real_mars_xvalues  -a_list[value]*fCO2(real_mars_xvalues[xmin_realdata:xmax_realdata]) -  b_list[value]) ))


# FUNCTION 1
# diff_vals=[]
# for value in range(len(a_list)):
#     val= fit_co2_ice_abundance_2([a_list[value],b_list[value]], real_mars_xvalues[xmin_realdata:xmax_realdata])
#     diff_vals.append(val)
#     print("a value:", a_list[value])
#     print(val)
# print("the standard deviation across the varying params is:", np.std(diff_vals))

# FUNCTION 2
diff_vals=[]
a_vals=[]
b_vals=[]
graph_diff=[]
for value in range(len(a_list)):
    val= fit_co2_ice_abundance([a_list[value],b_list[value]], real_mars_xvalues[xmin_realdata:xmax_realdata])
    # the real difference in values below
    new_func= np.sum(  np.abs( real_mars_yvalues[xmin_realdata:xmax_realdata]  -a_list[value]*fCO2(real_mars_xvalues[xmin_realdata:xmax_realdata])-b_list[value])  )
    graph_diff.
    diff_vals.append(val)
    a_vals.append(a_list[value])
    b_vals.append(b_list[value])
    print("a value:", a_list[value])
    print("b value:", b_list[value])
    print(val)
    #plt.plot(real_mars_xvalues[xmin_realdata:xmax_realdata],real_mars_yvalues[xmin_realdata:xmax_realdata],label='real mars')
    plt.plot(real_mars_xvalues[xmin_realdata:xmax_realdata],a_list[value]*fCO2(real_mars_xvalues[xmin_realdata:xmax_realdata])-b_list[value],label=f'CO2 MODEL, a:{a_list[value]}')
    plt.plot(co2_abs_x,-a_list[value]*co2_abs_y-b_list[value], label='CO2 CRISM')
    if b_list[value] <0:
        plt.plot(co2_abs_x,-a_list[value]*co2_abs_y-b_list[value], label='CO2 CRISM', color='green', linewidth=3)
plt.grid()
plt.plot(real_mars_xvalues[xmin_realdata:xmax_realdata],real_mars_yvalues[xmin_realdata:xmax_realdata],label='real mars', linewidth=4, color='red')
plt.legend()
plt.show()
print("the standard deviation across the varying params is:", np.std(diff_vals))
print(diff_vals)

values_frame = pd.DataFrame({'difference': diff_vals, 'a_values': a_vals, 'b_values': b_vals})
print(values_frame)

'''
plt.plot(real_mars_xvalues[xmin_realdata:xmax_realdata],real_mars_yvalues[xmin_realdata:xmax_realdata],label='real mars')
plt.plot(real_mars_xvalues[xmin_realdata:xmax_realdata],0.45*fCO2(real_mars_xvalues[xmin_realdata:xmax_realdata])-0.1,label='CO2 MODEL')
plt.plot(co2_abs_x,0.4*co2_abs_y-0.08, label='CO2 CRISM')
plt.legend()
plt.grid()
plt.show()
'''

'''
    
print(fit_co2_ice_abundance([150,-1], real_mars_xvalues[xmin_realdata:xmax_realdata]))

#'

# 😐😐😐😐😐😐😐😐

a,b = .38,0.05

plt.title('CO2 Model Data vs True Data')
plt.plot(co2_crism_data[xmin_co2crism:xmax_co2crism,0],a*co2_crism_data[xmin_co2crism:xmax_co2crism,1]-b, label= 'CRISM CO2')
plt.plot(real_mars_xvalues[xmin_realdata:xmax_realdata],real_mars[xmin_realdata:xmax_realdata]['0'], label='Mars True Reflectance')
plt.legend()
plt.grid()
plt.show()


#'
fCO2 = CubicSpline(co2_abs_x,co2_abs_y, bc_type='natural')
print(fCO2)




# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#######################__________ __________###################


# fitting the crism data to the real data
# fitting the lower resolution fit from cubic spline to the true data
print(fit_co2_ice_abundance([a,b],real_mars[xmin_r:xmax_r]['0']))



# a,b = .3,0
# plt.plot(co2_crism_data[xmin_co2crism:xmax_co2crism,0],a*co2_crism_data[xmin_co2crism:xmax_co2crism,1]-b, label= 'crism')
# plt.plot(wavelength_truedata[xmin_r:xmax_r],real_mars[xmin_r:xmax_r]['0'], label='Mars Reflectance')
# plt.legend()
# plt.show()



print(xmin_co2crism, xmax_co2crism)



####


# qr = np.arange(0,len(df[xmin:xmax,1]))

a, b = 0.3 , .5


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
plt.plot(wavelength_truedata[xmin_r:xmax_r], a*fCO2(wavelength_truedata[xmin_r:xmax_r]) -  b)
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

    

print(sc.optimize.minimize(fit_co2_ice_abundance ,[a,b], args = real_mars[xmin_r:xmax_r]['0']))



# THIS IS THE TRUE DATA SPECTRUM




def correct_baseline(params,yvalues):
    fHT = CubicSpline(wavelength_truedata,np.exp(-3*real_mars['0']), bc_type='natural')
    a, b  = params
    return np.sum(  np.abs( fHT(real_mars[xmin_r:xmax_r,0]) - yvalues  -  b * range(len(real_mars[xmin:xmax,0])))  )


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

