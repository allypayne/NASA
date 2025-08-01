import scipy as sc
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt


# READING IN CO2_ICE ABSORPTION DATA FROM PSG/CRISM
co2_crism_data= np.genfromtxt('crism_co2_ice.txt')
                             
# choose the SPECTRAL RANGE FOR THE CO2 ICE FEATURE/ portion of the data we are interested in fitting
um_min = 1.82 #microns
um_max = 2.2  #microns

# FUNCTION TO FIND THE DATA POINTS/ INDEX VALUES IN THE ARRAY OF REFLECTANCE DATA THAT ARE CLOSEST TO THE SPECIFIED UM_MIN/UM_MAX VALUES
def find_nearest_l(array, value):
    n = (np.abs(np.array(array)-value)).argmin()
    return n

# TRUE DATA FOR MARS ALBEDO/ REFLECTANCE VALUES
real_mars= pd.read_csv('organized_mars_albedo_data.csv')
del real_mars['Unnamed: 0']

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#######################__________DEFINING THE DATA POINTS TO BE CONSIDERED IN ANALYSIS __________###################
# giving the offsets an initial value to observe our plots/ data
a, b = 0.3 , .5


# REAL_MARS DATA
# --------------
# list of WAVELENGTH/x values (defining the x and y terms here makes things slightly easier moving forward)
real_mars_xvalues= real_mars['Wavelength']
real_mars_yvalues= real_mars['0']

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#######################_______ FIND INDEX VALUES FOR DATA__________###################


# REAL_MARS DATA
# --------------

# implement the find_nearest function to get the index values that define the boundary of our lower "resolution" data set (less data points over this range)
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

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#######################__________CUBIC SPLINE FUNCTION__________###################
# This will create a best fit to the data and interpolate it to ensure that each line has the same # of points/ spatial resolution amongst the data to do the spectral analysis

# CRISM CO2 ABSORPTION DATA
# --------------------------------

# This fits a function for the co2 crism data (this will allow us to interpolate the data for the points that we need (the location of the real_mars x values))
# this allows us to reduce the "resolution" of the CO2 model (the points that we fit with fco2) to match the lower resolution of the real_mars data set 
fCO2 = CubicSpline(co2_abs_x,co2_abs_y, bc_type='natural')

# unsure how to use optimize.minimize (maybe this will be useful later, but I didnt really try and created my own version of this)
#     print(sc.optimize.minimize(fit_co2_ice_abundance ,[a_list[value],b_list[value]], args = real_mars[xmin_r:xmax_r]['0']))


# sets up a local cubicspline function for the co2 data (fco2)
# The return statement of this function doesn't really work for whatever reason... (we move forward with just needing the fco2 portion to interpolate the data)
def fit_co2_ice_abundance(params,real_mars_xvals):
    fco2 = CubicSpline(co2_abs_x,co2_abs_y, bc_type='natural')
    params= a,b
    # switched to real mars y vals
    return np.sum(  np.abs( real_mars_yvalues[xmin_realdata:xmax_realdata]  -(a*fco2(real_mars_xvalues[xmin_realdata:xmax_realdata])) -  b)  )
    
# These lists of a (multiplies the resultant y value created by fco2 by this variable) and b (this is the vertical offset in the y value from fco2 data as shown below) 
# These lists were choosen to be close to the values that I knew the function should converge to (through some trail and error plot tests) and then choosing a large # of points in this range to improve the accuracy of the values that were chosen 
a_list= np.linspace(0.2,0.5,50)
b_list= np.linspace(-0.2,0.2,50)

### This function collects the offset values between each of the plots and it is stored in the "new_func" term
# the val/ diff_vals term stores the value associated with the original fit_co2_ice_abundance (this value seems to have an error as every answer appears exactly the same even though it should vary with the offsets)

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