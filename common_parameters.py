# Parameters for Models I, II, III, and IV, see the following paper:
# L. Moccia, and G. Laporte, 'Improved models for technology choice in a transit corridor with fixed demand'

import numpy as np

# PARAMETERS

# Parameters -- User, all Models
l   = 10.0  # Average journey length                                    (km)
v   = 4.0   # Access speed                                              (km/h)


# Parameters -- Values of time, all Models
Pa  = 12.5  # Value of access time savings              ($/h)
Pw  = 15.0  # Value of waiting time savings             ($/h)
Pv  = 10.0  # Value of in-vehicle time savings          ($/h)
mu  = 0.33  # ratio of Ph (value of time at home) and Pw 

#Alternative scenario, low income country
licf = 0.5       # low income country factor   
# Pa  = licf * Pa  # Value of access time savings              ($/h)
# Pw  = licf * Pw  # Value of waiting time savings             ($/h)
# Pv  = licf * Pv  # Value of in-vehicle time savings          ($/h)


# Parameters  -- Cost of waiting, all Models
epsilon	= 0.5   # rate of the average waiting time to the headway
fbar 	= 5.0   # Threshold frequency, for f < fbar user wait at home     (veh/h)
tw  	= 4.0   # Safety threshold time passengers spend waiting at station (min)
                # before the expected arrival if f >= fbar
t01 = tw / 60.0 # convert tw in hours                                       (h)
t11 = mu        #
t02 = 0.0       #                                                           (h)
t12 = 1.0       #
      
# Parameters -- Transit system, all Models
L   				= 20.0  	# Route length        (km) 
alpha 				= 0.35		# Passengers in the most loaded section, fraction                                                                  
nu  				= 0.9   	# Safety factor on capacity
min_demand_studied 	= 3000.0	# The cost function will be studied for values larger than this (pax)
max_demand_studied  = 60000.0	# The cost function will be studied for values smaller than this (pax)
step_demand         = 500.0     # The step for the demand studied


# Parameters --  Infrastructure cost, all Models

discount_rate           = 0.07             
land_cost_per_hectare   = 9.0                # Cost of land in million $ per hectare

infracost_m = [0.0,  10.0,    20.0,    35.0]  #Infrastructure cost      (million $/km)
width_m     = [0.0,  10.0,    10.0,    15.0]  # Width required by the infrastructure (m)
infralife_m = [0.0,  50.0,    100.0,   100.0] # Infrastructure technical life in years
inframaint_m= [0.0,  295.0,   1078.0, 1851.0] # Infrastructure maintenance cost ($/h)

# Rolling stock cost parameters
res_value_rs = 0.05                         # Residual value of the rolling stock
vehcost_m = [0.40,   0.62,    3.49,   2.77] # Vehicle cost           (million $/veh)
vehlife_m = [20.0,   20.0,   35.0,    35.0] # Vehicle technical life  (years)

# Stop cost, Models II, III,  IV, and V
stopcost_m  = [0.0, 0.125,    0.25,   0.5]  # Stop cost (million $/stop)


# Parameters -- Costs and mode specific coefficients 
# Source: Table 2 and A.1.2 Tirachini et al 2010
# in pos 0 Bus, followed by BRT, LRT, and HR
# HR adapted considering 3 vehicles per TU
# All models:
c1t_m   = [42.0,  42.0,    73.0,    130.0]# Unit operator cost per TU-hour, derived from Tirachini et al 2010   ($/TU-h)
#c1t_m   = [42.0,  42.0,   42.0,    42.0]# Unit operator cost per TU-hour, scenario of high labor productivity in rail systems   ($/TU-h)
#c1t_m   = [42.0 * licf,  42.0 * licf,    73.0 * licf,    130.0 * licf]# Unit operator cost per TU-hour, scenario low income country reduction   ($/TU-h)

#Model I
S_m     = [20.0,  30.0,    35.0,     40.0]# Speed,  with the same value of Tirachini et al 2010  (km/h)
#Model II, III, IV
Vmax_m      =  [22.0,  36.0,    44.5,     55.0]#  V_{max}. these values are such the avg running speed is equal to the S values of Tirachini et al 2010, results do not change significatively                      (km/h)
#Vmax_m     = [22.0,  44.5,    44.5,     55.0]# Scenario with equal Vmax for BRT and LRT

fmax_m  = [200.0,150.0,   80.0,     40.0]# Maximal frequency, as in Tirachini et 2010  (veh/h)
#fmax_m = [200.0,   60.0,   60.0,     40.0]# Maximal frequency, scenario with  urban segregation externality   (veh/h)

# Model I
d_m     = [ 0.4,   0.8,     1.0,      1.2]# Distance between stations  (as in Tirachini et al. 2010)       (km)
#d_m    = [ 0.4,   0.8,     0.8,      0.8]# Distance between stations (these values are used to highlight relevance of optimal stop spacing)   (km)


#Model I and II
c2_m    = [1.13,  1.42,    1.83,     3.32]# Unit operator cost per TU-km        ($/TU-km)
beta_m  = [ 4.0,  0.33,    0.33,     0.11]# Boarding and alighting time per pax (s-TU/pax)

#Model III,  IV, and V
c2v_m   = [1.13,  1.42,    1.83,     1.11]# Unit operator cost per vehicle-km   ($/veh-km)
kv_m    = [64.0, 101.0,   190.0,    250.0]# Vehicle capacity                    (pax/veh)
betav_m = [ 4.0,  0.33,    0.33,     0.33]# Boarding and alighting time per pax (s-veh/pax)


num_modes 			= len(fmax_m)  		  # Number of modes





# Parameters for lost time for a stop, Models II and following
a_m		= [1.2,    1.2,    1.4,      1.4] # Acceleration rate 					(m/s^2)
b_m		= [1.4,    1.4,    1.2,      1.1] # Deceleration rate 					(m/s^2)
t0s_m   = [2.0,    2.0,    2.0,      3.0] # Fixed time lost for a stop (doors)  (s)

# Computed parameters -- Time lost for a stop, Models II and following
def time_lost_stop(a, b, V, t0s):
# Compute time lost for acceleration, deceleration and opening and closing doors
# 	return time in hours
	return V/(2.0 * 12960.0) * (1/a + 1/b) + (t0s / 3600.0)

Tl_m = [time_lost_stop(a, b, V, t0s) for a, b, V, t0s in zip(a_m, b_m, S_m, t0s_m)]
#print 'Incremental time lost for Acc, Dec and fixed time at stop (h)', Tl_m

def stop_distance_min(a, b, V):
# Compute the shortest stop spacing to reach the design speed, return distance in km
	return V**2 / 25920.0 * (1/a + 1/b)

Dmin_m = [stop_distance_min(a, b, V) for a, b, V in zip(a_m, b_m, S_m)]
#print 'Minimal stop interdistance (km)', Dmin_m

# Parameters for fixed TU length, Models I and II
nfix_m  = [   1,    1,      1,      3]                      # Fixed number of vehicles per transit unit (veh/TU)
K_m     = [nfix_m[m] * kv_m[m] for m in range(num_modes)]   # TU capacity                           (pax/TU)

# Parameters for variable TU length, Models III and following
nmin_m  = [   1,    1,      1,      2]    # Minimal number of vehicles per transit unit (veh/TU)
nmax_m  = [   1,    1,      2,      5]    # Maximal number of vehicles per transit unit (veh/TU)


# Parameters -- Crowding penalty function, Models III and following
thetamin= 0.3                   # Avg occupancy rate up to there is no penalty
rho     = 1.0                   # Slope of the penalty  function 
xi      = 1.0 - thetamin * rho  # Derived parameter \xi



# Parameters -- peak and off-peak, Models IV and following
gamma   = 0.5          # Ratio between peak and off-peak hourly demand
chip    = 0.25          # ratio between peak hours and total service hours
chio    = 1.0 - chip    # ratio between off-peak hours and total service hours




# Style for plots, all Models
linestyle   = ['-', '--', '-.', ':']
colors      = ['brown', 'blue', 'green', 'red']
mode_label  = ['Bus', 'BRT', 'LRT', 'HR']

#Print  to stdout a latex table
from tabulate import tabulate
headers = ['Parameter', 'Unit'] + mode_label

row1    = ['Infrastructure cost', 'million $/km']   + ["%.0f" % infracost_m[m] for m in range(num_modes)] 
row2    = ['Width required by the infrastructure', 'm']   + ["%.0f" % width_m[m] for m in range(num_modes)] 
row3    = ['Stop cost', 'million $/stop']           + ["%.3f" % stopcost_m[m]    for m in range(num_modes)]
row4    = ['Infrastructure life', 'year']           + ["%.0f" % infralife_m[m] for m in range(num_modes)] 
row5    = ['Vehicle cost', 'million $/veh']           + ["%.2f" % vehcost_m[m]    for m in range(num_modes)] 
row6    = ['Vehicle technical life', 'year']           + ["%.0f" % vehlife_m[m]    for m in range(num_modes)] 

table = [row1, row2, row3, row4, row5, row6]
print 'The following is a table of the mode-specific parameters (almost) ready for latex (math mode must be corrected)'
print 'Parameters to compute fixed costs related to the infrastructure and rolling stock'
print tabulate(table, headers=headers, tablefmt="latex")


row1    = ['Crew', '$/TU-h']                        + ["%.0f" % c1t_m[m] for m in range(num_modes)]
row2    = ['Operating', '$/veh-km']                 + ["%.2f" % c2v_m[m] for m in range(num_modes)]
row3    = ['Infrastructure maintenance', '$/h']     + ["%.0f" % inframaint_m[m] for m in range(num_modes)] 
table   = [row1, row2, row3]
print 'Operating costs, including overhead'
print tabulate(table, headers=headers, tablefmt="latex")