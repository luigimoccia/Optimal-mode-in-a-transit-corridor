# Model III
# Version 23th July 2015
# The python code solves Model III for the optimization of a transit line 
# and for the strategic assessment of mode choice, see: 
# L. Moccia, and G. Laporte, 'Improved models for technology choice in a transit corridor with fixed demand', 2015
# For comments, please write to: 
# moccia@icar.cnr.it
# Extension of the Model from Tirachini et al 2010 published in TRE
# Optimal stop spacing considering vehicle dynamics
# Variable train consist
# Demand inelastic
# Minimization of user and agency cost,  f, d, n as variables


# Import formulae, and code shared by all models
from common_code import *

# Import parameters shared by all models
from common_parameters import *


# Computed cost parameters for Model III
hours_per_year  = 2947.0 # Equivalent number of peak hours per year - Tirachini et al 2010


c0l_m = [infra_fixed_cost_per_hour(m, hours_per_year, L, width_m, land_cost_per_hectare, \
        infracost_m, infralife_m, inframaint_m, discount_rate) \
        for m in range(num_modes)]

c0s_m = [stop_fixed_cost_per_hour(m, hours_per_year, stopcost_m, infralife_m, discount_rate) \
         for m in range(num_modes)]

c1v_m = [rolling_stock_cost_per_hour(m, hours_per_year, vehcost_m, vehlife_m, discount_rate, \
            res_value_rs) \
        for m in range(num_modes)]
        
R_m     = [2 * L / Vmax_m[m] for m in range(num_modes)] # Running time at speed Vmax (h)


def Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, betav, R, Tl, n, kv):
# User cost
    Ca = Pa * d * y / (2.0 * v)                         # Cost of access
    Cw = Cost_waiting(f, y, t01, t02, t11, t12, fbar, Pw, epsilon)      # Cost of waiting
    delta = penalty_crowding(f, n, y, l, L, kv, xi, rho, thetamin)
    ct = cycle_time(f, n, d, y, betav, R, L, Tl)
    Cv = Pv * l / (2.0 * L) * delta * ct * y            # Cost of in-vehicle time
    return Ca + Cw + Cv
    

def Co(c0l, c0s,c1t, c1v, y, betav, R, f, c2v, L, d, Tl, n):
# Agency cost
	cost0 = c0l + 2.0 * c0s* L /d
	cost1 = c1t * (y * betav / (3600.0 * float(n)) + (R + 2.0 * L / d * Tl) * f)  
	cost2 =  c1v * n * (y * betav / (3600.0 * float(n)) + (R + 2.0 * L / d * Tl) * f)
	cost3 =  2.0 * n * c2v * L * f
	cost = cost0 + cost1 + cost2 + cost3
	return cost

def Ctot(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, betav, R, c0l, c0s,c1t, c1v, c2v, Tl, n, kv):
# Total cost sum of user and agency cost
    return Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, betav, R, Tl, n, kv) + \
            Co(c0l, c0s,c1t, c1v, y, betav, R, f, c2v, L, d, Tl, n)



def Cop_of_f_and_d(vars):
# Operator cost as a function of f and d
    f, d = vars
    return Co(c0l, c0s,c1t, c1v, y, betav, R, f, c2v, L, d, Tl, n)
    
def Cuser_of_f_and_d(vars):
# user cost as a function of f and d
    f, d = vars
    return Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, betav, R, Tl, n, kv)

def Ctot_of_f_and_d(vars):
# Total cost as a function of f and d
    f, d = vars
    return Ctot(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, betav, R, c0l, c0s,c1t, c1v, c2v, Tl, n, kv)

def unpack_mode_specific_par(mode):
# return mode specific parameters
    return  c0l_m[mode], c0s_m[mode], c1t_m[mode], c1v_m[mode], c2v_m[mode], betav_m[mode], R_m[mode], Tl_m[m], \
            nmin_m[m], nmax_m[m], kv_m[m]

def coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, f, l, L, betav, R, c0l, c0s,c1t, c1v, c2v, Tl, n, kv, xi, rho):
# Coefficients a0, a1,...,a11
    if (f >= fbar):
        a0 = Pv * l * xi/ (2.0 * L) * R * y + c0l + c1v * betav/3600.0 * y   
        a4 = Pw * epsilon * y  
    else:
        a0 = Pw * tw / 60.0 * y + Pv * l * xi/ (2.0 * L) * R * y + c0l + c1v * betav/3600.0 * y  
        a4 = Pw * mu * epsilon * y 
    a1  = Pa * y / (2.0 * v)
    a2  = Pv * l * Tl * xi *  y  + 2.0 * c0s * L
    a3  = c1t * R  
    a5  = 2 * c1t * L * Tl
    a6  = c1t * betav/3600.0 * y
    a7  = c1v * R +  2.0 * c2v * L
    a8  = Pv  * l * y *  (xi * betav * y /(7200.0 * L)   +  rho * R * l * y/ (4.0 * L * L * kv) )
    a9  = 2 * c1v * L * Tl
    a10 = Pv * rho * l * l * y * y * Tl /(2.0 * L * kv)
    a11 = Pv * rho * l * l * y * y * y * betav / (14400.0 * L * L * kv)
    
    return a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11
    
    
def solve_heuristicaly(fmin, fmax, dmin, dmax, xi, rho):
# return heuristic solution of the lower convex approximation of the objective function
    def compute_my_f_d(myfmin, myfmax):
    # Compute the constrained approximated f and d in the range myfmin, myfmax
        from math import sqrt
        myd = sqrt((a2 + myfmin * (a5 + n * a9) + a10 / float(n * myfmax)) /float(a1))
        if myd < dmin:
            myd = dmin
        elif myd > dmax:
            myd = dmax
        myf = sqrt((a4 + a8/float(n) + a11/float(n * n * myfmax))/float(a3  +  n * a7))
        if myf < myfmin:
            myf = myfmin
        elif myf > myfmax:
            myf = myfmax
        myv = (a0 + a1 * myd + a2/float(myd) + a3 * myf + a4/float(myf) \
                + a5 * myfmin/float(myd) + a6 /float(n) + a7 * n * myf + a8/float(n * myf) \
                + a9 * n * myfmin /myd + a10 /float(n * myfmax * myd) + a11/float(n * n * myfmax * myf) ) / y
        return myv, myf, myd
            
    if (fmin >= fbar):
    # Case 1, only case b is relevant
        a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fmin, l, L, betav, R, \
                                c0l, c0s,c1t, c1v, c2v, Tl, n, kv, xi, rho)
        value, f_approx, d_approx = compute_my_f_d(fmin, fmax)
    else:
    #Case 2, if fmin < fbar, then also case a could be relevant
        #Case a
        fsupa = fbar
        if fmax < fbar:
            fsupa = fmax
        a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fmin, l, L, betav, R, \
                                c0l, c0s,c1t, c1v, c2v, Tl, n, kv, xi, rho)
        valuea, f_approxa, d_approxa = compute_my_f_d(fmin, fsupa)
        # Case b, only to be evaluated if fmax >= fbar
        if fmax >= fbar:
            a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fbar, l, L, betav, R, \
                                c0l, c0s,c1t, c1v, c2v, Tl, n, kv, xi, rho)
            valueb, f_approxb, d_approxb = compute_my_f_d(fbar, fmax)
            if valuea <= valueb:
                value       = valuea
                f_approx    = f_approxa
                d_approx    = d_approxa
            else:
                value       = valueb
                f_approx    = f_approxb
                d_approx    = d_approxb
            
        else:
            value       = valuea
            f_approx    = f_approxa
            d_approx    = d_approxa
              
    return value, f_approx, d_approx
  

#Print  to stdout a latex table
from tabulate import tabulate
headers = ['Parameter', 'Unit'] + mode_label
row1    = ['c_0l', '$/h']    +    ["%.0f" % c0l_m[m] for m in range(num_modes)] 
row2    = ['c_0s', '$/h']       + ["%.1f" % c0s_m[m] for m in range(num_modes)] 
row3    = ['c1t', '$/h']        + ["%.2f" % c1t_m[m] for m in range(num_modes)]
row4    = ['c1v', '$/h']        + ["%.2f" % c1v_m[m] for m in range(num_modes)]   
row5    = ['c2v', '$/h']        + ["%.2f" % c2v_m[m] for m in range(num_modes)] 
row6    = ['kv_m', 'pax/veh']   + ["%.0f" % kv_m[m]    for m in range(num_modes)] 
row7    = ['nmin_m', 'veh']     + ["%.0f" % nmin_m[m]    for m in range(num_modes)] 
row8    = ['nmax_m', 'veh']     + ["%.0f" % nmax_m[m]    for m in range(num_modes)]
row9    = ['betav', 's/veh']    + ["%.2f" % betav_m[m]    for m in range(num_modes)]

table = [row1, row2, row3, row4, row5, row6, row7, row8, row9]
print 'The following is a table of the mode-specific parameters (almost) ready for latex (math mode must be corrected)'
print ' These parameters are specific of Model III, other parameters are as in Models I and II'
print tabulate(table, headers=headers, tablefmt="latex")
    

# SOLVE
import numpy as np
from scipy.optimize import minimize
max_demand          = np.empty(num_modes)
y_range             = [np.array([]) for m in range(num_modes)]
min_f_range         = [np.array([]) for m in range(num_modes)]
f_min_range         = [np.array([]) for m in range(num_modes)]
min_d_range         = [np.array([]) for m in range(num_modes)]
min_avg_op_cost     = [np.array([]) for m in range(num_modes)]
min_avg_user_cost   = [np.array([]) for m in range(num_modes)]
min_avg_tot_cost    = [np.array([]) for m in range(num_modes)]
commercial_speed    = [np.array([]) for m in range(num_modes)]
avg_running_speed   =  [np.array([]) for m in range(num_modes)]
vehicle_per_TU      = [np.array([]) for m in range(num_modes)]
opt_gap             = [np.array([]) for m in range(num_modes)]
heur_gap            = [np.array([]) for m in range(num_modes)]
number_iterations   = [np.array([]) for m in range(num_modes)]
f_heur_range        = [np.array([]) for m in range(num_modes)]
f_max_range         = [np.array([]) for m in range(num_modes)]

for m in range(num_modes):
    max_demand[m]       = kv_m[m] * nmax_m[m] * fmax_m[m] * nu / alpha
    if max_demand[m] > max_demand_studied:
        max_demand[m] = max_demand_studied
    y_range[m]          = np.arange(min_demand_studied, max_demand[m], step=step_demand)
    min_f_range[m]      = np.empty_like(y_range[m])
    f_min_range[m]      = np.empty_like(y_range[m])
    min_d_range[m]      = np.empty_like(y_range[m])
    min_avg_tot_cost[m] = np.empty_like(y_range[m])
    min_avg_op_cost[m]  = np.empty_like(y_range[m])
    min_avg_user_cost[m]= np.empty_like(y_range[m])
    commercial_speed[m]	= np.empty_like(y_range[m])
    avg_running_speed[m]= np.empty_like(y_range[m])
    vehicle_per_TU[m]	= np.empty_like(y_range[m])
    opt_gap[m]	        = np.empty_like(y_range[m])
    heur_gap[m]	        = np.empty_like(y_range[m])
    number_iterations[m]= np.empty_like(y_range[m])
    f_heur_range[m]     = np.empty_like(y_range[m])
    f_max_range[m]      = np.empty_like(y_range[m])
    fmax                = fmax_m[m]
    dmin                = Dmin_m[m]
    dmax                = 2.0
    c0l, c0s,c1t, c1v, c2v, betav, R, Tl, nmin, nmax, kv = unpack_mode_specific_par(m)
    import time
    tot_timespent       = 0.0
    for i, y in enumerate(y_range[m]):
        best_value = 10000.0
        value_lower_foralln = 10000.0
        for n in range(nmin, nmax + 1):  
            fmin                    = alpha * y /(nu * kv * float(n)) # Minimal frequency
            ft = fmax
            max_ft = max_freq_theta(n, y, l, L, kv, thetamin)
            if max_ft > fmax:
                max_ft = fmax
            
            #Compute only if this value of n results in a feasible range for freq 
            if fmax >= fmin:  
                #First, evaluate for the linear part of the crowding penalty
                value_lower, f_heur, d_heur = solve_heuristicaly(fmin, max_ft, dmin, dmax, xi, rho)
                value_heur              =  Ctot_of_f_and_d((f_heur, d_heur))  / y
                #Second, evaluate without the crowding penalty if a feasible f range (max_ft, fmax) exists
                if fmax > max_ft:
                    value_lower2, f_heur2, d_heur2 = solve_heuristicaly(max_ft, fmax, dmin, dmax, 1.0, 0.0)
                    value_heur2              =  Ctot_of_f_and_d((f_heur2, d_heur2))  / y
                    if value_heur2 < value_heur:
                        value_lower = value_lower2
                        f_heur      = f_heur2
                        d_heur      = d_heur2
    
                bnds 					= ((fmin, fmax), (dmin, dmax))
                #starting point
                f_start = f_heur
                d_start = d_heur
                # Use heuristic solution as starting point for L-BFGS-B algorithm
                x0 						= np.array((f_start, d_start))
                timestart = time.clock()
                res = minimize(Ctot_of_f_and_d, x0=x0, bounds=bnds, method='L-BFGS-B')
                tot_timespent += time.clock() - timestart
                value = Ctot_of_f_and_d(res.x)  / y
                if value < best_value:
                    min_f_range[m][i]       = res.x[0]
                    min_d_range[m][i]       = res.x[1]
                    f_min_range[m][i]       = fmin
                    min_avg_tot_cost[m][i]  = value
                    min_avg_op_cost[m][i]   = Cop_of_f_and_d(res.x)   / y
                    min_avg_user_cost[m][i] = Cuser_of_f_and_d(res.x) / y
                    ct		                = cycle_time(res.x[0], n, res.x[1], y, betav, R, L, Tl) 
                    commercial_speed[m][i]  = 2.0 * L / (ct)
                    avg_running_speed[m][i] = 2.0 * L / (R  + 2.0 * L * Tl / res.x[1])
                    vehicle_per_TU[m][i]    = n
                    if value_lower < value_lower_foralln:
                        value_lower_foralln = value_lower
                    opt_gap[m][i]           = (value - value_lower_foralln) / value * 100
                    heur_gap[m][i]          = (value_heur - value) / value * 100
                    f_heur_range[m][i]      = f_heur
                    f_max_range[m][i]       = max_ft
                    number_iterations[m][i] = res.nit 
                    best_value = value
    # Print some computational details
    print_some_computational_details(mode_label[m], opt_gap[m], heur_gap[m], number_iterations[m], tot_timespent/nmax)
    



# Print a comparison of ratios of avg stop spacing
print_ratios_avg_stop_spacing(min_d_range, Tl_m)



# Print points of total cost intersection between  modes
print_itersection_points(y_range, min_avg_tot_cost, mode_label)

# PLOT
import matplotlib.pyplot as plt


model_id = 'fig_M3_'

name_outputfile1 = model_id + 'freq.pdf'
plot_frequency(y_range, min_f_range, f_min_range, 'Demand (pax/h)', 'Frequency (TU/h)', name_outputfile1, linestyle, colors, mode_label)

name_outputfile1max = model_id + 'freq_max.pdf'
plot_frequency_max(y_range, min_f_range, f_max_range, 'Demand (pax/h)', 'Frequency (TU/h)', name_outputfile1max, linestyle, colors, mode_label)


name_outputfile2 = model_id + 'tot_cost.pdf'
plot_single(y_range, min_avg_tot_cost, 'Demand (pax/h)', 'Average total cost ($/pax)', name_outputfile2, linestyle, colors, mode_label)

name_outputfile3 = model_id + 'op_cost.pdf'
plot_single(y_range, min_avg_op_cost, 'Demand (pax/h)', 'Average operator cost ($/pax)', name_outputfile3, linestyle, colors, mode_label)

name_outputfile4 = model_id + 'user_cost.pdf'
plot_single(y_range, min_avg_user_cost, 'Demand (pax/h)', 'Average user cost ($/pax)', name_outputfile4, linestyle, colors, mode_label)


# A figure with all plots related to the economical aspects
name_outputfile0e = model_id + 'all_eco_plots.pdf'
plot_three_economics(y_range, min_avg_op_cost, min_avg_user_cost, min_avg_tot_cost, name_outputfile0e, linestyle, colors, mode_label)


name_outputfile5 = model_id + 'stop_dist.pdf'
plot_single(y_range, min_d_range, 'Demand (pax/h)', 'Stop spacing (km)', name_outputfile5, linestyle, colors, mode_label)

name_outputfile6 = model_id + 'commercial_speed.pdf'
plot_single(y_range, commercial_speed, 'Demand (pax/h)', 'Commercial speed (km/h)', name_outputfile6, linestyle, colors, mode_label)

name_outputfile7 = model_id + 'avg_running_speed.pdf'
plot_single(y_range, avg_running_speed, 'Demand (pax/h)', 'Avg running speed (km/h)', name_outputfile7, linestyle, colors, mode_label)


name_outputfile8 = model_id + 'veh_TU.pdf'
fig, axes = plt.subplots(1,1, figsize=(xfigsize, yfigsize))   
axes.set_xlabel('Demand (pax/h)')
axes.set_ylabel('Vehicles per transit unit (veh/TU)')
from matplotlib.ticker import MaxNLocator
axes.yaxis.set_major_locator(MaxNLocator(integer=True))
#axes.grid()
for m in range(num_modes):
    if (nmax_m[m] > 1):
        axes.plot(y_range[m], vehicle_per_TU[m], ls=linestyle[m], lw=3,color=colors[m], label= mode_label[m])
axes.legend(loc='best', fancybox=True, framealpha=0.5, fontsize=11)
plt.savefig(name_outputfile8)
plt.close()



# With MAC open the pdf output files
import platform
if platform.system() == 'Darwin':
    import os
    command = 'open -a Preview.app ' \
        +        name_outputfile1  + ' ' + name_outputfile1max + ' ' +  name_outputfile2 \
        + ' ' +  name_outputfile3  + ' ' +  name_outputfile4 \
        + ' ' +  name_outputfile0e + ' ' +  name_outputfile5 \
        + ' ' +  name_outputfile6  + ' ' +  name_outputfile7  \
        + ' ' +  name_outputfile8
    os.system(command)






   
# Check specific values for a mode and a level y - TEST only
# m = 1
# c0l, c0s,c1, c2, betav, R, Tl = unpack_mode_specific_par(m)
# y = 20000
# print 'Check specific values for mode ', m, ' and  demand ', y, ' (pax/h)'
# fmin = alpha * y /(nu * K_m[m])
# fmax = fmax_m[m]
# bnds = ((fmin, fmax), (0.5, 2.0))
# x0 = np.array((fmax, 2.0))
# res = minimize(Ctot_of_f_and_d, x0=x0, bounds=bnds, method='L-BFGS-B')
# print 'max_demand[m]', m, ' opt freq', res.x[0], ' opt d', res.x[1]   


# Save on file some results
# from numpy import savetxt
# for m in range(num_modes):
#     name_test_file = 'mod3bistest' + str(m) + '.txt'
#     out_file = open(name_test_file,"w")
#     savetxt(out_file, np.c_[heur_gap[m], min_f_range[m], f_heur_range[m], f_min_range[m]], fmt='%.2f  %4d %4d  %4d')
#     out_file.close()
