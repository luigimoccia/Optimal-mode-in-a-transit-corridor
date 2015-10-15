# Model I
# Version 23th July 2015
# The python code solves Model I for the optimization of a transit line 
# and for the strategic assessment of mode choice, see: 
# L. Moccia, and G. Laporte, 'Improved models for technology choice in a transit corridor with fixed demand', 2015
# For comments, please write to: 
# moccia@icar.cnr.it
# Model from Tirachini et al 2010 published in TRE
# Demand fixed
# Minimization of total cost defined as the sum user and agency costs
# frequency is the only variable of this objective function



# Import formulae, and code shared by all models
from common_code import *

# Import parameters shared by all models
from common_parameters import *

# Computed cost parameters for Model I

hours_per_year  = 2947.0

c0_m = [infra_fixed_cost_per_hour(m, hours_per_year, L, width_m, land_cost_per_hectare, \
        infracost_m, infralife_m, inframaint_m, discount_rate) \
        for m in range(num_modes)]


c1_m = [nfix_m[m] * rolling_stock_cost_per_hour(m, hours_per_year, vehcost_m, vehlife_m, discount_rate, \
            res_value_rs) + c1t_m[m] \
        for m in range(num_modes)]


#Old values from Tirachini et al 2010
#c0_m    = [ 0.0, 9489.0,14866.0, 24910.0]# Unit fixed operator cost per hour   ($/h)
#c1_m    = [54.0,  63.0,   164.0,    354.0]# Unit operator cost per vehicle-hour ($/veh-h)


R_m     = [2 * L / S_m[m] for m in range(num_modes)] # Running time at speed S (h)



def Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, beta, R):
# User cost
    Ca = Pa * d * y / (2.0 * v)                                 # Cost of access
    Cw = Cost_waiting(f, y, t01, t02, t11, t12, fbar, Pw, epsilon)   # Cost of waiting
    Cv = Pv * l / (2.0 * L) * ((y/f) * beta / 3600.0  + R) * y  # Cost of in-vehicle time
    return Ca + Cw + Cv

def Co(c0, c1, y, beta, R, f, c2, L):
# Agency cost
    return c0 + c1 * (y * beta / 3600.0 + R * f) + 2.0 * c2 * L * f

def Ctot(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, beta, R, c0, c1, c2):
# Total cost sum of user and agency cost
    return Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, beta, R) + \
            Co(c0, c1, y, beta, R, f, c2, L)



def Cop_of_f(f):
# Operator cost as a function of f
    return Co(c0, c1, y, beta, R, f, c2, L)
    
def Cuser_of_f(f):
# user cost as a function of f
    return Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, beta, R)

def Ctot_of_f(f):
# Total cost as a function of f 
    return Ctot(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, beta, R, c0, c1, c2)
    
def fmin_closed(fmin, fmax):
# Compute optimal frequency by closed form
    def t_param_of_f(f):
        if (f >= fbar):
            t0 = t02
            t1 = t12
        else:
            t0 = t01
            t1 = t11
        return t0, t1
    def my_fopt():
    # Compute fopt with unpacked t1 value 
        from math import sqrt
        return sqrt((Pw * t1 * epsilon * y  +  Pv * l * y ** 2 * beta / (3600.0 * 2.0 * L)) 
                        / (c1 * R   +  2.0 * c2 * L))
    #Case 1, if f >= fbar relevant  only t02 and t12 are relevant  
    if (fmin >= fbar):
        t0, t1 = t_param_of_f(fmin)
        fopt = my_fopt()
        if fopt <= fmin:
            fopt = fmin
        elif fopt > fmax:
            fopt = fmax
    else:
        # First compute minimum for case t01 and t11
        t0, t1 = t_param_of_f(fmin)
        fopt1 = my_fopt()
        if fopt1 < fmin:
            fopt1 = fmin
        elif fopt1 > fmax:
            fopt1 = fmax
        if fopt1 > fbar:
            fopt1 = fbar - 0.01
        Copt1 = Ctot_of_f(fopt1)
        # Now compute minimum for case t02 and t12
        t0, t1 = t_param_of_f(fbar)
        fopt2 = my_fopt()
        if fopt2 < fmin:
            fopt2 = fmin
        elif fopt2 >  fmax:
            fopt2 = fmax
        if fopt2 < fbar:
            fopt2 = fbar
        Copt2 = Ctot_of_f(fopt2)
        # Determine fopt
        if Copt1 < Copt2:
            fopt = fopt1
        else:
            fopt = fopt2
            
    return fopt

def unpack_mode_specific_par(mode):
# return mode specific parameters
    return  c0_m[mode], c1_m[mode], c2_m[mode], beta_m[mode], d_m[mode], R_m[mode]

#Print  to stdout a latex table
from tabulate import tabulate
headers = ['Parameter', 'Unit'] + mode_label
row1    = ['c_0', '$/h']    + ["%.0f" % c0_m[m] for m in range(num_modes)] 
row2    = ['c_1', '$/h']    + ["%.1f" % c1_m[m] for m in range(num_modes)] 
row3    = ['c_2', '$/h']    + ["%.2f" % c2_m[m] for m in range(num_modes)] 
row4    = ['d', 'km']       + ["%.1f" % d_m[m]    for m in range(num_modes)] 
row5    = ['fmax', 'TU/h']   +["%.0f" % fmax_m[m]    for m in range(num_modes)] 
row6    = ['K', 'pax/TU']   + ["%.0f" % K_m[m]    for m in range(num_modes)] 
row7    = ['S', 'km/h']     + ["%.0f" % S_m[m]    for m in range(num_modes)] 
row8    = ['beta', 's/pax'] + ["%.2f" % beta_m[m]    for m in range(num_modes)] 
table = [row1, row2, row3, row4, row5, row6, row7, row8]
print 'The following is a table of the mode-specific parameters (almost) ready for latex (math mode must be corrected)'
print tabulate(table, headers=headers, tablefmt="latex")


# COMPUTE RESULTS
import numpy as np
from scipy.optimize import fminbound


max_demand          = np.empty(num_modes)
y_range             = [np.array([]) for m in range(num_modes)]
min_f_range         = [np.array([]) for m in range(num_modes)]
min_f_range2        = [np.array([]) for m in range(num_modes)]
f_min_range         = [np.array([]) for m in range(num_modes)]
min_avg_op_cost     = [np.array([]) for m in range(num_modes)]
min_avg_user_cost   = [np.array([]) for m in range(num_modes)]
min_avg_tot_cost    = [np.array([]) for m in range(num_modes)]

for m in range(num_modes):
    max_demand[m]       = K_m[m] * fmax_m[m] * nu / alpha
    if max_demand[m] > max_demand_studied:
        max_demand[m] = max_demand_studied
    y_range[m]          = np.arange(min_demand_studied, max_demand[m], step=step_demand)
    min_f_range[m]      = np.empty_like(y_range[m])
    min_f_range2[m]     = np.empty_like(y_range[m])
    f_min_range[m]      = np.empty_like(y_range[m])
    min_avg_tot_cost[m] = np.empty_like(y_range[m])
    min_avg_op_cost[m]  = np.empty_like(y_range[m])
    min_avg_user_cost[m]= np.empty_like(y_range[m])
    fmax                = fmax_m[m]
    c0, c1, c2, beta, d, R = unpack_mode_specific_par(m)
    for i, y in enumerate(y_range[m]): 
        fmin                     = alpha * y /(nu * K_m[m]) # Minimal frequency 
        min_f_range[m][i]        = fmin_closed(fmin, fmax)
        # Double check
        min_f_range2[m][i]       = fminbound(Ctot_of_f, x1=fmin, x2=fmax)
        from math import fabs
        if fabs(min_f_range[m][i] - min_f_range2[m][i]) > 0.0001:
             print min_f_range[m][i], min_f_range2[m][i]
        min_avg_tot_cost[m][i]  = Ctot_of_f(min_f_range[m][i])  / y
        min_avg_op_cost[m][i]   = Cop_of_f(min_f_range[m][i])   / y
        min_avg_user_cost[m][i] = Cuser_of_f(min_f_range[m][i]) / y
        f_min_range[m][i]       = fmin
    

 

# PLOTS
import matplotlib.pyplot as plt


# Print points of total cost intersection between  modes
print_itersection_points(y_range, min_avg_tot_cost, mode_label)




    # A plot with all 4 figures
fig, axes = plt.subplots(2,2) 
fig.set_size_inches(9,9)
i = 0
from math import floor
while (i < 4):
    ix = int(floor(i/2))
    iy = i - ix * 2
    axes[ix][iy].set_xlabel('Demand (pax/h)')
    axes[ix][iy].grid(ls=':', lw=0.5)
    for item in ([axes[ix][iy].title, axes[ix][iy].xaxis.label, axes[ix][iy].yaxis.label] +
             axes[ix][iy].get_xticklabels() + axes[ix][iy].get_yticklabels()):
        item.set_fontsize(10)
    for m in range(num_modes):
        if i == 0:
            axes[ix][iy].plot(y_range[m], min_f_range[m], ls=linestyle[m], lw=3, 
                                color=colors[m], label=mode_label[m])
            axes[ix][iy].set_ylabel('Frequency (veh/h)')
            axes[ix][iy].set_title('a) Optimal frequency')
        elif i == 1:
            axes[ix][iy].plot(y_range[m], min_avg_op_cost[m], ls=linestyle[m], lw=3, 
                                color=colors[m], label=mode_label[m]) 
            axes[ix][iy].set_ylabel('Avg operator cost ($/pax)')
            axes[ix][iy].set_title('b) Average operator cost')
        elif i == 2:
            axes[ix][iy].plot(y_range[m], min_avg_user_cost[m], ls=linestyle[m], lw=3, 
                                color=colors[m], label=mode_label[m])
            axes[ix][iy].set_ylabel('Avg user cost ($/pax)')
            axes[ix][iy].set_title('c) Average user cost')
        elif i == 3:
            axes[ix][iy].plot(y_range[m], min_avg_tot_cost[m], ls=linestyle[m], lw=3, 
                                color=colors[m], label=mode_label[m])
            axes[ix][iy].set_ylabel('Avg total cost ($/pax)')
            axes[ix][iy].set_title('d) Average total cost')
            
    axes[ix][iy].legend(loc='upper right')
    i += 1
name_outputfile0 = 'mod1_all_plots.pdf'
plt.savefig(name_outputfile0)
plt.close()

 

model_id = 'fig_M1_'

name_outputfile1 = model_id + 'freq.pdf'
plot_frequency(y_range, min_f_range, f_min_range, 'Demand (pax/h)', 'Frequency (TU/h)', name_outputfile1, linestyle, colors, mode_label)


name_outputfile2 = model_id + 'tot_cost.pdf'
plot_single(y_range, min_avg_tot_cost, 'Demand (pax/h)', 'Average total cost ($/pax)', name_outputfile2, linestyle, colors, mode_label)

name_outputfile3 = model_id + 'op_cost.pdf'
plot_single(y_range, min_avg_op_cost, 'Demand (pax/h)', 'Average operator cost ($/pax)', name_outputfile3, linestyle, colors, mode_label)

name_outputfile4 = model_id + 'user_cost.pdf'
plot_single(y_range, min_avg_user_cost, 'Demand (pax/h)', 'Average user cost ($/pax)', name_outputfile4, linestyle, colors, mode_label)


# A figure with all plots related to the economical aspects
name_outputfile0e = model_id + 'all_eco_plots.pdf'
plot_three_economics(y_range, min_avg_op_cost, min_avg_user_cost, min_avg_tot_cost, name_outputfile0e, linestyle, colors, mode_label)


 

    # If MAC OS open the pdf output files
import platform
if platform.system() == 'Darwin':
    import os
    command = 'open -a Preview.app ' + name_outputfile1 + ' ' +  name_outputfile2 \
        + ' ' +  name_outputfile3 + ' ' +  name_outputfile4 + ' ' + name_outputfile0 \
        + ' ' + name_outputfile0e
    os.system(command)




   # Check specific values for a mode and a level y
# m = 1
# c0, c1, c2, beta, d, R = unpack_mode_specific_par(m)
# y = 10000.0
# print 'Check specific values for mode ', m, ' and  demand ', y, ' (pax/h)'
# fmin = alpha * y /(nu * K_m[m])
# min_f = fmin_closed(fmin, max_demand[m])
# print 'max_demand[m]', max_demand[m], ' opt freq', min_f

      