# Model II
# Version 23th July 2015
# The python code solves Model II for the optimization of a transit line 
# and for the strategic assessment of mode choice, see: 
# L. Moccia, and G. Laporte, 'Improved models for technology choice in a transit corridor with fixed demand', 2015
# For comments, please write to: 
# moccia@icar.cnr.it
# Extension of the Model from Tirachini et al 2010 published in TRE
# Optimal stop spacing considering vehicle dynamics
# Demand inelastic
# Minimization of user and agency cost
# Variables:  f, d 

# Import formulae, and code shared by all models
from common_code import *

# Import parameters shared by all models
from common_parameters import *

# Computed cost parameters for Model II
hours_per_year  = 2947.0

c0l_m = [infra_fixed_cost_per_hour(m, hours_per_year, L, width_m, land_cost_per_hectare, \
        infracost_m, infralife_m, inframaint_m, discount_rate) \
        for m in range(num_modes)]

c0s_m = [stop_fixed_cost_per_hour(m, hours_per_year, stopcost_m, infralife_m, discount_rate) \
         for m in range(num_modes)]

c1_m = [nfix_m[m] * rolling_stock_cost_per_hour(m, hours_per_year, vehcost_m, vehlife_m, discount_rate, \
            res_value_rs) + c1t_m[m] \
        for m in range(num_modes)]

R_m     = [2 * L / Vmax_m[m] for m in range(num_modes)] # Running time at speed Vmax (h)

def Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, beta, R, Tl):
# User cost
    Ca = Pa * d * y / (2.0 * v)                                     # Cost of access
    Cw = Cost_waiting(f, y, t01, t02, t11, t12, fbar, Pw, epsilon)  # Cost of waiting
    Cv = Pv * l / (2.0 * L) * ((y/f) * beta / 3600.0  + R  + 2.0 * L / d * Tl) * y  # Cost of in-vehicle time
    return Ca + Cw + Cv


def Co(c0l, c0s, c1, y, beta, R, f, c2, L, d, Tl):
# Agency cost
    return c0l +  2.0 * c0s * L/d +  c1 * (y * beta / 3600.0 + (R + 2 * L / d * Tl) * f) + 2.0 * c2 * L * f

def Ctot(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, beta, R, c0l, c0s, c1, c2, Tl):
# Total cost sum of user and agency cost
    return Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, beta, R, Tl) + \
            Co(c0l, c0s, c1, y, beta, R, f, c2, L, d, Tl)

def coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, f, l, L, beta, R, c0l, c0s, c1, c2, Tl):
# Coefficients a0, a1, a2, a3, a4, a5
    if (f >= fbar):
        a0 = Pv * l / (2.0 * L) * R * y + c0l + c1 * beta/3600.0 * y   
        a4 = Pw * epsilon * y  +  Pv * l / (2.0 * L) * beta/3600.0 * y * y
    else:
        a0 = Pw * tw / 60.0 * y + Pv * l / (2.0 * L) * R * y + c0l + c1 * beta/3600.0 * y  
        a4 = Pw * mu * epsilon * y  +  Pv * l / (2.0 * L) * beta/3600.0 * y * y
    a1 = Pa * y / (2.0 * v)
    a2 = Pv * l * Tl * y + 2.0 * c0s * L
    a3 = c1 * R  +  2.0 * c2 * L
    a5 = 2 * c1 * L * Tl
    return a0, a1, a2, a3, a4, a5

def determinant(f, d, a0, a1, a2, a3, a4, a5):
# Compute determinant of the Hessian
    det = (4 * a4 * ((a5 * f + a2) / (d * d * d))) / (f * f * f) - (a5 * a5)/(d * d * d * d)
    return det

def Cop_of_f_and_d(vars):
# Operator cost as a function of f and d
    f, d = vars
    return Co(c0l, c0s, c1, y, beta, R, f, c2, L, d, Tl)
    
def Cuser_of_f_and_d(vars):
# user cost as a function of f and d
    f, d = vars
    return Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, beta, R, Tl)

def Ctot_of_f_and_d(vars):
# Total cost as a function of f and d
    f, d = vars
    return Ctot(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, beta, R, c0l, c0s, c1, c2, Tl)

def unpack_mode_specific_par(mode):
# return mode specific parameters
    return  c0l_m[mode], c0s_m[mode], c1_m[mode], c2_m[mode], beta_m[mode], R_m[mode], Tl_m[m]

def solve_heuristicaly(fmin, fmax, dmin, dmax):
# return heuristic solution of the lower convex approximation of the objective function
    def compute_my_f_d(myfmin, myfmax):
    # Compute the constrained approximated f and d in the range myfmin, myfmax
        from math import sqrt
        myd = sqrt((a2 + a5 * myfmin)/float(a1))
        if myd < dmin:
            myd = dmin
        elif myd > dmax:
            myd = dmax
        myf = sqrt(a4/float(a3))
        myf_unc = myf
        if myf < myfmin:
            myf = myfmin
        elif myf > myfmax:
            myf = myfmax
        myv = (a0 + a1 * myd + a2/float(myd) + a3 * myf + a4/float(myf) + a5 * myfmin/float(myd)) / y
        return myv, myf, myd, myf_unc
            
    if (fmin >= fbar):
    # Case 1, only case b is relevant
        a0, a1, a2, a3, a4, a5 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fmin, l, L, beta, R, c0l, c0s, c1, c2, Tl)
        value, f_approx, d_approx, f_unc = compute_my_f_d(fmin, fmax)
    else:
    #Case 2, if fmin < fbar then also case a could be relevant
        #Case a
        fsupa = fbar
        if fmax < fbar:
            fsupa = fmax
        a0, a1, a2, a3, a4, a5 = \
            coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fmin, l, L, beta, R, c0l, c0s, c1, c2, Tl)
        valuea, f_approxa, d_approxa, f_unca = compute_my_f_d(fmin, fsupa)
        # Case b
        if fmax >= fbar:
            a0, a1, a2, a3, a4, a5 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fbar, l, L, beta, R, c0l, c0s, c1, c2, Tl)
            valueb, f_approxb, d_approxb, f_uncb = compute_my_f_d(fbar, fmax)
            if valuea < valueb:
                value       = valuea
                f_approx    = f_approxa
                d_approx    = d_approxa
                f_unc       = f_unca
            else:
                value       = valueb
                f_approx    = f_approxb
                d_approx    = d_approxb
                f_unc       = f_uncb
        else:
            value       = valuea
            f_approx    = f_approxa
            d_approx    = d_approxa
            f_unc       = f_unca
            
    return value, f_approx, d_approx, f_unc
  

#Print  to stdout a latex table
from tabulate import tabulate
headers = ['Parameter', 'Unit'] + mode_label
row1    = ['bar a', 'm/s^2']    + ["%.1f" % a_m[m] for m in range(num_modes)] 
row2    = ['bar b', 'm/s^2']    + ["%.1f" % b_m[m] for m in range(num_modes)] 
row3    = ['c_0l', '$/h']    +    ["%.0f" % c0l_m[m] for m in range(num_modes)] 
row4    = ['c_0s', '$/h']       + ["%.1f" % c0s_m[m] for m in range(num_modes)] 
row5    = ['td', 's']           + ["%.0f" % t0s_m[m]    for m in range(num_modes)] 
row6    = ['Tl', 'h']           + ["%.4f" % Tl_m[m]    for m in range(num_modes)] 
row7    = ['Tl', 's']           + ["%.0f" % (Tl_m[m] * 3600)    for m in range(num_modes)]
row8    = ['dmin', 'km']        + ["%.2f" % (Dmin_m[m])    for m in range(num_modes)]

table = [row1, row2, row3, row4, row5, row6, row7, row8]
print 'The following is a table of the mode-specific parameters (almost) ready for latex (math mode must be corrected)'
print ' These parameters are specific of Model II, other parameters are as in Model I'
print tabulate(table, headers=headers, tablefmt="latex")



# SOLVE
import numpy as np
from scipy.optimize import minimize

    # data structures to track results
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
opt_gap             = [np.array([]) for m in range(num_modes)]
heur_gap            = [np.array([]) for m in range(num_modes)]
number_iterations   = [np.array([]) for m in range(num_modes)]
f_unc_range         = [np.array([]) for m in range(num_modes)]

for m in range(num_modes):
    max_demand[m]       = K_m[m] * fmax_m[m] * nu / alpha
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
    opt_gap[m]	        = np.empty_like(y_range[m])
    heur_gap[m]	        = np.empty_like(y_range[m])
    number_iterations[m]= np.empty_like(y_range[m])
    f_unc_range[m]      = np.empty_like(y_range[m])
    fmax                = fmax_m[m]
    dmin                = Dmin_m[m]
    dmax                = 2.0
    c0l, c0s, c1, c2, beta, R, Tl = unpack_mode_specific_par(m)
    import time
    tot_timespent       = 0.0
    for i, y in enumerate(y_range[m]): 
        fmin                    = alpha * y /(nu * K_m[m]) # Minimal frequency
        value_lower, f_heur, d_heur, f_unc = solve_heuristicaly(fmin, fmax, dmin, dmax)
        value_heur              =  Ctot_of_f_and_d((f_heur, d_heur))  / y
        bnds 					= ((fmin, fmax), (dmin, dmax))
        # x0 				    = np.array((fmax, dmax))
        # Use heuristic solution as starting point for L-BFGS-B algorithm
        x0 						= np.array((f_heur, d_heur))
        timestart = time.clock()
        res = minimize(Ctot_of_f_and_d, x0=x0, bounds=bnds, method='L-BFGS-B')
        tot_timespent += time.clock() - timestart
        minf = res.x[0]
        mind = res.x[1]
        value =  Ctot_of_f_and_d(res.x)  / y
        min_f_range[m][i]       = minf
        min_d_range[m][i]       = mind
        f_min_range[m][i]       = fmin
        min_avg_tot_cost[m][i]  = value
        min_avg_op_cost[m][i]   = Cop_of_f_and_d(res.x)   / y
        min_avg_user_cost[m][i] = Cuser_of_f_and_d(res.x) / y
        commercial_speed[m][i]  = 2.0 * L / (y/minf * beta / 3600.0 + R  + 2.0 * L * Tl / mind)
        avg_running_speed[m][i] = 2.0 * L / (R  + 2.0 * L * Tl / mind)
        opt_gap[m][i]           = (value - value_lower) / value * 100
        heur_gap[m][i]          = (value_heur - value) / value * 100
        number_iterations[m][i] = res.nit
        f_unc_range[m][i]       = f_unc
    # Print some computational details
    print_some_computational_details(mode_label[m], opt_gap[m], heur_gap[m], number_iterations[m], tot_timespent)
    

# Print a comparison of ratios of avg stop spacing
print_ratios_avg_stop_spacing(min_d_range, Tl_m)

# Print points of total cost intersection between  modes
print_itersection_points(y_range, min_avg_tot_cost, mode_label)

# Print ratio between f_min and f_unc
# for m in range(num_modes):
#     print 'Mode', m
#     print f_unc_range[m]/f_min_range[m]


# PLOT
import matplotlib.pyplot as plt



   # A figure with all plots
fig, axes = plt.subplots(3,2) 
fig.set_size_inches(8,12)
from math import floor
for i in range(6):
    ix = int(floor(i/2))
    iy = i - ix * 2
    axes[ix][iy].set_xlabel('Demand (pax/h)')
    axes[ix][iy].grid(ls=':', lw=0.5)
    for item in ([axes[ix][iy].title, axes[ix][iy].xaxis.label, axes[ix][iy].yaxis.label] +
             axes[ix][iy].get_xticklabels() + axes[ix][iy].get_yticklabels()):
        item.set_fontsize(8)
    for m in range(num_modes):
        if i == 0:
            axes[ix][iy].plot(y_range[m], min_f_range[m], ls=linestyle[m], lw=3, 
                                color=colors[m], label=mode_label[m])
            axes[ix][iy].set_ylabel('Frequency (veh/h)')
            axes[ix][iy].set_title('a) Optimal frequency')
        elif i == 1:
            axes[ix][iy].plot(y_range[m], min_d_range[m], ls=linestyle[m], lw=3,
            					color=colors[m], label=mode_label[m])
            axes[ix][iy].set_ylabel('Distance (km)')
            axes[ix][iy].set_title('b) Optimal stop distance')
        elif i == 2:
            axes[ix][iy].plot(y_range[m], commercial_speed[m], ls=linestyle[m], lw=3,
            					color=colors[m], label=mode_label[m])
            axes[ix][iy].set_ylabel('Commercial speed (km/h)')
            axes[ix][iy].set_title('c) Commercial speed')
        elif i == 3:
            axes[ix][iy].plot(y_range[m], min_avg_op_cost[m], ls=linestyle[m], lw=3, 
                                color=colors[m], label=mode_label[m]) 
            axes[ix][iy].set_ylabel('Avg operator cost ($/pax)')
            axes[ix][iy].set_title('d) Average operator cost')
        elif i == 4:
            axes[ix][iy].plot(y_range[m], min_avg_user_cost[m], ls=linestyle[m], lw=3, 
                                color=colors[m], label=mode_label[m])
            axes[ix][iy].set_ylabel('Avg user cost ($/pax)')
            axes[ix][iy].set_title('e) Average user cost')
        elif i == 5:
            axes[ix][iy].plot(y_range[m], min_avg_tot_cost[m], ls=linestyle[m], lw=3, 
                                color=colors[m], label=mode_label[m])
            axes[ix][iy].set_ylabel('Avg total cost ($/pax)')
            axes[ix][iy].set_title('f) Average total cost')   
    axes[ix][iy].legend(loc='upper right', fancybox=True, framealpha=0.5, fontsize=10)   
name_outputfile0 = 'mod2_all_plots.pdf'
plt.savefig(name_outputfile0)
plt.close()


model_id = 'fig_M2_'

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


name_outputfile5 = model_id + 'stop_dist.pdf'
plot_single(y_range, min_d_range, 'Demand (pax/h)', 'Stop spacing (km)', name_outputfile5, linestyle, colors, mode_label)


name_outputfile6 = model_id + 'commercial_speed.pdf'
plot_single(y_range, commercial_speed, 'Demand (pax/h)', 'Commercial speed (km/h)', name_outputfile6, linestyle, colors, mode_label)

name_outputfile7 = model_id + 'avg_running_speed.pdf'
plot_single(y_range, avg_running_speed, 'Demand (pax/h)', 'Avg running speed (km/h)', name_outputfile7, linestyle, colors, mode_label)


# With MAC open the pdf output files
import platform
if platform.system() == 'Darwin':
    import os
    command = 'open -a Preview.app ' \
        +        name_outputfile1  + ' ' +  name_outputfile2 \
        + ' ' +  name_outputfile3  + ' ' +  name_outputfile4 \
        + ' ' +  name_outputfile0e + ' ' +  name_outputfile5 \
        + ' ' +  name_outputfile6  + ' ' +  name_outputfile0 \
        + ' ' +  name_outputfile7
    os.system(command)



# #Check determinant Hessian
# for m in range(num_modes):
#     max_demand[m]       = K_m[m] * fmax_m[m] * nu / alpha
#     if max_demand[m] > max_demad_studied:
#         max_demand[m] = max_demad_studied
#     y_range[m]          = np.arange(min_demand_studied, max_demand[m], step=500)
#     fmax                = fmax_m[m]
#     c0l, c0s, c1, c2, beta, R, Tl = unpack_mode_specific_par(m)
#     for i, y in enumerate(y_range[m]): 
#         d = Dmin_m[m]
#         f = fmax
#         a0, a1, a2, a3, a4, a5 = \
#            coefficients(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, f, l, L, beta, R, c0l, c0s, c1, c2, Tl)
#         print 'a4, a5 ', a4, a5
#         det = determinant(f, d, a0, a1, a2, a3, a4, a5)
#         print 'det ', det, ' mode ', m, ' demand ', y
    
   
# Check specific values for a mode and a level y - TEST only
# m = 1
# c0l, c0s, c1, c2, beta, R, Tl = unpack_mode_specific_par(m)
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
#     name_test_file = 'mod2test' + str(m) + '.txt'
#     out_file = open(name_test_file,"w")
#     savetxt(out_file, min_avg_tot_cost[m], fmt='%.2f')
#     out_file.close()

