# Common code

import numpy as np


# Common formula, All Models
def Cost_waiting(f, y, t01, t02, t11, t12, fbar, Pw, epsilon):
# Cost of waiting
    if (f >= fbar):
        t0 = t02
        t1 = t12
    else:
        t0 = t01
        t1 = t11
    Cw = Pw * (t0 + t1 * (epsilon/f)) * y                      
    return Cw


# Common formulae, Models III and IV
def cycle_time(f, n, d, y, betav, R, L, Tl):
# Compute the cycle time
    tc = (y/f) * betav / (3600.0 * float(n)) + R  + 2.0 * L / d * Tl
    return tc

def max_freq_theta(n, y, l, L, kv, thetamin):
#Max freq to have  theta >= thetamin
    ft = l * y /(2.0 * L * thetamin * kv * float(n))
    return ft

def penalty_crowding(f, n, y, l, L, kv, xi, rho, thetamin):
# Return cost of crowding
    ft = max_freq_theta(n, y, l, L, kv, thetamin) 
    if f <= ft:
        delta = xi + rho * (l * y /(n * f * 2.0 * L * kv))  # Crowding penalty
    else:
        delta = 1.0
    return delta

# Common accounting formulae, All Models
def amortization_factor(discount_rate, technical_life):
    # Compute amortization factor
    from math import pow
    if technical_life > 0:
        af = discount_rate/(1.0 - 1.0/pow(1.0 + discount_rate, technical_life))
    else:
        af = 0.0
    return af

def infra_fixed_cost_per_hour(mode, hours_per_year, L, width_m, land_cost_per_hectare, infracost_m, infralife_m, inframaint_m, discount_rate):
    # Compute infrastructure fixed cost per equivalent hour of service in $/h
    area = L * 1000 * width_m[mode]                                 #(mq)
    value_land = area / 10000.0 * land_cost_per_hectare * 1000000.0 #($)
    value_line = L * infracost_m[mode] * 1000000.0                  #($)
    tot_value =   (value_land + value_line) 
    af = amortization_factor(discount_rate, infralife_m[mode])
    perhour = tot_value * af / hours_per_year
    perhour += inframaint_m[mode]
   # print m, value_line * af/hours_per_year, value_land * af/hours_per_year
    return perhour

def stop_fixed_cost_per_hour(mode, hours_per_year, stopcost_m, infralife_m, discount_rate):
    # Compute stop fixed cost per equivalent hour of service in $/h
    value_stop = stopcost_m[mode] * 1000000.0                     #($)
    af = amortization_factor(discount_rate, infralife_m[mode])
    perhour = value_stop * af / hours_per_year
    return perhour

def rolling_stock_cost_per_hour(mode, hours_per_year, vehcost_m, vehlife_m, discount_rate, res_value_rs):
    # Compute the rolling cost per equivalent hour of service in $/h
    af = amortization_factor(discount_rate, vehlife_m[mode])
    perhour = vehcost_m[mode] * 1000000.0 * (1.0 - res_value_rs) * af / hours_per_year
    return perhour


# Common utility code

def print_some_computational_details(mode_label, opt_gap, heur_gap, number_iterations, tot_timespent):
# Print some computational details
        print 'For mode %s the avg optimality gap is %.1f%% and the worst is %.1f%%' %  (mode_label, np.average(opt_gap), np.amax(opt_gap))
        print '         the median optimality gap is %.1f%% and the 80th percentile is %.1f%%' %  (np.percentile(opt_gap, 50), np.percentile(opt_gap, 80))
        print '     the avg heuristic gap is %.1f%% and the worst is %.1f%%' %  (np.average(heur_gap), np.amax(heur_gap))
        print '         the median heuristic gap is %.1f%% and the 80th percentile is %.1f%%' %  (np.percentile(heur_gap, 50), np.percentile(heur_gap, 80))
        print '     L-BFGS-B algorithm  iterated in avg %.1f times, using in total %.2f seconds' \
                     %  (np.average(number_iterations),tot_timespent)
  
 
def print_ratios_avg_stop_spacing(min_d_range, Tl_m): 
# print a comparison of ratios of avg stop spacing
    from math import sqrt
    num_modes = len(Tl_m)
    for m1 in range(num_modes - 1):
        Tl1 = Tl_m[m1]
        avg_d1 = np.average(min_d_range[m1])
        for m2 in range(m1 + 1, num_modes):
            Tl2 = Tl_m[m2]
            avg_d2 = np.average(min_d_range[m2])
            ratio = avg_d1/avg_d2
            ratio_approx = sqrt(Tl1/Tl2)
            error = (ratio - ratio_approx) / ratio * 100
            print 'Ratio of avd dist beween modes %2d %2d = %.2f  approx = %.2f accuracy %.0f%%' \
                        % (m1, m2, ratio, ratio_approx, error) 
                        

def print_itersection_points(y_range, min_avg_tot_cost, mode_label):
# Print points of total cost intersection between  modes
    from shapely.geometry import LineString
    num_modes = len(y_range)
    for m1 in range(num_modes - 1):
        for m2 in range(m1 + 1, num_modes):
            line1 = LineString([(a, b) for a , b in zip(y_range[m1], min_avg_tot_cost[m1])])
            line2 = LineString([(a, b) for a , b in zip(y_range[m2], min_avg_tot_cost[m2])])
            point = line1.intersection(line2)
            if not point.is_empty and point.geom_type == 'Point':
                print 'Point of intersection between modes %3s and %3s  (%5.0f, %.2f) '\
                     % (mode_label[m1], mode_label[m2], point.x, point.y)
            elif point.geom_type == 'MultiPoint':
                print 'Multiple points of intersection between modes %3s and %3s  '\
                     % (mode_label[m1], mode_label[m2])
                print '     the first is  (%5.0f, %.2f) '\
                     % (point[0].x, point[0].y)



# Common code for plots
import matplotlib.pyplot as plt
xfigsize = 6
yfigsize = 6

def plot_single(x_range, y_range, x_label, y_label, namefile, linestyle, colors, mode_label):
    fig, axes = plt.subplots(1,1, figsize=(xfigsize, yfigsize))  
    axes.set_xlabel(x_label)
    axes.set_ylabel(y_label)
    from matplotlib.ticker import MultipleLocator
    majorLocator   = MultipleLocator(10000)
    axes.xaxis.set_major_locator(majorLocator)
    #axes.grid()
    for m in range(len(x_range)):
        axes.plot(x_range[m], y_range[m], ls=linestyle[m], color=colors[m], label=mode_label[m], linewidth=3.0) 
        #axes.plot(x_range[m], y_range[m], 'o', color=colors[m], label=mode_label[m], linewidth=3.0) 
    #axes.legend(loc='upper right')
    axes.legend(loc='best',  framealpha=0.5,  prop={'size':18})
    plt.savefig(namefile)
    plt.close()	
    
    
def plot_frequency(y_range, min_f_range, f_min_range, x_label, y_label, namefile, linestyle, colors, mode_label):  
    fig, axes = plt.subplots(1,1, figsize=(xfigsize, yfigsize))  
    axes.set_xlabel(x_label)
    axes.set_ylabel(y_label)
    #axes.grid()
    for m in range(len(y_range)):
        axes.plot(y_range[m], min_f_range[m], ls=linestyle[m], color=colors[m], label=r'$\hat f$ ' + mode_label[m], linewidth=3.0) 
        m2 = m + 1
        if m2 >= len(y_range):
            m2 = 0
        axes.plot(y_range[m], f_min_range[m], ls=linestyle[m2], color=colors[m], label=r'$f_{min}$ '+ mode_label[m], linewidth=3.0)
    #axes.legend(loc='upper right')
    axes.legend(loc='best',  framealpha=0.5,  prop={'size':14})
    plt.savefig(namefile)
    plt.close()

def plot_frequency_max(y_range, min_f_range, f_max_range, x_label, y_label, namefile, linestyle, colors, mode_label):  
    fig, axes = plt.subplots(1,1, figsize=(xfigsize, yfigsize))  
    axes.set_xlabel(x_label)
    axes.set_ylabel(y_label)
    #axes.grid()
    for m in range(len(y_range)):
        axes.plot(y_range[m], min_f_range[m], ls=linestyle[m], color=colors[m], label=r'$\hat f$ ' + mode_label[m], linewidth=3.0) 
        m2 = m + 1
        if m2 >= len(y_range):
            m2 = 0
        axes.plot(y_range[m], f_max_range[m], ls=linestyle[m2], color=colors[m], label=r'$f: \theta = \theta_{min}$ '+ mode_label[m], linewidth=3.0)
    #axes.legend(loc='upper right')
    axes.legend(loc='best',  framealpha=0.5, fontsize=11)
    plt.savefig(namefile)
    plt.close()
    
    
def plot_three_economics(x_range, min_avg_op_cost, min_avg_user_cost, min_avg_tot_cost, name_outputfile0e, linestyle, colors, mode_label):
   # A figure with all plots related to the economical aspects
    fig, axes = plt.subplots(3,1) 
    fig.set_size_inches(8,12)
    from math import floor
    for i in range(3):
        axes[i].set_xlabel('Demand (pax/h)')
        #axes[i].grid(ls=':', lw=0.5)
        for item in ([axes[i].title, axes[i].xaxis.label, axes[i].yaxis.label] +
                 axes[i].get_xticklabels() + axes[i].get_yticklabels()):
            item.set_fontsize(8)
        for m in range(len(x_range)):
            if i == 0:
                axes[i].plot(x_range[m], min_avg_op_cost[m], ls=linestyle[m], lw=3, 
                                    color=colors[m], label=mode_label[m]) 
                axes[i].set_ylabel('Operator cost ($/pax)')
                axes[i].set_title('a) Operator cost')
            elif i == 1:
                axes[i].plot(x_range[m], min_avg_user_cost[m], ls=linestyle[m], lw=3, 
                                    color=colors[m], label=mode_label[m])
                axes[i].set_ylabel('User cost ($/pax)')
                axes[i].set_title('b) User cost')
            elif i == 2:
                axes[i].plot(x_range[m], min_avg_tot_cost[m], ls=linestyle[m], lw=3, 
                                    color=colors[m], label=mode_label[m])
                axes[i].set_ylabel('Total cost ($/pax)')
                axes[i].set_title('c) Total cost')   
        axes[i].legend(loc='upper right', fancybox=True, framealpha=0.5, fontsize=10)   
    plt.savefig(name_outputfile0e)
    plt.close()
