# Model IV
# Version 23th July 2015
# The python code solves Model IV for the optimization of a transit line 
# and for the strategic assessment of mode choice, see: 
# L. Moccia, and G. Laporte, 'Improved models for technology choice in a transit corridor with fixed demand', 2015
# For comments, please write to: 
# moccia@icar.cnr.it
# Extension of the Model from Tirachini et al 2010 published in TRE
# Optimal stop spacing considering vehicle dynamics
# Variable train consist
# Two periods, peak and off-peak
# Demand inelastic
# Minimization of user and agency cost,  
# Variables: fp, fo, d, npeak, noff 



# Import formulae, and code shared by all models
from common_code import *

# Import parameters shared by all models
from common_parameters import *


# Computed cost parameters for Model IV
hours_per_year  = 5500.0

c0l_m = [infra_fixed_cost_per_hour(m, hours_per_year, L, width_m, land_cost_per_hectare, \
        infracost_m, infralife_m, inframaint_m, discount_rate) \
        for m in range(num_modes)]

c0s_m = [stop_fixed_cost_per_hour(m, hours_per_year, stopcost_m, infralife_m, discount_rate) \
         for m in range(num_modes)]

c1v_m = [rolling_stock_cost_per_hour(m, hours_per_year, vehcost_m, vehlife_m, discount_rate, \
            res_value_rs) \
        for m in range(num_modes)]


R_m     = [2 * L / Vmax_m[m] for m in range(num_modes)] # Running time at speed Vmax (h)


# Formulae for Model IV

def Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, fp, fo, l, L, betav, R, Tl, npeak, noff, kv):
# User cost
    yo = y * gamma
    Ca = Pa * d * y / (2.0 * v)  * (chip  +  gamma * chio)         # Cost of access
    Cw = chip * Cost_waiting(fp, y, t01, t02, t11, t12, fbar, Pw, epsilon) \
       + chio * Cost_waiting(fo, yo, t01, t02, t11, t12, fbar, Pw, epsilon)  # Cost of waiting
    deltap  = penalty_crowding(fp, npeak, y, l, L, kv, xi, rho, thetamin)
    deltao  = penalty_crowding(fo, noff, yo, l, L, kv, xi, rho, thetamin)
    ctp     = cycle_time(fp, npeak, d, y, betav, R, L, Tl)
    cto     = cycle_time(fo, noff, d, yo, betav, R, L, Tl)
    Cv = Pv * l * y / (2.0 * L) * (chip * deltap * ctp  +  chio * gamma * deltao * cto) # Cost of in-vehicle time
    return Ca + Cw + Cv


def Co(c0l, c0s,c1t, c1v, y, betav, R, fp, fo, c2v, L, d, Tl, npeak, noff):
# Agency cost
    yo = gamma * y
    ctp = cycle_time(fp, npeak, d,  y, betav, R, L, Tl)
    cto = cycle_time(fo,  noff, d, yo, betav, R, L, Tl)
    dim_fleet_p = npeak * fp *  ctp
    dim_fleet_o = noff  * fo *  cto
    # check if fleet peak is larger tha off-peak
    if dim_fleet_p >= dim_fleet_o:
        largest_fleet = dim_fleet_p
    else:
        largest_fleet = dim_fleet_o
    cost_fixed      =   c0l  +  2.0 * c0s * L /d 
    cost_fixed_veh  =   c1v * largest_fleet 
    cost_crew_p     =   c1t * chip * fp *  ctp  
    cost_crew_o     =   c1t * chio * fo *  cto  
    cost_veh_km     =   2.0 * c2v * L * (chip * npeak * fp   +  chio * noff * fo)
    cost = cost_fixed + cost_fixed_veh + cost_crew_p + cost_crew_o + cost_veh_km
    return cost

def Ctot(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, fp, fo, l, L, betav, R, c0l, c0s,c1t, c1v, c2v, Tl, npeak, noff, kv):
# Total cost sum of user and agency cost
    return Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, fp, fo, l, L, betav, R, Tl, npeak, noff, kv) + \
            Co(c0l, c0s,c1t, c1v, y, betav, R, fp, fo, c2v, L, d, Tl, npeak, noff)


def coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fp, fo, l, L, betav, R, c0l, c0s,c1t, c1v, c2v, Tl,  kv):
# Coefficients a0, a1,...,a11
    if (fp >= fbar):
        a0 = Pv * chip *  l * xi/ (2.0 * L) * R * y + c0l + c1v * betav/3600.0 * y   
        a4 = Pw * chip *  epsilon * y  
    else:
        a0 = Pw * chip *  tw / 60.0 * y + Pv * chip *  l * xi/ (2.0 * L) * R * y + c0l + c1v * betav/3600.0 * y  
        a4 = Pw * chip *  mu * epsilon * y 
    if (fo >= fbar):
        a0 += Pv * chio *  l * xi/ (2.0 * L) * R * y * gamma
        a13 = Pw * chio *  epsilon * y * gamma
    else:
        a0 += Pw * chio *  tw / 60.0 * y * gamma
        a13 = Pw * chio *  mu * epsilon * y * gamma
    a1  = Pa  *  y / (2.0 * v) * (chip + gamma * chio)
    a2  = Pv * chip *  l * Tl * xi *  y  + 2.0 * c0s * L + Pv * chio *  l * Tl * xi *  y * gamma
    a3  = c1t * chip *  R  
    a5  = 2 * c1t * chip *  L * Tl
    a6  = c1t * chip *  betav/3600.0 * y
    a7  = c1v * R   +  2.0 * c2v * chip *  L
    a8  = Pv * chip * l * y *  (xi * betav * y /(7200.0 * L)   +  rho * R * l * y/ (4.0 * L * L * kv) )
    a9  = 2 * c1v * L * Tl
    a10 = Pv * chip *  rho * l * l * y * y * Tl /(2.0 * L * kv)
    a11 = Pv * chip *  rho * l * l * y * y * y * betav / (14400.0 * L * L * kv)
    
    a12 = c1t * chio *  R
    a14 = 2 * c1t * chio *  L * Tl
    a15 = c1t * chio *  betav/3600.0 * y * gamma
    a16 = 2.0 * c2v * chio *  L
    a17 = Pv * chio * l * y * gamma * (xi * betav * y * gamma/(7200.0 * L)   +  rho * R * l * y * gamma/ (4.0 * L * L * kv) )
    a18 = Pv * chio *  rho * l * l * y * y * gamma * gamma * Tl /(2.0 * L * kv)
    a19 = Pv * chio *  rho * l * l * y * y * y * gamma * gamma * gamma * betav / (14400.0 * L * L * kv)
    
    return a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19
    
    
def solve_heuristicaly(fminp, fmaxp, fmino, fmaxo, dmin, dmax):
# return heuristic solution of the lower convex approximation of the objective function
    def compute_my_f_d(myfminp, myfmaxp, myfmino, myfmaxo):
    # Compute the constrained approximated f and d in the range myfmin, myfmax
        from math import sqrt
        num = a2 + myfminp * (a5 + npeak * a9) + a10/float(npeak * myfmaxp) \
                 + myfmino * a14               + a18/float(noff  * myfmaxo)
        den = float(a1)
        myd = sqrt(num/den)
        if myd < dmin:
            myd = dmin
        elif myd > dmax:
            myd = dmax
        num = a4 + a8/float(npeak) + a11/float(npeak * npeak * myfmaxp)
        den = float(a3  +  npeak * a7)
        if den > 0:
            myfp = sqrt(num/den)
        else:
            myfp = myfminp
        if myfp < myfminp:
            myfp = myfminp
        elif myfp > myfmaxp:
            myfp = myfmaxp
        
        num = a13 + a17/float(noff) + a19/float(noff * noff * myfmaxo)
        den = float(a12  +  noff * a16)
        if den > 0:
            myfo = sqrt(num/den)
        else:
            myfo = myfmino
        if myfo < myfmino:
            myfo = myfmino
        elif myfo > myfmaxo:
            myfo = myfmaxo
        
        myv = (a0 + a1 * myd + a2/float(myd) + a3 * myfp + a4/float(myfp) \
                + a5 * myfminp/float(myd) + a6 /float(npeak) + a7 * npeak * myfp + a8/float(npeak * myfp) \
                + a9 * npeak * myfminp/myd + a10 /float(npeak * myfmaxp * myd) + a11/float(npeak * npeak * myfmaxp * myfp) \
                + a12 * myfo + a13/float(myfo) + a14 * myfmino/float(myd) + a15/float(noff) \
                + a16 *noff * myfo + a17/float(noff * myfo) + a18/float(noff * myfmaxo * myd) \
                + a19/float(noff * noff * myfmaxo * myfo) \
                ) / y
        return myv, myfp, myfo, myd
    
    
    if (fmino >= fbar and fminp >= fbar):
    # Case 1, only case b is relevant
        a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fminp, fmino, l, L, betav, R, \
                                c0l, c0s,c1t, c1v, c2v, Tl, kv)
        value, fp_approx, fo_approx, d_approx = compute_my_f_d(fminp, fmaxp, fmino, fmaxo)
    if (fmino < fbar and fminp >= fbar):
    # Case a could be relevant for fo
        #Case a
        fsupoa = min(fbar, fmaxo)
        a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fminp, fmino, l, L, betav, R, \
                                c0l, c0s,c1t, c1v, c2v, Tl, kv)
        value, fp_approx, fo_approx, d_approx = compute_my_f_d(fminp, fmaxp, fmino, fsupoa)
        # Case b, only to be evaluated if fmaxo >= fbar
        if (fmaxo >= fbar):
            a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19 = \
                    coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fminp, fbar, l, L, betav, R, \
                                    c0l, c0s,c1t, c1v, c2v, Tl, kv)
            valueb, fp_approxb, fo_approxb, d_approxb = compute_my_f_d(fminp, fmaxp, fbar, fmaxo)
            if (valueb < value):
                value       = valueb
                fp_approx   = fp_approxb
                fo_approx   = fo_approxb
                d_approx    = d_approxb
    if (fmino >= fbar and fminp < fbar):
    # Case a could be relevant for fp. Observe that this case wouldn't occur if the constraint fp >= fo is enforced,
    #   but, we consider it here to allow test with demand larger in off-peak hours -- Test only
        #Case a
        fsuppa = min(fbar, fmaxp)
        a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fminp, fmino, l, L, betav, R, \
                                c0l, c0s,c1t, c1v, c2v, Tl, kv)
        value, fp_approx, fo_approx, d_approx = compute_my_f_d(fminp, fsuppa, fmino, fmaxo)
        # Case b, only to be evaluated if fmaxp >= fbar
        if (fmaxp >= fbar):
            a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19 = \
                    coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fbar, fmino, l, L, betav, R, \
                                    c0l, c0s,c1t, c1v, c2v, Tl, kv)
            valueb, fp_approxb, fo_approxb, d_approxb = compute_my_f_d(fbar, fmaxp, fmino, fmaxo)
            if (valueb < value):
                value       = valueb
                fp_approx   = fp_approxb
                fo_approx   = fo_approxb
                d_approx    = d_approxb
    if (fmino < fbar and fminp < fbar):
    # Case a could be relevant for both fp and fo
    # There are 4 cases
        #Case 1, a relevant for fo only
        fsupoa = min(fbar, fmaxo)
        new_minfp = min(fbar, fmaxp)
        a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, new_minfp, fmino, l, L, betav, R, \
                                c0l, c0s,c1t, c1v, c2v, Tl, kv)
        value1, fp_approx1, fo_approx1, d_approx1 = compute_my_f_d(new_minfp, fmaxp, fmino, fsupoa)  
        #Case 2, a relevant for fp only
        fsuppa = min(fbar, fmaxp)
        new_minfo = min(fbar, fmaxo)
        a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fminp, new_minfo, l, L, betav, R, \
                                c0l, c0s,c1t, c1v, c2v, Tl, kv)
        value2, fp_approx2, fo_approx2, d_approx2 = compute_my_f_d(fminp, fsuppa, new_minfo, fmaxo)  
        #Case 3, a relevant for both fp and fo
        fsuppa = min(fbar, fmaxp)
        fsupoa = min(fbar, fmaxo)
        a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fminp, fmino, l, L, betav, R, \
                                c0l, c0s,c1t, c1v, c2v, Tl, kv)
        value3, fp_approx3, fo_approx3, d_approx3 = compute_my_f_d(fminp, fsuppa, fmino, fsupoa)  
        #Case 4, a not relevant for both
        if (fmaxp >= fbar and fmaxo >= fbar):
            a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19 = \
                coefficients(Pa, Pw, Pv, v, y, tw, mu, epsilon, fbar, fbar, l, L, betav, R, \
                                c0l, c0s,c1t, c1v, c2v, Tl, kv)
            value4, fp_approx4, fo_approx4, d_approx4 = compute_my_f_d(fbar, fmaxp, fbar, fmaxo)
        else:
            value4 = 2 * value3  #a shortcut to eliminate it from comparison
        # Take the best
        if (value1 == min(value1, value2, value3, value4)):
            value       = value1
            fp_approx   = fp_approx1
            fo_approx   = fo_approx1
            d_approx    = d_approx1
        elif (value2 == min(value1, value2, value3, value4)):
            value       = value2
            fp_approx   = fp_approx2
            fo_approx   = fo_approx2
            d_approx    = d_approx2
        elif (value3 == min(value1, value2, value3, value4)):
            value       = value3
            fp_approx   = fp_approx3
            fo_approx   = fo_approx3
            d_approx    = d_approx3
        elif (value4 == min(value1, value2, value3, value4)):
            value       = value4
            fp_approx   = fp_approx4
            fo_approx   = fo_approx4
            d_approx    = d_approx4
    
    return value, fp_approx, fo_approx, d_approx




def Cop_of_f_and_d(vars):
# Operator cost as a function of f and d
    fp, fo, d = vars
    return Co(c0l, c0s,c1t, c1v, y, betav, R, fp, fo, c2v, L, d, Tl, npeak, noff)
    
def Cuser_of_f_and_d(vars):
# user cost as a function of f and d
    fp, fo, d = vars
    return Cu(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, fp, fo, l, L, betav, R, Tl, npeak, noff, kv)

def Ctot_of_f_and_d(vars):
# Total cost as a function of f and d
    fp, fo, d = vars
    return Ctot(Pa, Pw, Pv, d, v, y, tw, mu, epsilon, fp, fo, l, L, betav, R, c0l, c0s,c1t, c1v, c2v, Tl, npeak, noff, kv)

def unpack_mode_specific_par(mode):
# return mode specific parameters
    return  c0l_m[mode], c0s_m[mode], c1t_m[mode], c1v_m[mode], c2v_m[mode], betav_m[mode], R_m[mode], Tl_m[m], \
            nmin_m[m], nmax_m[m], kv_m[m]


    

# SOLVE
import numpy as np
from scipy.optimize import minimize
max_demand          = np.empty(num_modes)
y_range             = [np.array([]) for m in range(num_modes)]
min_fp_range        = [np.array([]) for m in range(num_modes)]
min_fo_range        = [np.array([]) for m in range(num_modes)]
fp_min_range        = [np.array([]) for m in range(num_modes)]
fo_min_range        = [np.array([]) for m in range(num_modes)]
min_d_range         = [np.array([]) for m in range(num_modes)]
min_avg_op_cost     = [np.array([]) for m in range(num_modes)]
min_avg_user_cost   = [np.array([]) for m in range(num_modes)]
min_avg_tot_cost    = [np.array([]) for m in range(num_modes)]
commercial_speedp   = [np.array([]) for m in range(num_modes)]
avg_running_speed   =  [np.array([]) for m in range(num_modes)]
vehicle_per_TUp     = [np.array([]) for m in range(num_modes)]
commercial_speedo   = [np.array([]) for m in range(num_modes)]
vehicle_per_TUo     = [np.array([]) for m in range(num_modes)]
opt_gap             = [np.array([]) for m in range(num_modes)]
heur_gap            = [np.array([]) for m in range(num_modes)]
number_iterations   = [np.array([]) for m in range(num_modes)]
fp_heur_range       = [np.array([]) for m in range(num_modes)]
fo_heur_range       = [np.array([]) for m in range(num_modes)]
fp_max_range        = [np.array([]) for m in range(num_modes)]
fo_max_range        = [np.array([]) for m in range(num_modes)]

for m in range(num_modes):
    max_demand[m]       = kv_m[m] * nmax_m[m] * fmax_m[m] * nu / alpha
    if max_demand[m] > max_demand_studied:
        max_demand[m] = max_demand_studied
    y_range[m]          = np.arange(min_demand_studied, max_demand[m], step=step_demand)
    min_fp_range[m]     = np.empty_like(y_range[m])
    min_fo_range[m]     = np.empty_like(y_range[m])
    fp_min_range[m]     = np.empty_like(y_range[m])
    fo_min_range[m]     = np.empty_like(y_range[m])
    min_d_range[m]      = np.empty_like(y_range[m])
    min_avg_tot_cost[m] = np.empty_like(y_range[m])
    min_avg_op_cost[m]  = np.empty_like(y_range[m])
    min_avg_user_cost[m]= np.empty_like(y_range[m])
    commercial_speedp[m]= np.empty_like(y_range[m])
    vehicle_per_TUp[m]	= np.empty_like(y_range[m])
    commercial_speedo[m]= np.empty_like(y_range[m])
    avg_running_speed[m]= np.empty_like(y_range[m])
    vehicle_per_TUo[m]	= np.empty_like(y_range[m])
    opt_gap[m]	        = np.empty_like(y_range[m])
    heur_gap[m]	        = np.empty_like(y_range[m])
    number_iterations[m]= np.empty_like(y_range[m])
    fp_heur_range[m]    = np.empty_like(y_range[m])
    fo_heur_range[m]    = np.empty_like(y_range[m])
    fp_max_range[m]     = np.empty_like(y_range[m])
    fo_max_range[m]     = np.empty_like(y_range[m])
    fmax                = fmax_m[m]
    dmin                = Dmin_m[m]
    dmax                = 2.0
    c0l, c0s,c1t, c1v, c2v, betav, R, Tl, nmin, nmax, kv = unpack_mode_specific_par(m)
    import time
    tot_timespent       = 0.0
    count_larger_fleet  = 0
    count_num_sol       = 0
    for i, y in enumerate(y_range[m]):
        best_value = 10000.0
        if i == 0: 
            nmin_p = nmin
            nmin_o = nmin
        else:
            nmin_p = int(vehicle_per_TUp[m][i - 1])
            nmin_o = int(vehicle_per_TUo[m][i - 1])
        
        for npeak in range(nmin_p, nmax + 1):
            for noff in range(nmin_o, npeak + 1):  
                fminp   = alpha * y /(nu * kv * float(npeak)) # Minimal frequency peak
                fmino   = alpha * y * gamma /(nu * kv * float(noff)) # Minimal frequency off-peak
                
                #Compute only if these value of npeak and noff result in a feasible range for freq 
                if fmax >= fminp and fmax >= fmino:  
                    max_ftp = max_freq_theta(npeak, y, l, L, kv, thetamin)
                    max_fto = max_freq_theta(noff, y * gamma, l, L, kv, thetamin)
                    if max_ftp > fmax:
                        max_ftp = fmax
                    if max_fto > fmax:
                        max_fto = fmax   
                    value_lower, fp_heur, fo_heur, d_heur = solve_heuristicaly(fminp, max_ftp, fmino, max_fto, dmin, dmax)
                    value_heur              =  Ctot_of_f_and_d((fp_heur, fo_heur, d_heur))  / y
                    bnds 					= ((fminp, fmax), (fmino, fmax),(dmin, dmax))
                    #Is the fleet at peak larger than at off-peak? Approximate condition
                    if fp_heur* npeak > fo_heur * noff:
                        fleet_p_larger = True
                    else:
                        fleet_p_larger = False
                    #starting point
                    fp_start = fp_heur
                    fo_start = fo_heur
                    d_start  = d_heur
                    sol_start = (fp_start, fo_start, d_start)
                    x0 	     = np.array(sol_start)
                    timestart = time.clock()
                    res = minimize(Ctot_of_f_and_d, x0=x0, bounds=bnds, method='L-BFGS-B')
                    bestfp = res.x[0]
                    bestfo = res.x[1]
                    
                    sol = (bestfp, bestfo, res.x[2])
                    tot_timespent += time.clock() - timestart
                    #sol = sol_start
                    value = Ctot_of_f_and_d(sol) / y  
                    if value < best_value:
                        min_fp_range[m][i]      = bestfp
                        min_fo_range[m][i]      = bestfo
                        min_d_range[m][i]       = res.x[2]
                        fp_min_range[m][i]      = fminp
                        fo_min_range[m][i]      = fmino
                        #print m, fminp, min_fp_range[m][i], ftp, fmino, min_fo_range[m][i], fto, value, y
                        min_avg_tot_cost[m][i]  = value
                        min_avg_op_cost[m][i]   = Cop_of_f_and_d(sol)   / y
                        min_avg_user_cost[m][i] = Cuser_of_f_and_d(sol) / y
                        #print '     ', min_avg_op_cost[m][i], min_avg_user_cost[m][i]
                        commercial_speedp[m][i]  = 2.0 * L / (cycle_time(bestfp, npeak, res.x[2], y, betav, R, L, Tl))
                        commercial_speedo[m][i]  = 2.0 * L / (cycle_time(bestfo, npeak, res.x[2], gamma * y, betav, R, L, Tl))
                        avg_running_speed[m][i] = 2.0 * L / (R  + 2.0 * L * Tl / res.x[2])
                        vehicle_per_TUp[m][i]    = npeak
                        vehicle_per_TUo[m][i]    = noff
                        opt_gap[m][i]           = (value - value_lower) / value * 100
                        heur_gap[m][i]          = (value_heur - value) / value * 100
                        fp_heur_range[m][i]     = fp_heur
                        fo_heur_range[m][i]     = fo_heur
                        fp_max_range[m][i]      = max_freq_theta(npeak, y, l, L, kv, thetamin)
                        fo_max_range[m][i]      = max_freq_theta(noff, y * gamma, l, L, kv, thetamin)
                        number_iterations[m][i] = res.nit 
                        best_value = value
                        count_num_sol += 1
                        if fleet_p_larger:
                            count_larger_fleet += 1
            # if opt_gap[m][i] > 3.0:
#                 print m,  y, opt_gap[m][i],  min_fp_range[m][i], fp_heur_range[m][i], min_fo_range[m][i], fo_heur_range[m][i], vehicle_per_TUp[m][i], vehicle_per_TUo[m][i]
        #print m, y, min_fp_range[m][i], fp_max_range[m][i]
    print 'The following optimality gaps are exact only when the unconstrained solution belongs to the domain'
    print 'Non-linear constraint approximately satisfied (%): ', count_larger_fleet/float(count_num_sol)*100
    print_some_computational_details(mode_label[m], opt_gap[m], heur_gap[m], number_iterations[m], tot_timespent/nmax)
    #print m, min_fp_range[m]



# Print a comparison of ratios of avg stop spacing
print_ratios_avg_stop_spacing(min_d_range, Tl_m)



# Print points of total cost intersection between  modes
print_itersection_points(y_range, min_avg_tot_cost, mode_label)
		

# PLOT
import matplotlib.pyplot as plt

model_id = 'fig_M4_'


name_outputfile1 = model_id + 'freq_peak.pdf'
plot_frequency(y_range, min_fp_range, fp_min_range, 'Demand (pax/h)', 'Peak frequency (TU/h)', name_outputfile1, linestyle, colors, mode_label)

name_outputfile1max = model_id + 'freq_p_max.pdf'
plot_frequency_max(y_range, min_fp_range, fp_max_range, 'Demand (pax/h)', 'Peak frequency (TU/h)', name_outputfile1max, linestyle, colors, mode_label)

name_outputfile1maxo = model_id + 'freq_o_max.pdf'
plot_frequency_max(y_range, min_fo_range, fo_max_range, 'Demand (pax/h)', 'Off-peak frequency (TU/h)', name_outputfile1maxo, linestyle, colors, mode_label)


name_outputfile1off = model_id + 'freq_off.pdf'
plot_frequency(y_range, min_fo_range, fo_min_range, 'Demand (pax/h)', 'Off-peak frequency (TU/h)', name_outputfile1off, linestyle, colors, mode_label)


name_outputfile2 = model_id + 'tot_cost.pdf'
plot_single(y_range, min_avg_tot_cost, 'Demand (pax/h)', 'Average total cost ($/pax)', name_outputfile2, linestyle, colors, mode_label)

name_outputfile3 = model_id + 'op_cost.pdf'
plot_single(y_range, min_avg_op_cost, 'Demand (pax/h)', 'Operator cost ($/pax)', name_outputfile3, linestyle, colors, mode_label)

name_outputfile4 = model_id + 'user_cost.pdf'
plot_single(y_range, min_avg_user_cost, 'Demand (pax/h)', 'User cost ($/pax)', name_outputfile4, linestyle, colors, mode_label)


# A figure with all plots related to the economical aspects
name_outputfile0e = model_id + 'all_eco_plots.pdf'
plot_three_economics(y_range, min_avg_op_cost, min_avg_user_cost, min_avg_tot_cost, name_outputfile0e, linestyle, colors, mode_label)


name_outputfile5 = model_id + 'stop_dist.pdf'
plot_single(y_range, min_d_range, 'Demand (pax/h)', 'Stop spacing (km)', name_outputfile5, linestyle, colors, mode_label)


name_outputfile6 = model_id + 'commercial_speed.pdf'
fig, axes = plt.subplots(1,1, figsize=(xfigsize, yfigsize))   
axes.set_xlabel('Demand (pax/h)')
axes.set_ylabel('Commercial speed (km/h)')
#axes.grid()
for m in range(num_modes):
    axes.plot(y_range[m], commercial_speedp[m], ls='--', lw=3,color=colors[m], label='peak ' + mode_label[m])
    axes.plot(y_range[m], commercial_speedo[m], ls=':', lw=3,color=colors[m], label='off-peak ' + mode_label[m])
axes.legend(loc='best', fancybox=True, framealpha=0.5, fontsize=11)
plt.savefig(name_outputfile6)
plt.close()

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
        axes.plot(y_range[m], vehicle_per_TUp[m], ls=linestyle[m], lw=3,color=colors[m], label=r'$n^p$ ' +mode_label[m])
        m2 = m - 2
        if m2 < 0: 
            m2 = 0
        axes.plot(y_range[m], vehicle_per_TUo[m], ls=linestyle[m2], lw=3,color=colors[m], label=r'$n^o$ ' + mode_label[m])
axes.legend(loc='best', fancybox=True, framealpha=0.5, fontsize=11)
plt.savefig(name_outputfile8)
plt.close()



# With MAC open the pdf output files
import platform
if platform.system() == 'Darwin':
    import os
    command = 'open -a Preview.app ' \
        +        name_outputfile1  + ' ' +  name_outputfile1max  + ' ' + name_outputfile2 \
        + ' ' +  name_outputfile3  + ' ' +  name_outputfile4    + ' ' + name_outputfile1maxo\
        + ' ' +  name_outputfile0e + ' ' +  name_outputfile5 \
        + ' ' +  name_outputfile6  + ' ' +  name_outputfile7 \
        + ' ' + name_outputfile1off + ' ' +  name_outputfile8
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


# Save on file some results  - TEST only
# from numpy import savetxt
# for m in range(num_modes):
#     name_test_file = 'mod3bistest' + str(m) + '.txt'
#     out_file = open(name_test_file,"w")
#     savetxt(out_file, np.c_[heur_gap[m], min_fp_range[m], f_heur_range[m], fp_min_range[m]], fmt='%.2f  %4d %4d  %4d')
#     out_file.close()