# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 16:14:38 2018

@author: Kaustav
"""

import numpy as np
from collections import deque
from GARCH_CFA_mixing import GARCH_CFA_mixing
from Monte_mixing_pw import Monte_mixing_pw
from ImpVolBrent_pw import ImpVolBrent_pw
# from DeltaStrikes_pw import DeltaStrikes_pw
import copy as cp


# Computes the price of a European put option in the GARCH diffusion model 
# with piecewise-constant parameters (rho = 0 a.e.):
    
# dS = (rd - rf)Sdt + sqrt{V}dW,
# dV = kap(the - V)dt + lam*V*dB,

# via both the Monte-Carlo and closed form approximation mixing solution methodology, 
# and compares the prices, implied vols and runtimes.


# model_params (list): [kap_ar, the_ar, lam_ar] is a list of np.arrays, 
#   where e.g. kap_ar = np.array([kap2, kap1]), piecewise parameters are 
#   specified backward in time.
# globals_params (list): [S0, V0, rd_deque, rf_deque, Strk, dt].
# S0 (float): initial spot.
# V0 (float): initial variance/volatility.
# rd_deque (deque): domestic interest rate, given backward, e.g., rd_deque = deque([rd2, rd1]).
# rf_deque (deque): foreign interest rate, given backward, e.g., rf_deque = deque([rf2, rf1]).
# Strk (float): Strike of the contract.
# dt (deque): deque of time increments over which each parameter is 'alive',
#   given backward, e.g., dt = deque([dt2, dt1]). Note sum(dt) gives option maturity T.
# N_PATH (int): # of Monte-Carlo paths.
# N_TIME (int): # of time discretisation points for SDE integration.


def GARCH_compare_pw(model_params, global_params, N_PATH, N_TIME): 
    
    
    S0 = global_params[0]
    # V0 = global_params[1]
    rd_deque = global_params[2]
    rf_deque = global_params[3]
    Strk = global_params[4]
    dt = global_params[5]
    
    # Need to append rho = 0 array into the model_params list in order for it 
    # to work with Monte_mixing_pw
    rho0append = np.zeros([len(rd_deque)])
    model_params_monte = cp.deepcopy(model_params)
    model_params_monte.append(rho0append)
    
    
    # Compute prices
    (GPrice_approx, elapsed_approx) = GARCH_CFA_mixing(model_params, global_params)
    (GPrice_Monte, __, elapsed_MC)  = Monte_mixing_pw(model_params_monte, global_params, N_PATH, N_TIME, 'g', 'Put')
    
    # Compute implied volatilities  
    ImpVolG_Monte = ImpVolBrent_pw(S0, Strk, rd_deque, rf_deque, dt, 'Put', GPrice_Monte)
    ImpVolG_approx = ImpVolBrent_pw(S0, Strk, rd_deque, rf_deque, dt, 'Put', GPrice_approx) 
    
    # Compute implied volatilities in cent
    ImpVolG_approx_cent = 100*ImpVolG_approx
    ImpVolG_Monte_cent = 100*ImpVolG_Monte
    
    # Absolute and relative error of option price
    Abs_Err = np.abs(GPrice_approx - GPrice_Monte)
    Abs_Err_cent = 100*Abs_Err
    Rel_Err_cent = 100*Abs_Err/GPrice_Monte
        
    # Absolute and relative error of Implied vols
    ImpVol_AbsErr = np.abs(ImpVolG_approx - ImpVolG_Monte)
    ImpVol_AbsErr_bp = 10000*ImpVol_AbsErr
    ImpVol_RelErr_cent = 100*ImpVol_AbsErr/ImpVolG_Monte if ImpVolG_Monte != 0 else 0
    
    # speedup
    speedup = elapsed_MC/elapsed_approx
    
    
    
    # Print some things!
    print("\nModel", "GARCH diffusion", "\n")
    
    print("Strike", "{:.7f}".format(Strk), "\n")
              
    print("Put prices")
    print("MC", "{:.7f}".format(GPrice_Monte), "CF", "{:.7f}".format(GPrice_approx), \
          "AbsErr","{:.3f}".format(Abs_Err_cent)+"%", "RelErr","{:.3f}".format(Rel_Err_cent)+"%", "\n")
    
    print("Implied Vols")
    print("MC", "{:.3f}".format(ImpVolG_Monte_cent)+"%", "CF", "{:.3f}".format(ImpVolG_approx_cent)+"%", \
          "AbsErr","{:.3f}".format(ImpVol_AbsErr_bp)+"bp", "RelErr","{:.3f}".format(ImpVol_RelErr_cent)+"%", "\n")
        
    print("Elapsed")
    print("MC", "{:.3f}".format(elapsed_MC)+"s", "CF", "{:.5f}".format(elapsed_approx)+"s", \
          "Speed up","{:.2f}".format(speedup)+"x")

    
    return

    







# Example. Uncomment the following code then run


# Common parameters
S0 = 100
V0 = 0.0036       # remember V0 is initial variance for Heston and GARCH
rd3 = 0.01
rd2 = 0.01
rd1 = 0.02
rf3 = 0.00
rf2 = 0.00
rf1 = 0.00
dt3 = 6/12
dt2 = 3/12
dt1 = 3/12

rd_deque = deque([rd3, rd2, rd1])
rf_deque = deque([rf3, rf2, rf1])
dt = deque([dt3, dt2, dt1])

# Delta = 0.5
# sig = np.sqrt(V0)
# Strk = DeltaStrikes_pw(S0, sig, rd_deque, rf_deque, dt, Delta, 'Put')
Strk = S0*1.01

# # Model parameters
kap3 = 5.0
kap2 = 5.1
kap1 = 4.9
the3 = 0.019
the2 = 0.020
the1 = 0.018
lam3 = 0.4142
lam2 = 0.4342
lam1 = 0.3942

kap_ar = np.array([kap3, kap2, kap1])
the_ar = np.array([the3, the2, the1])
lam_ar = np.array([lam3, lam2, lam1])

# Monte-Carlo parameters
T = sum(dt)
N_PATH = 100000
Steps_day_N = 24
N_TIME = int(Steps_day_N*253*T)

global_params = [S0, V0, rd_deque, rf_deque, Strk, dt]
model_params = [kap_ar, the_ar, lam_ar]


GARCH_compare_pw(model_params, global_params, N_PATH, N_TIME)
    