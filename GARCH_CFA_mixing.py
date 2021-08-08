# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 15:34:44 2017

@author: Dr Kaustav Das (kaustav.das@monash.edu)
"""


import numpy as np
import copy as cp
from collections import deque
from Omg import Omg1, Omg2, Omg3, Omg4, Omg5
from parPutBS_pw import parPutBS_pw
# from DeltaStrikes_pw import DeltaStrikes_pw
from timeit import default_timer as timer 



# Computes the price of a European put option in the GARCH diffusion model 
# with piecewise-constant parameters (rho = 0 a.e.):

# dS = (rd - rf)Sdt + sqrt{V}dW,
# dV = kap(the - V)dt + lam*V*dB, 

# via the closed form approximation mixing solution methodology.


# model_params (list): [kap, the, lam] is a list of np.arrays, 
#   where e.g. kap = np.array([kap2, kap1]), piecewise parameters are 
#   specified backward in time.
# globals_params (list): [S0, V0, rd_deque, rf_deque, Strk, dt].
# S0 (float): initial spot.
# V0 (float): initial variance/volatility.
# rd_deque (deque): domestic interest rate, given backwards , e.g., rd_deque = deque([rd2, rd1]).
# rf_deque (deque): foreign interest rate, given backwards, e.g., rf_deque = deque([rf2, rf1]).
# Strk (float): Strike of the contract.
# dt (deque): deque of time increments over which each parameter is 'alive',
#   given backwards, e.g., dt = deque([dt2, dt1]). Note sum(dt) gives option maturity T.


def GARCH_CFA_mixing(model_params, global_params):
    
    
    # These are actually all arrays, but we have not put the subscript 'ar'
    # for simplicity
    kap = model_params[0]
    the = model_params[1]
    lam = model_params[2]

    S0 = global_params[0]
    V0 = global_params[1]
    _rd_deque = global_params[2]
    _rf_deque = global_params[3]
    Strk = global_params[4]
    _dt = global_params[5]
    
    # Start timer
    start = timer()
    
    # Copy deques
    rd_deque = cp.copy(_rd_deque)
    rf_deque = cp.copy(_rf_deque)
    dt = cp.copy(_dt)
    
    # Initialise deques for the SDE parameters
    d_kap = deque([])
    d_mkap = deque([])
    d_kapthe = deque([])
    d_lamsq = deque([])
    d_mlamsqplkap = deque([])
    d_ones = deque([])
    
    # Calculate inputs for Omg functions  
    mkap = -1*kap
    kapthe = kap*the
    lamsq = lam**2
    mlamsqplkap = -1*lamsq+kap
    
    # Convert SDE np.arrays to deques
    for i in range(0, len(kap)):
        d_kap.append(kap[i])
        d_mkap.append(mkap[i])
        d_kapthe.append(kapthe[i])
        d_lamsq.append(lamsq[i])
        d_mlamsqplkap.append(mlamsqplkap[i])
        d_ones.append(1)
        
    # Compute point y and PutBS partial derivatives    
    y = V0*Omg1(d_mkap, d_ones, dt) + Omg2(d_mkap, d_kap, d_ones, d_kapthe, dt)    
    (PutBS, xxPutBS, yyPutBS, xyPutBS) = parPutBS_pw(S0, y, Strk, rd_deque, rf_deque, dt)
    
    bigterm = ((V0**2)*Omg3(d_mkap, d_mkap, d_lamsq, d_ones, d_ones, d_lamsq, dt)+\
    2*V0*Omg4(d_mkap, d_mkap, d_lamsq, d_mlamsqplkap, d_ones, d_ones, d_lamsq, d_kapthe, dt)+\
    2*Omg5(d_mkap, d_mkap, d_lamsq, d_mlamsqplkap, d_kap, d_ones, d_ones, d_lamsq, d_kapthe, d_kapthe, dt))
                    
    yyterm = yyPutBS*bigterm
    
    price = PutBS + yyterm
    
    end = timer()
    
    elapsed = end - start

    return (price, elapsed)
    





# Example code.

if __name__ == '__main__':


    # Common parameters
    S0 = 100
    V0 = 0.0036       # remember V0 is initial variance for Heston and GARCH
    rd3 = 0.02
    rd2 = 0.03
    rd1 = 0.01
    rf3 = 0.00
    rf2 = 0.00
    rf1 = 0.00
    dt3 = 6/12
    dt2 = 3/12
    dt1 = 3/12
    
    rd_deque = deque([rd3, rd2, rd1])
    rf_deque = deque([rf3, rf2, rd1])
    dt = deque([dt3, dt2, dt1])
    
    # Delta = 0.5
    # sig = np.sqrt(V0)
    # Strk = DeltaStrikes_pw(S0, sig, rd_deque, rf_deque, dt, Delta, 'Put')
    Strk = S0*1.01
    
    # Model parameters
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
    
    global_params = [S0, V0, rd_deque, rf_deque, Strk, dt]
    model_params = [kap_ar, the_ar, lam_ar]
    
    
    print(GARCH_CFA_mixing(model_params, global_params))