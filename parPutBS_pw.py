# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 12:58:18 2017

@author: kaustavd
"""


from math import sqrt, exp, log
import copy as cp
from collections import deque
from scipy.stats import norm


# Computes the PutBS function as well as second order Greeks.

# x and y are the first two arguments of PutBS.
# Strk (float): Strike of the contract.
# rd_deque (deque): domestic interest rate, given backward, e.g., rd_deque = deque([rd2, rd1]).
# rf_deque (deque): foreign interest rate, given backward, e.g., rf_deque = deque([rf2, rf1]).
# dt (deque): deque of time increments over which each parameter is 'alive',
#   given backward, e.g., dt = deque([dt2, dt1]). Note sum(dt) gives option maturity T.


def parPutBS_pw(x, y, Strk, _rd_deque, _rf_deque, _dt):
    
    
    # Copy deques
    rd_deque = cp.copy(_rd_deque)
    rf_deque = cp.copy(_rf_deque)
    dt = cp.copy(_dt)
    
    sqrty = sqrt(y)
    
    # We now compute discretised versions of int (rd - rf)dt, e^(-int rd dt)
    # and e^(-int rf dt), as well as T.
    rsumdt = 0
    expmrd = 1
    expmrf = 1
    T = 0
    
    lastlayer = deque([])
    while dt != lastlayer:    
        DT = dt.popleft()
        RD = rd_deque.popleft()
        RF = rf_deque.popleft()
        R = RD - RF
        rsumdt += R*DT
        expmrd *= exp(-DT*RD)
        expmrf *= exp(-DT*RF)
        T += DT
        
        
        
        
    dpl = (log(x/Strk) + rsumdt)/sqrty + 0.5*sqrty
    dm = dpl - sqrty
    
    expmrf_phidpl = expmrf*norm.pdf(dpl)
    
    PutBS  = Strk*expmrd*norm.cdf(-1.0*dm) - x*expmrf*norm.cdf(-1.0*dpl)
    xxPutBS = expmrf_phidpl/(x*sqrty)
    yyPutBS = x*expmrf_phidpl/(4.0*sqrty**3)*(dpl*dm - 1.0)
    xyPutBS = -1.0*expmrf_phidpl*dm/(2.0*y)
    
    res = (PutBS, xxPutBS, yyPutBS, xyPutBS)
    
    return res







# Example. Uncomment the following code then run


# x = 100
# y = 10**(-6)
# Strk = 100

# rd2 = 0.02
# rd1 = 0.01
# rf2 = 0.00
# rf1 = 0.00
# dt2 = 1/12.0
# dt1 = 1/12.0

# rd_deque = deque([rd2, rd1])
# rf_deque = deque([rf2, rf1])
# dt = deque([dt2, dt1])

# print(parPutBS_pw(x, y, Strk, rd_deque, rf_deque, dt))