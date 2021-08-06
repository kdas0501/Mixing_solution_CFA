# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 15:20:57 2017

@author: Kaustav
"""

import copy as cp
from math import exp
from collections import deque


# Computes the Phi function from the mixing solution paper over T_i to T_{i+1}
# for general arguments.

# params (deque): deque([Kn, K(n-1), ..., K1, DT]) - the K_i make up
#   the 'first' argument of Phi, and DT = T_{i+1} - T_i.
# ints (deque): deque([pn,  p(n-1), ...,  p2, p1]) - 'second' argument of Phi.
# eps (float): constant numeric precision for Taylor expansions. eps = 10**(-3)
#   is fine!


def PhiG(_ints, _params, eps):
    
    
    # Copy deques
    ints = cp.copy(_ints)
    params = cp.copy(_params)
    
    origints = cp.copy(_ints)
    origparams = cp.copy(_params) 
    
    p = ints.popleft()
    K = params.popleft()
    DT = params[-1]
    
    lastlayer = ints == deque([])
    
    # For first parameter Kn approximately 0
    if abs(K) < 10**(-8):
        if lastlayer:
            res = DT/(p+1.0)
        else:
            newints = cp.copy(ints)
            newints[0] += p + 1.0
            res = (PhiG(ints, params, eps) - PhiG(newints,params,eps) )*DT/(p + 1.0)
            
    # For first int pn = 0
    elif p == 0:
        # If smallish K, Taylor expansion for the exp() term, as division by 
        # K maybe (or maybe not...) causes numerical errors!
        if abs(K) < 10**(-2):
            if lastlayer:
                res = PhiG(deque([0]), deque([0,DT]), eps)
                aux = 1.0
                m = 1
                while aux > eps:
                    aux *= K*DT/m
                    res += aux*PhiG(deque([m]), deque([0,DT]), eps)
                    m += 1
            else: 
                origparams[0] = 0.0
                res = PhiG(origints, origparams, eps)
                aux = 1.0
                m = 1
                while aux > eps:
                    aux *= K*DT/m
                    origints[0] += 1
                    res += aux*PhiG(origints, origparams, eps)
                    m += 1
        # If K not small, calculate normally
        else:
            if lastlayer:
                res = (exp(K*DT)-1.0)/K
            else:
                newparams = cp.copy(params)
                newparams[0] += K
                res = (PhiG(ints, params, eps)*exp(K*DT) - PhiG(ints, newparams, eps))/K
                
    # From here pn > 0 and Kn > 0
            
    # If smallish Kn, Taylor expansion for the exp() term, as division by 
    # Kn maybe (or maybe not...) causes numerical errors!
    elif abs((K)**(p + 1.0)) < 10**(-2):
        if lastlayer:
            res = PhiG(deque([p]), deque([0,DT]), eps)
            aux = 1.0
            m = 1
            while aux > eps:
                aux *= K*DT/m
                res += aux*PhiG(deque([p+m]), deque([0,DT]), eps)
                m += 1
        else:
            origparams[0] = 0.0
            res = PhiG(origints, origparams, eps)
            aux = 1.0
            m = 1
            while aux > eps:
                aux *= K*DT/m
                origints[0] += 1
                res += aux*PhiG(origints, origparams, eps)
                m += 1
    # If Kn not small, calculate normally
    else:
        if lastlayer:
            res = (exp(K*DT) - p/DT*PhiG(deque([p-1]), deque([K, DT]), eps))/K
        else:
            newints = cp.copy(ints)
            newparams = cp.copy(params)
            newints[0] += p
            newparams[0] += K
            origints[0] += -1
            res = (exp(DT*K)*PhiG(ints, params, eps) \
                   - p/DT*PhiG(origints, origparams, eps) \
                       - PhiG(newints, newparams, eps))/K
          
    return res










# Example. Uncomment the following code then run



# ints   = deque([3,0,2,0,1])
# params = deque([-1.0, -1.0, 0.17, 0.8, 1.0, 1/12.0])


# print(PhiG(ints, params, 10**(-3.0)))









