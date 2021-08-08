# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 00:16:43 2021

@author: Dr Kaustav Das (kaustav.das@monash.edu)
"""


import copy as cp
from collections import deque
from math import exp
from PhiG import PhiG


# Computes e^{k dot dt}, where k and dt are deques.

def Expo(_k, _dt):
    
    
    # Copy deques
    k = cp.copy(_k)
    dt = cp.copy(_dt)
        
    lastlayer = k == deque([])
    
    if lastlayer:
        return 1
    else:
        K = k.popleft()
        DT = dt.popleft()
    
    res = Expo(k, dt)*exp(DT*K)
    
    return res


# Example code.
    
# if __name__ == '__main__':

    
#     dt = deque([1, 2])
#     k = deque([0.5, 3])
    
#     print(Expo(k, dt))










    
# Precision for Taylor expansions in PhiG
eps = 10**(-3)


# 1D Omega function from mixing solution paper. 
# k1 and l1 are deques of piecewise-constant parameters and dt is the 
# time discretisation deque

def Omg1(_k1, _l1, _dt):
    
    
    # Copy deques
    k1 = cp.copy(_k1)
    l1 = cp.copy(_l1)
    dt = cp.copy(_dt)
    
    lastlayer = k1 == deque([])
    if lastlayer:
        return 0
    else:
        K1 = k1.popleft()
        L1 = l1.popleft()
        DT = dt.popleft()
         
        res = Omg1(k1, l1, dt)+L1*Expo(k1, dt)*PhiG(deque([0]), deque([K1,DT]), eps)
        
    return res


# Example code.
    
if __name__ == '__main__':

    
    k1 = deque([0.3,0.2,0.1,0.1])
    l1 = deque([0.5,0.3,0.1,0.1])
    dt = deque([0.3,0.05,0.2,0.01])
        
    print(Omg1(k1,l1,dt))














# 2D Omega function from mixing solution paper. 
# k1, k2 and l1, l2 are deques of piecewise-constant parameters and dt is the 
# time discretisation deque

def Omg2(_k2, _k1, _l2, _l1, _dt):
    
    
    # Copy deques
    k2 = cp.copy(_k2)
    k1 = cp.copy(_k1)
    l2 = cp.copy(_l2)
    l1 = cp.copy(_l1)
    dt = cp.copy(_dt)
    
    lastlayer = k2 == deque([])
    if lastlayer:
        return 0
    else:
        K2 = k2.popleft()
        K1 = k1.popleft()
        L2 = l2.popleft()
        L1 = l1.popleft()
        DT = dt.popleft()
        
        Exp2 = Expo(k2,dt)
        Exp21 = Exp2*Expo(k1,dt)
        
        res = Omg2(k2,k1,l2,l1,dt)+ \
        L2*Exp2*Omg1(k1,l1,dt)*PhiG(deque([0]), deque([K2,DT]),10**(-8.0))+\
        L2*L1*Exp21*PhiG(deque([0,0]), deque([K2,K1,DT]),eps)
        
    return res


# Example code.

# if __name__ == '__main__':

    
#     k2 = deque([0.3,0.5,0.4,0.1])
#     k1 = deque([0.2,2.3,0.3,2])
#     l2 = deque([0.9,0.5,0.1,0.5])
#     l1 = deque([0.3,1,3,0.3])
#     dt = deque([0.5,0.2,0.3,0.5])
    
#     print(Omg2(k2,k1,l2,l1,dt))














# 3D Omega function from mixing solution paper. 
# k1, k2, k3 and l1, l2, l3 are deques of piecewise-constant parameters 
# and dt is the time discretisation deque

def Omg3(_k3, _k2, _k1, _l3, _l2, _l1, _dt): 
    
    
    # Copy deques
    k3 = cp.copy(_k3)
    k2 = cp.copy(_k2)
    k1 = cp.copy(_k1)
    l3 = cp.copy(_l3)
    l2 = cp.copy(_l2)
    l1 = cp.copy(_l1)
    dt = cp.copy(_dt)
    
    lastlayer = k1 == deque([])
    
    if lastlayer:
        return 0
    else:
        K3 = k3.popleft()
        K2 = k2.popleft()
        K1 = k1.popleft()
        L3 = l3.popleft()
        L2 = l2.popleft()
        L1 = l1.popleft()
        DT = dt.popleft()
        
        L3L2 = L3*L2
        L3L2L1 = L3L2*L1
        Exp3 = Expo(k3, dt)
        Exp32 = Exp3*Expo(k2,dt)
        Exp321 = Exp32*Expo(k1,dt)
        
        res = Omg3(k3,k2,k1,l3,l2,l1,dt)+\
        L3*Exp3*Omg2(k2,k1,l2,l1,dt)*PhiG(deque([0]), deque([K3,DT]),eps)+\
        L3L2*Exp32*Omg1(k1,l1,dt)*\
        PhiG(deque([0,0]), deque([K3,K2,DT]),eps)+L3L2L1*Exp321*PhiG(deque([0,0,0]), deque([K3,K2,K1,DT]),eps) 
    
    return res
    

# Example code.

# if __name__ == '__main__':


#     k3 = deque([0.2,0.5,0.1,0.5])    
#     k2 = deque([0.5,0.3,0.2,0.5])
#     k1 = deque([0.3,0.2,0.1,0.1])
#     l3 = deque([0.2,0.4,1,1.2])
#     l2 = deque([0.3,0.2,0.1,0.5])
#     l1 = deque([0.5,0.3,0.1,0.1])
#     dt = deque([0.3,0.05,0.2,0.01])
    
#     print(Omg3(k3,k2,k1,l3,l2,l1,dt))














# 4D Omega function from mixing solution paper. 
# k1, k2, k3, k4 and l1, l2, l3, l4 are deques of piecewise-constant parameters
# and dt is the time discretisation deque

def Omg4(_k4, _k3, _k2, _k1,_l4, _l3, _l2, _l1, _dt):
    
    
    # Copy deques
    k4 = cp.copy(_k4)
    k3 = cp.copy(_k3)
    k2 = cp.copy(_k2)
    k1 = cp.copy(_k1)
    l4 = cp.copy(_l4)
    l3 = cp.copy(_l3)
    l2 = cp.copy(_l2)
    l1 = cp.copy(_l1)
    dt = cp.copy(_dt)
    
    lastlayer = k1 == deque([])
    
    if lastlayer:
        return 0
    else:
        K4 = k4.popleft()
        K3 = k3.popleft()
        K2 = k2.popleft()
        K1 = k1.popleft()
        L4 = l4.popleft()
        L3 = l3.popleft()
        L2 = l2.popleft()
        L1 = l1.popleft()
        DT = dt.popleft()
        
        L4L3 = L4*L3
        L4L3L2 = L4L3*L2
        L4L3L2L1 = L4L3L2*L1
        Exp4 = Expo(k4,dt)
        Exp43 = Exp4*Expo(k3,dt)
        Exp432 = Exp43*Expo(k2,dt)
        Exp4321 = Exp432*Expo(k1,dt)
    
        res = Omg4(k4,k3,k2,k1,l4,l3,l2,l1,dt)+\
        L4*Exp4*Omg3(k3,k2,k1,l3,l2,l1,dt)*PhiG(deque([0]), deque([K4,DT]),eps)+\
        L4L3*Exp43*Omg2(k2,k1,l2,l1,dt)*PhiG(deque([0,0]), deque([K4,K3,DT]),eps)+\
        L4L3L2*Exp432*Omg1(k1,l1,dt)*PhiG(deque([0,0,0]), deque([K4,K3,K2,DT]),eps)+\
        L4L3L2L1*Exp4321*PhiG(deque([0,0,0,0]), deque([K4,K3,K2,K1,DT]),eps) 
    
    return res


# Example code.

# if __name__ == '__main__':


#     k4 = deque([0.7,0.1,0.9,0.2])     
#     k3 = deque([0.2,0.5,0.1,0.5])    
#     k2 = deque([0.5,0.3,0.2,0.5])
#     k1 = deque([0.3,0.2,0.1,0.1])
#     l4 = deque([0.1,0.9,1.5,2])
#     l3 = deque([0.2,0.4,1,1.2])
#     l2 = deque([0.3,0.2,0.1,0.5])
#     l1 = deque([0.5,0.3,0.1,0.1])
#     dt = deque([0.3,0.05,0.2,0.01])
    
#     print(Omg4(k4,k3,k2,k1,l4,l3,l2,l1,dt))














# 5D Omega function from mixing solution paper. 
# k1, k2, k3, k4, k5 and l1, l2, l3, l4, l5 are deques of piecewise-constant 
# parameters and dt is the time discretisation deque

def Omg5(_k5, _k4, _k3, _k2, _k1,_l5, _l4, _l3, _l2, _l1, _dt):  
    
    
    # Copy deques
    k5 = cp.copy(_k5)
    k4 = cp.copy(_k4)
    k3 = cp.copy(_k3)
    k2 = cp.copy(_k2)
    k1 = cp.copy(_k1)
    l5 = cp.copy(_l5)
    l4 = cp.copy(_l4)
    l3 = cp.copy(_l3)
    l2 = cp.copy(_l2)
    l1 = cp.copy(_l1)
    dt = cp.copy(_dt)
    
    lastlayer = k1 == deque([])
    
    if lastlayer:
        return 0
    else:
        K5 = k5.popleft()
        K4 = k4.popleft()
        K3 = k3.popleft()
        K2 = k2.popleft()
        K1 = k1.popleft()
        L5 = l5.popleft()
        L4 = l4.popleft()
        L3 = l3.popleft()
        L2 = l2.popleft()
        L1 = l1.popleft()
        DT = dt.popleft()
        
        L5L4 = L5*L4
        L5L4L3 = L5L4*L3
        L5L4L3L2 = L5L4L3*L2
        L5L4L3L2L1 = L5L4L3L2*L1
        Exp5 = Expo(k5, dt)
        Exp54 = Exp5*Expo(k4,dt)
        Exp543 = Exp54*Expo(k3,dt)
        Exp5432 = Exp543*Expo(k2,dt)
        Exp54321 = Exp5432*Expo(k1,dt)
    
        res = Omg5(k5,k4,k3,k2,k1,l5,l4,l3,l2,l1,dt)+\
        L5*Exp5*Omg4(k4,k3,k2,k1,l4,l3,l2,l1,dt)*PhiG(deque([0]), deque([K5,DT]),eps)+\
        L5L4*Exp54*Omg3(k3,k2,k1,l3,l2,l1,dt)*PhiG(deque([0,0]), deque([K5,K4,DT]),eps)+\
        L5L4L3*Exp543*Omg2(k2,k1,l2,l1,dt)*PhiG(deque([0,0,0]), deque([K5,K4,K3,DT]),eps)+\
        L5L4L3L2*Exp5432*Omg1(k1,l1,dt)*PhiG(deque([0,0,0,0]), deque([K5,K4,K3,K2,DT]),eps)+\
        L5L4L3L2L1*Exp54321*PhiG(deque([0,0,0,0,0]), deque([K5,K4,K3,K2,K1,DT]),eps) 
    
    return res


# Example code.

# if __name__ == '__main__':


#     k5 = deque([0.7,0.1,0.9,0.2])         
#     k4 = deque([0.7,0.1,0.9,0.2])     
#     k3 = deque([0.2,0.5,0.1,0.5])    
#     k2 = deque([0.5,0.3,0.2,0.5])
#     k1 = deque([0.3,0.2,0.1,0.1])
#     l5 = deque([0.5,0.3,0.1,0.1])
#     l4 = deque([0.1,0.9,1.5,2])
#     l3 = deque([0.2,0.4,1,1.2])
#     l2 = deque([0.3,0.2,0.1,0.5])
#     l1 = deque([0.5,0.3,0.1,0.1])
#     dt = deque([0.3,0.05,0.2,0.01])
    
#     print(Omg5(k5,k4,k3,k2,k1,l5,l4,l3,l2,l1,dt))