# Mixing_solution_CFA
This repository contains working code for pricing European put options with piecewise-constant parameters in the Heston model and GARCH diffusion model (GARCH diffusion model rho = 0 a.e..) using the closed-form approximation formulas from the article 'Closed-form approximations with respect to the mixing solution for option pricing under stochastic volatility'. 

Title: Closed-form approximations with respect to the mixing solution for option pricing under stochastic volatility

Authors: Kaustav Das (Monash University) and Nicolas Langrene (CSIRO, Data61)

email addresses: kaustav.das@monash.edu and nicolas.langrene@csiro.au




### Quickstart for readers of article

Simply run Heston_compare_pw.py, which computes and compares the price and implied volatility of a European put option in the Heston model with piecewise-constant parameter inputs, where the price is obtained via the mixing solution closed-form approximation method and mixing solution Monte-Carlo method.



### Main files 
This repository contains the following .py files:

- **GARCH_CFA_mixing.py:**
  Computes the price of a European put option in the GARCH diffusion model with piecewise-constant parameter inputs via the mixing solution closed-form approximation method.
  
- **Heston_CFA_mixing.py:**
    Computes the price of a European put option in the Heston model with piecewise-constant parameter inputs via the mixing solution closed-form approximation method.
    
- **Omg.py:**
  Contains functions Omg1, Omg2, Omg3, Omg4, Omg5 which compute the integral operator (1 to 5 fold) for piecewise-constant parameter inputs.
  
- **PhiG.py:**
  Computes the n-fold function Phi for any n.
  
- **parPutBS_pw.py:**
  Computes the second-order partial derivatives of the function Put_BS for piecewise-constant parameter inputs.
  




### Auxiliary files
The rest of the .py files are auxiliary files that are not required for the closed-form approximation formula.

  - **Monte_mixing_pw.py:**
    Computes the price of a European put option in the Heston, GARCH diffusion, Ornstein-Uhlenbeck, Inverse-Gamma, and Verhulst models with piecewise-constant parameter inputs via the mixing solution Monte-Carlo method.
  
- **GARCH_compare_pw.py:**
    Compares the price and implied volatility of a European put option in the GARCH diffusion model with piecewise-constant parameter inputs, where the price is obtained via the mixing solution closed-form approximation method and mixing solution Monte-Carlo method.
   
- **Heston_compare_pw.py:**
    Compares the price and implied volatility of a European put option in the Heston model with piecewise-constant parameter inputs, where the price is obtained via the mixing solution closed-form approximation method and mixing solution Monte-Carlo method.
    
- **BSform_pw.py:** 
  Computes the usual Black-Scholes price of a European Put/Call option for piecewise-constant parameter inputs.
  
- **ImpVolBrent_pw.py:** 
  Computes the implied volatility of a European put/call option with piecewise-constant parameter inputs via Brent's method. 
  
- **DeltaStrikes_pw.py:**
  Computes the strike of an option contract corresponding to a given European put/call option Delta with piecewise-constant parameter inputs.
