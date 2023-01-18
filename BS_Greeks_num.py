# Tien Dat, Truong
# Group 14
# C-Exercise 24

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import scipy.integrate as integrate

def BS_Price_Int(S0, r, sigma, T, f):
    # define integrand as given in the exercise
    def integrand(x):
        return 1 / math.sqrt(2 * math.pi) * f(
            S0 * math.exp((r - 0.5 * math.pow(sigma, 2)) * T + sigma * math.sqrt(T) * x)) * math.exp(-r * T) * math.exp(
            -1 / 2 * math.pow(x, 2))

    # perform integration
    I = integrate.quad(integrand, -np.inf, np.inf)
    # return value of the integration
    return I[0]

def BS_Greeks_num(r, sigma, S0, T, g ,eps):
    V0= BS_Price_Int(S0, r, sigma, T, g)
    Delta= (BS_Price_Int(S0+eps*S0,r,sigma,T,g) - V0)/(eps*S0)
    Gamma= (BS_Price_Int(S0+eps*S0,r,sigma,T,g) - 2*V0+ BS_Price_Int(S0-eps*S0,r,sigma,T,g))/(eps*S0)**2
    Vega= (BS_Price_Int(S0,r,sigma+eps*sigma,T,g)- BS_Price_Int(S0,r,sigma,T,g))/(eps*sigma)
    return Delta,Vega,Gamma
#Input parameters

def g(x):
    return np.maximum(x-100,0)
r=0.03
sigma=0.2
S0=np.arange(60,141,1)
T=1
eps=0.001
delta = np.empty(len(S0))
for i in range (0,len(S0),1):
    delta[i]= BS_Greeks_num(r, sigma, S0[i], T, g,eps)[0]
#Plotting
plt.clf()
plt.plot(S0,delta)
plt.show()
