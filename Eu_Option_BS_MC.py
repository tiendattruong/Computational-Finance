# Tien Dat, Truong
# Group 14
# C-Exercise 12

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

# Create function for a call option
def f(x):
    y= np.maximum(0,x-K)
    return y
# Create main function
def Eu_Option_BS_MC (S0, r, sigma, T, K, M, f):
    # Create X is N(0;1)-distributed random variable M times
    X= np.random.normal(0,1,M)
    # Using formula
    ST=S0*np.exp((r-0.5*sigma**2)*T+sigma*np.sqrt(T)*X)
    # Follow note 2.1-2.6, we have
    H = np.exp(-r * T) * f(ST)
    V0 = np.mean(H)
    eps = 1.96 * 1 / (np.sqrt(M)) * np.sqrt(np.var(H))
    c1 = V0 - eps
    c2 = V0 + eps
    return V0,c1,c2
# Computes the price of a call in the Black-Scholes model
def BlackScholes_EuPut(t, S_t, r, sigma, T, K):
    d_1 = (math.log(S_t / K) + (r + 1 / 2 * math.pow(sigma, 2)) * (T - t)) / (sigma * math.sqrt(T - t))
    d_2 = d_1 - sigma * math.sqrt(T - t)
    C = S_t * scipy.stats.norm.cdf(d_1) - K * math.exp(-r * (T - t)) * scipy.stats.norm.cdf(d_2)
    return C
# Test parameters
S0 = 100
r = 0.05
sigma = 0.3
T = 1
M = 10000
K =90
print(Eu_Option_BS_MC (S0, r, sigma, T, K, M, f))
print(BlackScholes_EuPut(0, S0, r, sigma, T, K))
