# Tien Dat, Truong
# Group 14
# C-Exercise 14

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

def BS_EuCall_MC_IS (S0, r, sigma, K, T, mu, N, alpha):
    # Create function for a call option
    def f(x):
        y = np.maximum(0, x - K)
        return y
    # generate  N samples
    Y = np.random.normal(mu, 1, N)
    ST = np.empty(len(Y))
    X= np.empty(len(Y))
    # Compute ST, V0
    for i in range(0, len(Y)):
        ST[i] = S0 * np.exp((r - 0.5 * sigma ** 2) * T + sigma * np.sqrt(T) * Y[i])
        X[i]= f(ST[i])
    V0= np.mean(np.exp(-r*T-Y*mu+mu**2/2)*X)
    #Compute confidence interval
    eps= scipy.stats.norm.ppf(1-(1-alpha)/2)*np.sqrt(np.var(np.exp(-r*T-Y*mu+mu**2/2)*X)/N)
    CIl = V0 - eps
    CIr = V0 + eps
    return V0, CIl, CIr

# Test parameters
S0 = 100
r = 0.05
sigma = 0.3
K = 220
T = 1
N = 10000
alpha = 0.95

# BS function
def BlackScholes_EuCall (t, S_t, r, sigma, T, K):
    #calculate d1,d2
    d1= (np.log(S_t/K) + (r + (sigma**2)/2)*(T-t))/(sigma*np.sqrt(T-t))
    d2= d1 - sigma*np.sqrt(T-t)
    EU_Call = S_t * scipy.stats.norm.cdf(d1) - K * np.exp(-r * (T - t)) * scipy.stats.norm.cdf(d2)
    return EU_Call
# Choose mu
d= (np.log(K/S0)-(r-0.5*sigma**2)*T)/(sigma*np.sqrt(T))
mu = np.arange(0,d+abs(d),0.01)
V0= np.empty(len(mu))
for i in range (0,len(mu)):
    V0[i]= BS_EuCall_MC_IS (S0, r, sigma, K, T, mu[i], N, alpha)[0]- BlackScholes_EuCall (0, S0, r, sigma, T, K)

# Plotting
plt.clf()
plt.plot(mu,V0)
plt.show()

