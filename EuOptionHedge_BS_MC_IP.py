# Tien Dat, Truong
# Group 14
# C-Exercise 30

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import cmath
import scipy.integrate as integrate
import scipy.interpolate

#create call function
def g(x):
    return np.max(x-90,0)

def EuOptionHedge_BS_MC_IP (St, r, sigma, g, T, t, N):
    #create matrix
    X= np.random.normal(0,1,N)
    ST=np.empty(N)
    Z_prime = np.empty(N)
    #fill in matrix
    for i in range (0,N,1):
        ST[i]= St*np.exp((r-(sigma**2)/2)*(T-t)+sigma*np.sqrt(T-t)*X[i])
        Z_prime[i]= np.exp(-(sigma**2)/2*(T-t)+sigma*np.sqrt(T-t)*X[i])*(g(ST[i])>0)
    phi1_t= np.mean(Z_prime)
    return phi1_t
#parameters
t = 0; St = 100; r = 0.05; sigma = 0.2; T = 1; N = 10000;

print("Hedging positon phi1_t is:" + str(EuOptionHedge_BS_MC_IP (St, r, sigma, g, T, t, N)))
#Comparing hedging value
X=(np.log(St/90)+r*(T-t)+(sigma**2)/2*(T-t))/(sigma*np.sqrt(T-t))
norm_cdf=scipy.stats.norm.cdf(X,0,1)
print("Formula 3.30 value is: " +str(norm_cdf))










