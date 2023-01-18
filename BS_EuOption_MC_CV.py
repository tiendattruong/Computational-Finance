# Tien Dat, Truong
# Group 14
# C-Exercise 15

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

def BS_EuOption_MC_CV (S0, r, sigma, T, K, M):
    def g(x, K):
        return max(x - K, 0)
    # generate M samples
    X = np.random.normal(0, 1, M)
    Y = np.random.normal(0, 1, M)
    ST = np.empty(len(X),dtype=float)
    ST_2 = np.empty(len(X),dtype=float)
    self_quanto = np.empty(len(X),dtype=float)
    Normal_Call = np.empty(len(X),dtype=float)
    # compute ST and Y for each sample
    for i in range(0, len(X)):
        ST[i] = S0 * math.exp((r - 0.5 * math.pow(sigma, 2)) * T + sigma * math.sqrt(T) * X[i])
        self_quanto[i] = g(ST[i], K)*ST[i]
        ST_2[i] = S0 * math.exp((r - 0.5 * math.pow(sigma, 2)) * T + sigma * math.sqrt(T) * Y[i])
        Normal_Call[i]= g(ST_2[i],K)
    #Calculate beta
    beta= np.mean((self_quanto-np.mean(self_quanto))*(Normal_Call-np.mean(Normal_Call)))/np.var(Normal_Call)
    # return result
    V0 = (np.mean(self_quanto-beta*Normal_Call)+ np.mean(beta*Normal_Call))*np.exp(-r*T)
    return V0
# Plain BS MC
def Eu_Option_BS_MC(S0, r, sigma, T, K, M, f):
    # generate M samples
    X = np.random.normal(0, 1, M)
    ST = np.empty(len(X), dtype=float)
    Y = np.empty(len(X), dtype=float)

    # compute ST and Y for each sample
    for i in range(0, len(X)):
        ST[i] = S0 * math.exp((r - 0.5 * math.pow(sigma, 2)) * T + sigma * math.sqrt(T) * X[i])
        Y[i] = f(ST[i], K)*ST[i]

    # calculate V0
    V0 = math.exp(-r * T) * np.mean(Y)

    # compute confidence interval
    epsilon = 1.96 * math.sqrt(np.var(Y) / M)
    c1 = V0 - epsilon
    c2 = V0 + epsilon
    return V0

def g(x, K):
    return max(x - K, 0)


# Test parameters
S0 = 100
r = 0.05
sigma = 0.3
T = 1
K = 110
M = 100000
print(BS_EuOption_MC_CV (S0, r, sigma, T, K, M),Eu_Option_BS_MC(S0, r, sigma, T, K, M, g))
