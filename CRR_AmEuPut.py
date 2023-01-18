# Tien Dat, Truong
# Group 14
# C-Exercise 09

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

#a
def CRR_AmEuPut(S_0, r, sigma, T, M, K, EU):
    # compute values of u, d and q
    delta_t = T / M
    alpha = math.exp(r * delta_t)
    beta = 1 / 2 * (1 / alpha + alpha * math.exp(math.pow(sigma, 2) * delta_t))
    u = beta + math.sqrt(math.pow(beta, 2) - 1)
    d = 1 / u
    q = (math.exp(r * delta_t) - d) / (u - d)
    # create return matrix S
    S = np.empty((M + 1, M + 1))
    # calculate value and fill in the matrix
    for i in range(1, M + 2, 1):
        for j in range(1, i + 1, 1):
            S[j - 1, i - 1] = S_0 * math.pow(u, j - 1) * math.pow(d, i - j)
    # create put return matrix
    V = np.empty((M + 1, M + 1))
    for j in range(1, M + 2):
        V[j - 1, M] = np.maximum(0, (K - S[j - 1, M]))
    # condition for European and American Put Options
    if (EU == 1):
        for i in range(M + 1, 1, -1):
            for j in range(1, i + 1, 1):
                V[j - 2, i - 2] = np.exp(-r * delta_t) * (q * V[j - 1, i - 1] + (1 - q) * V[j - 2, i - 1])
        return V[0,0]
    elif (EU == 0):
        for i in range(M + 1, 1, -1):
            for j in range(1, i + 1, 1):
                V[j - 2, i - 2] = np.maximum(np.maximum(K-S[j-2,i-2],0),(np.exp(-r * delta_t) * (q * V[j - 1, i - 1] + (1 - q) * V[j - 2, i - 1])))
        return V[0,0]
    else:
        print('Error: variable EU must be 0 or 1')
#b
def BlackScholes_EuPut (t, S_t, r, sigma, T, K):
    # calculate d1,d2
    d1= (np.log(S_t/K) + (r + (sigma**2)/2)*(T-t))/(sigma*np.sqrt(T-t))
    d2= d1 - sigma*np.sqrt(T-t)
    EU_Put = K * np.exp(-r * (T - t)) * scipy.stats.norm.cdf(-d2) - S_t * scipy.stats.norm.cdf(-d1)
    return EU_Put
#c
# Input parameters
S_0 = 100
r = 0.05
sigma= 0.3
T = 1
K = 120
t=0
S_t=100
# Create result matrix and fill in
M=np.arange(10,501,1)
N=len(M)
EU_Put= np.zeros(N)
BS_Put= np.zeros(N)
for i in range (0,N):
    EU_Put[i]=CRR_AmEuPut(S_0, r, sigma, T, M[i], K,1)
    BS_Put[i]=BlackScholes_EuPut (t, S_t, r, sigma, T, K)
#Plotting
plt.clf()
plt.plot(M,EU_Put)
plt.plot(M,BS_Put)
plt.show()
#Print result
print(CRR_AmEuPut(S_0, r, sigma, T, 500, K,0))

