# Tien Dat, Truong
# Group 14
# C-Exercise 08

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
#a
def CRR_EuCall (S_0, r, sigma, T, M, K):
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
    # create EU call return matrix
    V= np.empty((M + 1, M + 1))
    for j in range (1, M+2):
        V[j-1,M]= np.maximum(0,(S[j-1,M]-K))
    for i in range(M+1,1,-1):
        for j in range(1, i + 1, 1):
            V[j - 2, i - 2]= np.exp(-r*delta_t)*(q*V[j-1,i-1] + (1-q)*V[j-2,i-1])
    return V[0,0]
#b
def CRR_EuCall_altcond (S_0, r, sigma, T, M, K):
    # compute values of u, d and q
    delta_t = T / M
    alpha = math.exp(r * delta_t)
    gamma= (K/S_0)**(2/M)
    beta= 1 / 2 * (gamma / alpha + alpha * math.exp(math.pow(sigma, 2) * delta_t))
    u = beta + math.sqrt(math.pow(beta, 2) - gamma)
    d = gamma / u
    q = (math.exp(r * delta_t) - d) / (u - d)
    # create return matrix S
    S = np.empty((M + 1, M + 1))
    # calculate value and fill in the matrix
    for i in range(1, M + 2, 1):
        for j in range(1, i + 1, 1):
            S[j - 1, i - 1] = S_0 * math.pow(u, j - 1) * math.pow(d, i - j)
    #create EU call return matrix
    V= np.empty((M + 1, M + 1))
    for j in range (1, M+2):
        V[j-1,M]= np.maximum(0,(S[j-1,M]-K))
    for i in range(M+1,1,-1):
        for j in range(1, i + 1, 1):
            V[j - 2, i - 2]= np.exp(-r*delta_t)*(q*V[j-1,i-1] + (1-q)*V[j-2,i-1])
    return V[0,0]
#c
def BlackScholes_EuCall (t, S_t, r, sigma, T, K):
    #calculate d1,d2
    d1= (np.log(S_t/K) + (r + (sigma**2)/2)*(T-t))/(sigma*np.sqrt(T-t))
    d2= d1 - sigma*np.sqrt(T-t)
    EU_Call = S_t * scipy.stats.norm.cdf(d1) - K * np.exp(-r * (T - t)) * scipy.stats.norm.cdf(d2)
    return EU_Call
#d
# Input parameters
S_0 = 100
S_t=100
r = 0.03
sigma = 0.3
T = 1
M = 100
t=0
# Create result matrix and fill in
K=np.arange(70,201,1)
N=len(K)
BS_CRR= np.empty(N)
BS_CRR_alcond= np.empty(N)
for i in range (0,N):
    BS_CRR[i]= CRR_EuCall (S_0, r, sigma, T, M, K[i]) - BlackScholes_EuCall (t, S_t, r, sigma, T, K[i])
    BS_CRR_alcond[i]= CRR_EuCall_altcond (S_0, r, sigma, T, M, K[i]) - BlackScholes_EuCall (t, S_t, r, sigma, T, K[i])

# Plotting
plt.clf()
plt.plot(K,BS_CRR)
plt.plot(K,BS_CRR_alcond)
plt.show()
#e
def CRR_DownOutCall (S_0, r, sigma, T, M, K, B):
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
            # condition for down and out call
            if S[j - 1, i - 1] > B:
                S[j - 1, i - 1] = S[j - 1, i - 1]
            else:
                S[j - 1, i - 1]=0
    #create EU call return matrix
    V= np.empty((M + 1, M + 1))
    for j in range (1, M+2):
        V[j-1,M]= np.maximum(0,(S[j-1,M]-K))
    for i in range(M+1,1,-1):
        for j in range(1, i + 1, 1):
            V[j - 2, i - 2]= np.exp(-r*delta_t)*(q*V[j-1,i-1] + (1-q)*V[j-2,i-1])
    return V[0,0]

#f
# Input parameters
S_0 = 100
r = 0.03
sigma= 0.3
T = 1
M = 1000
K = 100
B = 90
# Test parameters
print(CRR_DownOutCall (S_0, r, sigma, T, M, K, B))