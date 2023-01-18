# Tien Dat, Truong
# Group 14
# C-Exercise 35

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import cmath
import scipy.integrate as integrate
import scipy.interpolate

def BS_EuCall_FiDi_Explicit (r, sigma, a, b, m, nu_max, T, K):
    #input
    q = 2*r/sigma**2
    deltax_tilde = (b-a)/m
    deltat_tilde = (sigma**2)*T/(2*nu_max)
    lamda= deltat_tilde/deltax_tilde**2
    #create matrix
    x =np.arange(a, b, deltax_tilde)
    w= np.empty((m,nu_max))
    S= np.empty(m)
    V0= np.empty(m)
    #fill in matrix
    for i in range (0,m,1):
        w[i,0]= np.maximum(np.exp(x[i]/2*(q+1))-np.exp(x[i]/2*(q-1)), 0)

    for i in range (1,nu_max,1):
        for j in range (1,m-1,1):
            w[j,i] = lamda * w[j-1,i-1] +(1-2*lamda) * w[j,i-1] + lamda* w[j+1,i-1]
    for i in range (0,m,1):
        S[i]= K * np.exp(x[i])
        V0[i]= K*w[i,nu_max-1]*np.exp((-x[i]/2)*(q-1)-sigma**2/2*T*((q-1)**2/4+q))
    return S,V0
# Input parameters
r = 0.05
sigma = 0.2
a = -0.7
b = 0.4
m = 100
nu_max = 2000
T = 1
K = 100

result = BS_EuCall_FiDi_Explicit(r, sigma, a, b, m, nu_max, T, K)


def EuCall_BlackScholes(t, S_t, r, sigma, T, K):
    d_1 = (math.log(S_t / K) + (r + 1 / 2 * math.pow(sigma, 2)) * (T - t)) / (sigma * math.sqrt(T - t))
    d_2 = d_1 - sigma * math.sqrt(T - t)
    Call = S_t * scipy.stats.norm.cdf(d_1) - K * math.exp(-r * (T - t)) * scipy.stats.norm.cdf(d_2)
    return Call

# result of BS-method
V0_BS= np.empty(len(result[0]))
for i in range (0,len(result[0]),1):
    V0_BS[i]=EuCall_BlackScholes(0, result[0][i], r, sigma, T, K)

#Plotting
plt.clf()
plt.plot(result[0],result[1],label= 'EU_Call_Fidi')
plt.plot(result[0],V0_BS, label= 'EU_Call_BS')
plt.xlabel('S')
plt.ylabel('VO')
plt.legend(loc='best')
plt.show()