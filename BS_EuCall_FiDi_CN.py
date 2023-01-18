# Tien Dat, Truong
# Group 14
# C-Exercise 38

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import cmath
import scipy.integrate as integrate
import scipy.interpolate

# recall function
def TriagMatrix_Inversion(alpha,beta,gamma,b):
    #create matrix
    n=len(alpha)
    alpha_hat=np.empty(n)
    b_hat= np.empty(n)
    x= np.empty(n)
    #first value
    alpha_hat[0]=alpha[0]
    b_hat[0]=b[0]
    #fill in matrix
    for i in range (1,n):
        alpha_hat[i]= alpha[i]-(gamma[i]/alpha_hat[i-1])*beta[i-1]
        b_hat[i]= b[i]-(gamma[i]/alpha_hat[i-1])*b_hat[i-1]

    x[n-1]= b_hat[n-1]/alpha_hat[n-1]
    for i in range (n-2,-1,-1):
        x[i]= (1/alpha_hat[i])* (b_hat[i]-beta[i]*x[i+1])
    return x

def BS_EuCall_FiDi_CN (r, sigma, a, b, m, nu_max, T, K):
    # input
    q = 2 * r / sigma ** 2
    deltax_tilde = (b - a) / m
    deltat_tilde = (sigma ** 2) * T / (2 * nu_max)
    lamda = deltat_tilde / deltax_tilde ** 2

    # create matrix
    x = np.arange(a, b, deltax_tilde)
    w = np.empty((m, nu_max))
    S = np.empty(m)
    V0 = np.empty(m)

    # fill in matrix
    for i in range(0, m, 1):
        w[i,0]= np.maximum(np.exp(x[i]/2*(q+1))-np.exp(x[i]/2*(q-1)), 0)

    #create matrix A
    alpha= np.ones(m - 2)*(1 + lamda)
    beta = np.ones(m - 2) * (-lamda/2)
    gamma = np.ones(m - 2) * (-lamda/2)

    #fill in matrix
    for i in range (1,nu_max,1):
        for j in range (1,m-1,1):
            w[j,i] = lamda/2 * w[j-1,i-1] +(1-lamda) * w[j,i-1] + lamda/2* w[j+1,i-1]
        w[m-2,i] = w[m-2,i] + lamda/2 * (np.exp((q+1)/2*b+(q+1) ** 2 / 4 * i * deltat_tilde) - np.exp((q - 1) / 2 * b + (q - 1) ** 2 / 4 * i * deltat_tilde))
        w[1:-1,i] = TriagMatrix_Inversion(alpha, beta, gamma, w[1:-1, i])
        w[m-1, i] = np.exp((q + 1) / 2 * b + (q + 1) ** 2 / 4 * i * deltat_tilde) - np.exp((q - 1) / 2 * b + (q - 1) ** 2 / 4 * i * deltat_tilde)
    for i in range (0,m,1):
        S[i]= K * np.exp(x[i])
        V0[i]= K*w[i,nu_max-1]*np.exp((-x[i]/2)*(q-1)-sigma**2/2*T*((q-1)**2/4+q))
    return [S, V0]

# Input parameters
r = 0.05
sigma = 0.2
a = -0.7
b = 0.4
m = 100
nu_max = 2000
T = 1
K = 100


result = BS_EuCall_FiDi_CN(r, sigma, a, b, m, nu_max, T, K)


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
plt.plot(result[0],result[1],label= 'EU_Call_CN')
plt.plot(result[0],V0_BS, label= 'EU_Call_BS')
plt.xlabel('S')
plt.ylabel('VO')
plt.legend(loc='best')
plt.show()