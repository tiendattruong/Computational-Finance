# Tien Dat, Truong
# Group 14
# C-Exercise 31

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import cmath
import scipy.integrate as integrate
import scipy.interpolate


def SimPath_Ito_Euler(X0, a, b, T, m, N):
    delta_T= T/m
    delta_W= np.random.normal(0,np.sqrt(delta_T),(N,m))
    #Create matrix
    Y= np.empty((N,m))
    #Fill in matrix
    for i in range(0,N-1,1):
        Y[i, 0] = X0
    for i in range(1,m,1):
        Y[:, i]= Y[:, i-1]+ a(Y[:, i-1],(i-1)*delta_T)*delta_T+ b(Y[:, i-1],(i-1)*delta_T)*delta_W[:,i-1]
    return Y

#Parameters
N = 10000
m = 100
nu0 = 0.3**2
lamda = 2.5
kappa = 0.3**2
sigma_tilde = 0.2
T = 1

# Create "sufficiently smooth” coefficients a(x,t),b(x,t)
def a(x,t):
    return kappa-lamda*x

def b(x,t):
    return np.sqrt(x)*sigma_tilde

# Plotting
t=np.arange(0,1,T/m)
plt.clf()
plt.plot(t,SimPath_Ito_Euler(nu0, a, b, T, m, N)[0,:])
plt.xlabel("t")
plt.ylabel("X(t)")
plt.show()
