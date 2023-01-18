# Tien Dat, Truong
# Group 14
# C-Exercise 36

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import cmath
import scipy.integrate as integrate
import scipy.interpolate

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

#Input; re-order the input equation according to tridiagonal matrices order
alpha=np.array([1,3,2])
beta= np.array([2,1,0])
gamma= np.array([0,1,1])
b=np.array([3,1,3])

# Result
x=TriagMatrix_Inversion(alpha,beta,gamma,b)
print(x)

