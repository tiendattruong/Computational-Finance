# Tien Dat, Truong
# Group 14
# C-Exercise 11

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats



def Sample_TruncExp_AR (a, b,lamda,N):
    beta=1/lamda
    x= np.random.uniform(0,1,N)
    y = np.random.exponential(beta,N)
    def f(x):
        return (lamda*np.exp(-lamda*x))/(np.exp(-lamda*a)-np.exp(-lamda*b))
    fx= np.empty(N)
    sample=[]
    for i in range (0,N-1):
        if fx[i]<=1 and fx[i]>=0:
            sample.append(fx[i])
    return x,fx

a = 1
b = 3
lamda = 2
N= 10000

c = Sample_TruncExp_AR (a, b, lamda, N)

plt.clf()
plt.plot(c[0],c[1])
plt.show()