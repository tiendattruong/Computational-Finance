# Tien Dat, Truong
# Group 14
# C-Exercise 32

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import cmath
import scipy.integrate as integrate
import scipy.interpolate

def Sim_Paths_GeoBM(X0, mu, sigma, T, N):
    delta_T= T/N
    delta_W= np.random.normal(0,np.sqrt(delta_T),N)
    #Create matrix
    X_exact = np.empty(N)
    X_Euler = np.empty(N)
    X_Milshtein= np.empty(N)
    #Fill in matrix
    for i in range (0,N-1,1):
        X_exact[i] = X0
        X_Euler[i] = X0
        X_Milshtein[i] = X0
    for i in range (1,N,1):
        X_exact[i] = X_exact[i-1] * np.exp((mu - sigma**2/2)*delta_T + sigma * delta_W[i-1])
        X_Euler[i] = X_Euler[i-1] * (1 + mu * delta_T + sigma * delta_W[i-1])
        X_Milshtein[i] = X_Milshtein[i-1] * (1+mu*delta_T + sigma*delta_W[i-1] + sigma**2/2*(delta_W[i-1]**2- delta_T))
    return X_exact,X_Euler,X_Milshtein

#Parameters
X0=100
mu=0.1
sigma=0.3
T=1
N=(10,100,1000,10000)

#Plotting
t1=np.arange(0,1,T/N[0])
t2=np.arange(0,1,T/N[1])
t3=np.arange(0,1,T/N[2])
t4=np.arange(0,1,T/N[3])

plt.clf()
plt.subplot(2,2,1)
plt.plot(t1,Sim_Paths_GeoBM(X0, mu, sigma, T, N[0])[0])
plt.plot(t1,Sim_Paths_GeoBM(X0, mu, sigma, T, N[0])[1])
plt.plot(t1,Sim_Paths_GeoBM(X0, mu, sigma, T, N[0])[2])

plt.subplot(2,2,2)
plt.plot(t2,Sim_Paths_GeoBM(X0, mu, sigma, T, N[1])[0])
plt.plot(t2,Sim_Paths_GeoBM(X0, mu, sigma, T, N[1])[1])
plt.plot(t2,Sim_Paths_GeoBM(X0, mu, sigma, T, N[1])[2])

plt.subplot(2,2,3)
plt.plot(t3,Sim_Paths_GeoBM(X0, mu, sigma, T, N[2])[0])
plt.plot(t3,Sim_Paths_GeoBM(X0, mu, sigma, T, N[2])[1])
plt.plot(t3,Sim_Paths_GeoBM(X0, mu, sigma, T, N[2])[2])

plt.subplot(2,2,4)
plt.plot(t4,Sim_Paths_GeoBM(X0, mu, sigma, T, N[3])[0])
plt.plot(t4,Sim_Paths_GeoBM(X0, mu, sigma, T, N[3])[1])
plt.plot(t4,Sim_Paths_GeoBM(X0, mu, sigma, T, N[3])[2])
plt.show()