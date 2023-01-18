# Tien Dat, Truong
# Group 14
# C-Exercise 26

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import cmath
import scipy.integrate as integrate

def Heston_PCall_Laplace (S0, r, nu0, kappa, lamda, sigma_tilde,T, K, R, p):
    #Create f_tilde, chi and intergrand function
    def f_tilde(z):
        return (K**(1-z))/(z*(z-1))
    def chi(u):
        a = cmath.exp(complex(0,u*(np.log(S0)+r*T)))
        d_u = np.sqrt(lamda**2+(sigma_tilde**2)*complex(u**2,u))
        s= np.sinh(d_u*T/2)
        c= np.cosh(d_u*T/2)
        m= a*((np.exp(lamda*T/2)/(c+lamda*s/d_u))**(2*kappa/(sigma_tilde**2)))\
           *np.exp((-nu0*((complex(u**2,u))*s/d_u)/(c+lamda*s/d_u)))
        return m
    def intergrand(x):
        return np.exp(-r*T)/math.pi*np.real(f_tilde(complex(R,x))*chi(complex(x,-R)))
    I = integrate.quad(intergrand,0,50)
    return I[0]

#Input Parameters
r=0.05
S0= np.arange(50,151,1)
nu0= 0.3**2
kappa=0.3**2
lamda = 2.5
sigma_tilde=0.2
T=1
K=100
p=1
R=1.1

#Caculating
V0=np.empty(len(S0))
for i in range (0,len(S0),1):
    V0[i]= Heston_PCall_Laplace (S0[i], r, nu0, kappa, lamda, sigma_tilde,T, K, R, p)

#Plotting
plt.clf()
plt.plot(S0,V0)
plt.xlabel("Initial stock price")
plt.ylabel("Fair price of power call")
plt.title(" Fair price of the power call in the Heston model using the Laplace transform")
plt.show()
