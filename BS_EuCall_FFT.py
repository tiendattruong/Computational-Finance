# Tien Dat, Truong
# Group 14
# C-Exercise 27


import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
import cmath
import scipy.integrate as integrate
import scipy.interpolate


def BS_EuCall_FFT (S0, r, sigma, T, K, R, N, M, kappa1):
    #Create f_tidle_0, chi, g function
    def f_tidle_0(z):
        return 1/(z*(z-1))

    def chi(u):
        return cmath.exp(complex(0,1)*u*(np.log(S0)+r*T)-(np.complex(0,1)*u+u**2)*(sigma**2)/2*T)

    def g(u):
        return f_tidle_0(R+complex(0,1)*u)*chi(u-complex(0,1)*R)
    #Create requried matrix
    delta= M/N
    m= np.arange(1,N+1,1)
    kappa = np.empty(N)
    x= np.empty(N,dtype=complex)
    V_km= np.empty(N)
    #Fill in matrix
    for i in range (0,N,1):
        kappa[i]= kappa1+(m[i]-1)*2*np.pi/M
        x[i]= g((m[i]-0.5)*delta)*delta*cmath.exp(-complex(0,1)*(m[i]-1)*delta*kappa1)
    x_hat= np.fft.fft(x,axis=-1)
    for i in range (0,N,1):
        V_km[i]= 1/np.pi*(np.exp(-r*T+(1-R)*kappa[i])) *np.real(x_hat[i]*cmath.exp(-complex(0,1)/2*delta*kappa[i]))
    # Create 1d function based on kappa and V_km
    f= scipy.interpolate.interp1d(kappa,V_km)
    # Input K to generate V0
    V0=f(np.log(K))
    return V0

def EuCall_BlackScholes(t, S_t, r, sigma, T, K):
    d_1 = (math.log(S_t / K) + (r + 1 / 2 * math.pow(sigma, 2)) * (T - t)) / (sigma * math.sqrt(T - t))
    d_2 = d_1 - sigma * math.sqrt(T - t)
    phi = scipy.stats.norm.cdf(d_1)
    C = S_t * phi - K * math.exp(-r * (T - t)) * scipy.stats.norm.cdf(d_2)
    return C

#Input parameters
S0=100
r=0.05
sigma=0.2
T=1
K= np.arange(80,181,1)
R=1.1
N=2**11
M=50
kappa1=np.log(80)

# Result
BSCall=np.empty(len(K))
for i in range(0,len(K),1):
    BSCall[i]=EuCall_BlackScholes(0, S0, r, sigma, T, K[i])
V0=BS_EuCall_FFT (S0, r, sigma, T, K, R, N, M, kappa1)

#Plotting
plt.plot(K,V0)
plt.plot(K,BSCall)
plt.xlabel("Strike Price")
plt.ylabel("Initial Price of European call option")
plt.title(" Initial price of a European call options in the Black-Scholes model via the fast Fourier transform approach")
plt.show()






