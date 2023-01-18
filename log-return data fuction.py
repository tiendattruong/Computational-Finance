# Tien Dat, Truong
# Group 14
# C-Exercise 04

import math
import numpy as np
import matplotlib.pyplot as plt

#a
#Create log-return data fuction
def log_return(data):
    return np.diff(np.log(data))

#b
#Import data
dax = np.genfromtxt('time_series_dax_2021.csv', delimiter=';', usecols=4, skip_header=1)
dax = np.flip(dax)
#plot log return data

logreturn = log_return(dax)

plt.plot(logreturn)
plt.title('Daily log returns')
plt.xlabel('Trading Day')
plt.ylabel('DAX log returns')
plt.show()

#mean and standard deviation
N=len(dax)
mu= np.mean(logreturn)
sigma=math.sqrt(np.var(logreturn))

#annualized empirical mean and standard deviation
anu_mu= (250/(N-1))*sum(logreturn)
anu_sigma=np.sqrt((250/(N-1))*sum((logreturn-mu/250)**2))
print(anu_mu,anu_sigma)

#c
#stimulate time series data
logreturn_sim= np.random.normal(mu,sigma,N-1)
#plotting 2 time series data
plt.plot(logreturn,'r')
plt.plot(logreturn_sim,'b')
plt.title('Daily log returns')
plt.xlabel('Trading Day')
plt.ylabel('DAX log returns')
plt.show()

#d
#The simulated data is well-matched with the original data. However, it can not simulate the extreme data case, which means it can not simualate the extreme gain or extreme loss in the original return data