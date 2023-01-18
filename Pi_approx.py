# Tien Dat, Truong
# Group 14
# C-Exercise 10

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

#Create function
def Pi_approx(N):
    inside=0
    #Creat null matrix for positon (x,y)
    x_inside = []
    y_inside = []
    #Fill in x, y matrix
    for _ in range(N):
        x = np.random.uniform(0.0, 1.0)
        y = np.random.uniform(0.0, 1.0)
        if x ** 2 + y ** 2 <= 1:
            inside += 1
            x_inside.append(x)
            y_inside.append(y)
    # Estimating pi
    pi = 4 * inside / N
    return pi,(x_inside,y_inside)
# Approximate pi with N=10,100,1000,10000
print(Pi_approx(10)[0],Pi_approx(100)[0],Pi_approx(1000)[0],Pi_approx(10000)[0])

# Plotting
a=(Pi_approx(10))
b=(Pi_approx(100))
c=(Pi_approx(1000))
d=(Pi_approx(10000))

plt.subplot(2,2,1)
plt.scatter(a[1][0],a[1][1])

plt.subplot(2,2,2)
plt.scatter(b[1][0],b[1][1])

plt.subplot(2,2,3)
plt.scatter(c[1][0],c[1][1])

plt.subplot(2,2,4)
plt.scatter(d[1][0],d[1][1])

plt.show()
