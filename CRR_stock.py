# Tien Dat, Truong
# Group 14
# C-Exercise 03

import math
import numpy as np
import matplotlib.pyplot as plt

def CRR_stock(S_0, r, sigma, T, M):
    # estimate u, d and q
    delta_t = T / M
    alpha = math.exp(r * delta_t)
    beta = 1 / 2 * (1 / alpha + alpha * math.exp(math.pow(sigma, 2) * delta_t))
    u = beta + math.sqrt(math.pow(beta, 2) - 1)
    d = 1 / u
    q = (math.exp(r * delta_t) - d) / (u - d)

    # create return matrix S
    S = np.empty((M + 1, M + 1))

    # calculate value and fill in the matrix
    for i in range(1, M + 2, 1):
        for j in range(1, i + 1, 1):
            S[j - 1, i - 1] = S_0 * math.pow(u, j - 1) * math.pow(d, i - j)
    return S[M,M]

# test parameters
print(CRR_stock(100,0.05,0.3,1,500))
