# Tien Dat, Truong
# Group 14
# C-Exercise 23

import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

def CRR_AmOption (S_0, T, M):
    r=0
    sigma= np.sqrt(2)
    # compute values of u, d and q
    delta_t = T / M
    alpha = math.exp(r * delta_t)
    beta = 1 / 2 * (1 / alpha + alpha * math.exp(math.pow(sigma, 2) * delta_t))
    u = beta + math.sqrt(math.pow(beta, 2) - 1)
    d = 1 / u
    q = (math.exp(r * delta_t) - d) / (u - d)

    # allocate matrix S
    S = np.empty((M + 1, M + 1))

    # fill matrix S with stock prices
    for i in range(1, M + 2, 1):
        for j in range(1, i + 1, 1):
            S[j - 1, i - 1] = S_0 * math.pow(u, j - 1) * math.pow(d, i - j)

    # g_ST will contain American option payoff
    g_ST = np.empty((M + 1, M + 1))

    for j in range(1, M + 2):
        if S[j - 1, M]<1:
            g_ST[j - 1, M] = S[j - 1, M]**(3/4)
        else:
            g_ST[j - 1, M] = 3 * np.sqrt(S[j - 1, M])+ S[j - 1, M]**(3/2)

    for i in range(M + 1, 1, -1):
        for j in range(1, i + 1, 1):
            if S[j - 2, i - 2]<1:
                g_ST[j - 2, i - 2] = np.maximum((S[j - 2, i-2]**(3/4)),(np.exp(-r * delta_t) * (q * g_ST[j - 1, i - 1] + (1 - q) * g_ST[j - 2, i - 1])))
            else:
                g_ST[j - 2, i - 2] = np.maximum((3 * np.sqrt(S[j - 2, i-2])+ S[j - 2, i-2]**(3/2)), (np.exp(-r * delta_t) * (q * g_ST[j - 1, i - 1] + (1 - q) * g_ST[j - 2, i - 1])))
    return g_ST[0,0]

print(CRR_AmOption (1, 1, 500))

def CRR_AmOptionDirect (S_0, T):
    t=0
    M=500
    r=0
    sigma= np.sqrt(2)
    # compute values of u, d and q
    delta_t = T / M
    alpha = math.exp(r * delta_t)
    beta = 1 / 2 * (1 / alpha + alpha * math.exp(math.pow(sigma, 2) * delta_t))
    u = beta + math.sqrt(math.pow(beta, 2) - 1)
    d = 1 / u
    q = (math.exp(r * delta_t) - d) / (u - d)

    # allocate matrix S
    S = np.empty((M + 1, M + 1))

    # fill matrix S with stock prices
    for i in range(1, M + 2, 1):
        for j in range(1, i + 1, 1):
            S[j - 1, i - 1] = S_0 * math.pow(u, j - 1) * math.pow(d, i - j)

    # Create and fill g_ST matrix
    g_ST = np.empty((M + 1, M + 1))

    for j in range(1, M + 2):
        if S[j - 1, M]<1:
            g_ST[j - 1, M] = S[j - 1, M]**(3/4)
        else:
            g_ST[j - 1, M] = 3 * np.sqrt(S[j - 1, M])+ S[j - 1, M]**(3/2)

    for i in range(M + 1, 1, -1):
        for j in range(1, i + 1, 1):
            if S[j - 2, i - 2]<1:
                g_ST[j - 2, i - 2] = np.maximum((S[j - 2, i-2]**(3/4)),(np.exp(-r * delta_t) * (q * g_ST[j - 1, i - 1] + (1 - q) * g_ST[j - 2, i - 1])))
            else:
                g_ST[j - 2, i - 2] = np.maximum((3 * np.sqrt(S[j - 2, i-2])+ S[j - 2, i-2]**(3/2)), (np.exp(-r * delta_t) * (q * g_ST[j - 1, i - 1] + (1 - q) * g_ST[j - 2, i - 1])))

    # Create and fill V2 matrix
    V2= np.empty((M + 1, M + 1))
    for j in range(1, M + 2):
        if S[j - 1, M]< np.exp(-(T-t)):
            V2[j - 1, M] = g_ST[j - 1, M]
        else:
            V2[j - 1, M] = np.exp((-1/4)*(T-t))*3*np.sqrt(S[j - 1, M]) + np.exp((3/4)*(T-t))*S[j - 1, M]**(3/2)

    for i in range(M + 1, 1, -1):
        for j in range(1, i + 1, 1):
            if S[j - 2, i - 2] < np.exp(-(T-t)):
                V2[j - 2, i - 2] = np.maximum(g_ST[j - 2, i - 2], (np.exp(-r * delta_t) * (q * V2[j - 1, i - 1] + (1 - q) * V2[j - 2, i - 1])))
            else:
                V2[j - 2, i - 2] = np.maximum((np.exp((-1/4)*(T-t))*3*np.sqrt(S[j - 2, i-2]) + np.exp((3/4)*(T-t))*S[j - 2, i-2]**(3/2)), (np.exp(-r * delta_t) * (q * V2[j - 1, i - 1] + (1 - q) * V2[j - 2, i - 1])))
    return V2[0,0]

print(CRR_AmOptionDirect (1, 1))