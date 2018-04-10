#! /usr/bin/env python
## Soil-Plant-Atmospheric Continuum calculation emphasizing on surface energy balance
## Developed initially by Ardiansyah (UNSOED), http://ardiansyah.net, email : ardi.plj@gmail.com , cc : ard@unsoed.ac.id

import os as os ##operating system function
import numpy as np ##numerical python library
import matplotlib.pyplot as plt
import scipy as sc ##scientific python library
import math as mth ##matemathical operation library

def thomasAlgorithm(i1, A, B, C, D): ##tridiagonal matrix solution
    m = n-2 ## 15 nodes to evaluate, n = 17, including 0 and 17th
    for i in xrange(i1, m):
        C[i] = C[i]/B[i] ## update C 
        D[i] = D[i]/B[i] ## update D
        B[i+1] = B[i+1] - A[i+1]*C[i] ## update B
        D[i+1] = D[i+1] - A[i+1]*D[i] ## update D

    D[m] = D[m]/B[m]
    for i in xrange(m-1, i1-1, -1):
        D[i] = D[i] - C[i] * D[i+1]

    return D

def fourier(delta, j, z, a0, a, b): ##create fourier function for air temperature and air humidity
    i = np.arange(1, 5) ##create i array for sigma 
    w = 2 * np.pi / 24
    sigma = np.exp(-z*(i*w/(2*delta))**0.5) * (a * np.cos(i*w*time[j] - z*(i*w/(2*delta))**0.5)) + (b * np.sin(i*w*time[j] - z*(i*w/(2*delta))**0.5))
    fourier = a0 + sum(sigma)
    return fourier
