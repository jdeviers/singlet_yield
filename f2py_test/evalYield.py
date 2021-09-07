#!/usr/bin/env python3

import numpy as np
#import eYield_f90_staticbounds as f90
import eYield_f90_dynakernel as f90

# System setup
Sx = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
Sy = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
Sz = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])

tmp1 = np.array([Sx,Sy,Sz])
tmp2 = np.array([Sz,Sy,Sx])

lambda1 = [1.,2.,3.]
lambda2 = [4.,8.,12.]

d1 = tmp1.shape[0]
d2 = tmp2.shape[0]

Sxyz1 = tmp1.copy()
Sxyz2 = tmp2.copy()

k = 1.

print(Sxyz1)
#help(f90)

# Testing f90 routine
f90.evalyield_offdiag2p(d1, d2, k, Sxyz1, lambda1, Sxyz2, lambda2)
#print(evalYield_offdiag2p(d1, d2, k, Sxyz1, lambda1, Sxyz2, lambda2))

