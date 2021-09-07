#!/usr/bin/env python3

import numpy as np
import eYield_f90_staticbounds as f90_1
import eYield_f90_dynakernel as f90_2
from numpy import asfortranarray

# System setup
Sx = asfortranarray([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
Sy = asfortranarray([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
Sz = asfortranarray([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])

tmp1 = asfortranarray([Sx,Sy,Sz])
tmp2 = asfortranarray([Sz,Sy,Sx])

lambda1 = asfortranarray([1.,2.,3.])
lambda2 = asfortranarray([4.,8.,12.])

d1 = tmp1.shape[0]
d2 = tmp2.shape[0]

Sxyz1 = tmp1.copy()
Sxyz2 = tmp2.copy()

k = 1.

print("\nInput spin matrices:\n")
print(Sxyz1)
print("\n------------STATIC BOUNDS DOCTYPE:----------------\n")
print(f90_1.evalyield_offdiag2p.__doc__)
print("\n---------KERNEL DYNAMIC BOUNDS DOCTYPE:-----------\n")
print(f90_2.evalyield_offdiag2p.__doc__)


# Testing f90 routine

f90_1.evalyield_offdiag2p(k, Sxyz1, lambda1, Sxyz2, lambda2)
print("\nSTATIC BOUNDS RESULTS:\t\t",f90_1.evalyield_offdiag2p(k, Sxyz1, lambda1, Sxyz2, lambda2) )
print("KERNEL DYNAMIC BOUNDS RESULTS:\t",f90_2.evalyield_offdiag2p(k, Sxyz1, lambda1, Sxyz2, lambda2),"\n" )
