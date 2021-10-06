#!/usr/bin/env python3

import numpy as np
from numpy import asfortranarray

import py3_f90_evalYield
from py3_f90_evalYield import evalyield as f90_mod


def evalYield_full_py(k, Sxyz1, lambda1, Sxyz2, lambda2):
    # Sxyz = [Sx, Sy, Sz]
    d1 = Sxyz1.shape[-1]
    d2 = Sxyz2.shape[-1]
    z = d1 * d2 // 4
    v = 0.0
    counter = 0
    k2 = k*k
    for a1 in range(d1):
        for a2 in range(d1):
            dl1 = lambda1[a1] - lambda1[a2]
            sA = Sxyz1[:,a1,a2]
            for b1 in range(d2):
                for b2 in range(d2):
                    dl2 = lambda2[b1] - lambda2[b2]
                    sB = Sxyz2[:,b1,b2]
                    counter += 1
                    v += np.abs(sA[0]*sB[0] + sA[1]*sB[1] + sA[2]*sB[2])**2 / (k2 + (dl1 + dl2)**2)
    v *= k2/z
#    print(counter)
    return 1/4 + v

def evalYield_diag_py(k, Sxyz1, lambda1, Sxyz2, lambda2):
    # Sxyz = [Sx, Sy, Sz]
    # Sxyz1 = np.transpose(Sxyz1, (1,2,0))
    # Sxyz2 = np.transpose(Sxyz2, (1,2,0))
    d1 = Sxyz1.shape[0]
    d2 = Sxyz2.shape[0]
    z = d1 * d2 // 4
    v = 0.0
    for a1 in range(d1):
        sA = Sxyz1[a1,a1,:]
        for b1 in range(d2):
            sB = Sxyz2[b1,b1,:]
            #print(np.abs(sA[0]*sB[0] + sA[1]*sB[1] + sA[2]*sB[2])**2)
            v += np.abs(sA[0]*sB[0] + sA[1]*sB[1] + sA[2]*sB[2])**2
    v /= z
    return 1/4 + v



def evalYield_offdiag2p_kernel_F_py(k2, a1, Sxyz1_a1, lambda1, Sxyz2, lambda2):
    d1 = Sxyz1_a1.shape[-1]
    d2 = Sxyz2.shape[-1]
    lambda1_a1 = lambda1[a1]
    y = 0.0
    for b1 in range(d2):
        Sxyz2_b1 = Sxyz2[:,:,b1].copy()
        a2 = a1
        b2 = b1
        sAx = Sxyz1_a1[0,a2]
        sAy = Sxyz1_a1[1,a2]
        sAz = Sxyz1_a1[2,a2]
        dl1 = lambda1_a1 - lambda1[a2]
        while True:
            b2 += 1
            if b2 == d2:
                b2 = 0
                a2 += 1
                if a2 == d1:
                    break
                sAx = Sxyz1_a1[0,a2]
                sAy = Sxyz1_a1[1,a2]
                sAz = Sxyz1_a1[2,a2]
                dl1 = lambda1_a1 - lambda1[a2]
            sBx = Sxyz2_b1[0,b2]
            sBy = Sxyz2_b1[1,b2]
            sBz = Sxyz2_b1[2,b2]
            dl2 = lambda2[b1] - lambda2[b2]
            y += np.abs(sAx*sBx + sAy*sBy + sAz*sBz)**2 / (k2 + (dl1 + dl2)**2)
    return y


def evalYield_offdiag2p_py(k, Sxyz1, lambda1, Sxyz2, lambda2):
    d1 = Sxyz1.shape[-1]
    d2 = Sxyz2.shape[-1]
    z = (d1 * d2) // 4
    v = 0.0
    k2 = k*k
    for a1 in range(d1):
        Sxyz1_a1 = Sxyz1[:,:,a1].copy()
        v += evalYield_offdiag2p_kernel_F_py(k2, int(a1), Sxyz1_a1, lambda1, Sxyz2, lambda2)
    v *= 2*k2/z
    return v

# System setup
#Sx = asfortranarray([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
#Sy = asfortranarray([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
#Sz = asfortranarray([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])

N = 10

Sx = asfortranarray(np.random.rand(N,N))
symm_x = np.tril(Sx) + np.tril(Sx, -1).T

Sy = asfortranarray(np.random.rand(N,N))
symm_y = np.tril(Sy) + np.tril(Sy, -1).T

Sz = asfortranarray(np.random.rand(N,N))
symm_z = np.tril(Sz) + np.tril(Sz, -1).T

tmp1 = asfortranarray([Sx,Sy,Sz])
tmp2 = asfortranarray([Sz,Sy,Sx])

lambda1 = asfortranarray(np.random.rand(N))
lambda2 = asfortranarray(np.random.rand(N))

#lambda1 = asfortranarray([1.,2.,3.])
#lambda2 = asfortranarray([4.,8.,12.])

d1 = tmp1.shape[0]
d2 = tmp2.shape[0]

Sxyz1 = tmp1.copy()
Sxyz2 = tmp2.copy()

Sxyz1T = np.transpose(Sxyz1, (1,2,0)).copy()
Sxyz2T = np.transpose(Sxyz2, (1,2,0)).copy()

d = d1*d2
RS_scaling = (d**2 - d)

k = 1.

print("\nInput spin matrices:\n")
print(Sxyz1)

print("\n---------evalYield MODULE DOCTYPE:-----------\n")
print(f90_mod.__doc__)


# Testing f90 routine

print("Sxyzi shape:",Sxyz1.shape)
print("SxyziT shape:",Sxyz1T.shape)



print("\nPYTHON CODES:")
print("\nCodes that use Sxyzi: evalYield_full_py, evalYield_offdiag2p")
print("Codes that use SxyziT: evalYield_diag_py")
print("\nFULL YIELD:\t\t\t\t\t",evalYield_full_py(k, Sxyz1, lambda1, Sxyz2, lambda2) )
print("DIAG CONTRIB: FULL SAMPLING:\t\t\t",evalYield_diag_py(k, Sxyz1T, lambda1, Sxyz2T, lambda2) )
print("OFFDIAG CONTRIB: FULL SAMPLING (SERIAL):\t",evalYield_offdiag2p_py(k, Sxyz1, lambda1, Sxyz2, lambda2) )

print("\nFORTRAN CODES:")
print("\nCodes that use Sxyzi: evalyield_full, evalYield_offdiag2p*")
print("Codes that use SxyziT: evalyield_diag,evalyield_offdiag_random")
print("\nFULL YIELD:\t\t\t\t\t",f90_mod.evalyield_full(k, Sxyz1, lambda1, Sxyz2, lambda2) )
print("DIAG CONTRIB: FULL SAMPLING:\t\t\t",f90_mod.evalyield_diag(k, Sxyz1T, Sxyz2T) )
print("OFFDIAG CONTRIB: FULL SAMPLING (SERIAL):\t",f90_mod.evalyield_offdiag2p_serial(k, Sxyz1, lambda1, Sxyz2, lambda2) )
print("OFFDIAG CONTRIB: FULL SAMPLING (PARALLEL):\t",f90_mod.evalyield_offdiag2p(k, Sxyz1, lambda1, Sxyz2, lambda2),"\n" )

print("For such a small system, offdiag random is very sensitive to nr_draws. This algo however reproduces the full sampling methods very well for larger systems (e.g d1,d2 = 250).")
print("OFFDIAG CONTRIB: RANDOM SAMPLING:", RS_scaling * f90_mod.evalyield_offdiag_random(k, Sxyz1T, lambda1, Sxyz2T, lambda2, 5))


