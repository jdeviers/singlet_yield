import numpy as np

def evalYield(k, Sxyz1, lambda1, Sxyz2, lambda2):
    # Sxyz = [Sx, Sy, Sz]
    d1 = Sxyz1.shape[-1]
    d2 = Sxyz2.shape[-1]
    z = d1 * d2 // 4
    v = 0.0
    k2 = k*k
    for a1 in range(d1):
        for a2 in range(d1):
            dl1 = lambda1[a1] - lambda1[a2]
            sA = Sxyz1[:,a1,a2]
            for b1 in range(d2):
                for b2 in range(d2):
                    dl2 = lambda2[b1] - lambda2[b2]
                    sB = Sxyz2[:,b1,b2]
                    v += np.abs(sA[0]*sB[0] + sA[1]*sB[1] + sA[2]*sB[2])**2 / (k2 + (dl1 + dl2)**2)
    v *= k2/z
    return 1/4 + v
 
def evalYield_offdiag2p_kernel_F(k2, a1, Sxyz1_a1, lambda1, Sxyz2, lambda2):
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
#@jit(nopython=True,fastmath=True,parallel=True)
def evalYield_offdiag2p(k, Sxyz1, lambda1, Sxyz2, lambda2):
    d1 = Sxyz1.shape[-1]
    d2 = Sxyz2.shape[-1]
    z = (d1 * d2) // 4
    v = 0.0
    k2 = k*k
    for a1 in range(d1):
        Sxyz1_a1 = Sxyz1[:,:,a1].copy()
        v += evalYield_offdiag2p_kernel_F(k2, int(a1), Sxyz1_a1, lambda1, Sxyz2, lambda2)
    v *= 2*k2/z
    return v

Sx = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
Sy = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
Sz = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])

tmp1 = np.array([Sx,Sy,Sz])
tmp2 = np.array([Sz,Sy,Sx])

lambda1 = [1.,2.,3.]
lambda2 = [4.,8.,12.]

Sxyz1 = tmp1.copy()
Sxyz2 = tmp2.copy()

k = 1.

# It seems normal to obtain different results
print(evalYield(k, Sxyz1, lambda1, Sxyz2, lambda2))
print(evalYield_offdiag2p(k, Sxyz1, lambda1, Sxyz2, lambda2))

