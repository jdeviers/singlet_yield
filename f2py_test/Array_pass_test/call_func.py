#!/usr/bin/env python3

import numpy as np
import arr
from numpy import asfortranarray

#S = np.array([[1.,2.,3.], [4.,5.,6.], [7.,8.,9.]])

S = asfortranarray([[1.,2.,3.], [4.,5.,6.], [7.,8.,9.]])
print(S.flags.f_contiguous)
print(S)

a1 = S.shape[0]
a2 = S.shape[1]


#help(arr)
print(arr.passingarray.arrmanip.__doc__)

#arr.passingarray.arrmanip(a1,a2,S)
