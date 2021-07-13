#!/bin/bash

# Doesn't work yet 

# f2py
# -c                          -> compile 
# --fcompiler='gnu95'         -> compiler name
# -m eYield_f90               -> python module name
# eYield_f90_standalone.f90   -> Fortran source file

f2py -c --fcompiler='gnu95' -m eYield_f90 eYield_f90_standalone.f90 --f90flags=fopenmp
