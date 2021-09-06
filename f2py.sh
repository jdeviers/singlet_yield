#!/bin/bash


# f2py
# -c                          -> compile 
# --fcompiler='gnu95'         -> compiler name
# -m eYield_f90               -> python module name
# eYield_f90_standalone.f90   -> Fortran source file


# All array bounds are declared at compile time
f2py --fcompiler='gnu95' --f90flags="-fopenmp" -c eYield_f90_staticbounds.f90 -m eYield_f90_staticbounds

# Array bounds delared at: compile time for the function, runtime for the kernel
f2py --fcompiler='gnu95' --f90flags="-fopenmp" -c eYield_f90.f90              -m eYield_f90
