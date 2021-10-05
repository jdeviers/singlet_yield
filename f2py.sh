#!/bin/bash


# f2py
# -c                          -> compile 
# --fcompiler='gnu95'         -> compiler name
# -m eYield_f90               -> python module name
# eYield_f90_standalone.f90   -> Fortran source file


# All array bounds are declared at compile time
f2py3 -h interface_SK.pyf -m py3_f90_staticKernel eYield_f90_staticKernel.f90 --overwrite-signature
f2py3 --fcompiler='gnu95' --f90flags="-fopenmp" interface_SK.pyf -c eYield_f90_staticKernel.f90

# Array bounds declared at: compile time for the function, runtime for the kernel
f2py3 -h interface_DK.pyf -m py3_f90_dynamicKernel eYield_f90_dynamicKernel.f90 --overwrite-signature
f2py3 --fcompiler='gnu95' --f90flags="-fopenmp" interface_DK.pyf -c eYield_f90_dynamicKernel.f90

# Random sampling
f2py3 -h interface_RS.pyf -m py3_f90_randomSampling eYield_f90_randomSampling.f90 --overwrite-signature
f2py3 --fcompiler='gnu95' --f90flags="-fopenmp" interface_RS.pyf -c eYield_f90_randomSampling.f90
