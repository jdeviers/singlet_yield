#!/bin/bash


# f2py
# -c                          -> compile 
# --fcompiler='gnu95'         -> compiler name
# -m eYield_f90               -> python module name
# eYield_f90_standalone.f90   -> Fortran source file


# Random sampling
f2py3 -h interface_eY.pyf -m py3_f90_evalYield eYield_f90mod.f90 --overwrite-signature
f2py3 --fcompiler='gnu95' --f90flags="-fopenmp" interface_eY.pyf -c eYield_f90mod.f90
