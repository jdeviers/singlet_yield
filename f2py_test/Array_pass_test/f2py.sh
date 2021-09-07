#!/bin/bash

# f2py
# -c                          -> compile 
# --fcompiler='gnu95'         -> compiler name
# -m eYield_f90               -> python module name
# eYield_f90_standalone.f90   -> Fortran source file

f2py3 func.f90 -h func.pyf
f2py3 --fcompiler='gnu95' -c func.f90 -m arr -DF2PY_REPORT_ON_ARRAY_COPY=1
