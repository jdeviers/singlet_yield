#!/bin/bash

# -- Vars --
ScriptLoc=$(readlink -f "$0")
FC="gfortran"
FCFLAGS=" -fbacktrace -Wall -Wno-tabs -fcheck=all "
MPFLAGS=" -fopenmp "

# -- Fcts --
compile() {
  ${FC} -c mod_rwfile.f90 random_sampling.f90 ${FCFLAGS};
  ${FC} -o sy_rs rs_caller.f90 mod_rwfile.o random_sampling.o ${FCFLAGS} ${MPFLAGS}
}

clean() {
  rm *.o *.mod sy_rs
}

# -- Script --
read -p "Press 1 to compile rs_caller, 2 to clean it, 3 to do both and run the program: " y

if [[ "${y}" == "1" ]]; then
  compile
elif [[ "${y}" == "2" ]]; then
  clean
elif [[ "${y}" == "3" ]]; then
  clean
  compile
  ./sy_rs
else
  echo "Wrong input."
  exec "${Scriptloc}"
fi

