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

plot() {
  gnuplot plot*.plt
  evince *.pdf
}

# -- Script --
echo -e "Clean, compile and run the random_sampling program. To use, press:"
echo -e "  * 1 to compile rs_caller;\n  * 2 to clean rs_caller files;\n  * 3 to clean then compile and run;\n  * 4 to plot results."
read -p "Type here: " y

case $y in
  1) compile
     ;;
  2) clean
     ;;
  3) clean 
     compile
     ./sy_rs
     ;;
  4) plot
     ;;
  *) echo "Wrong input" 
     exec "${Scriptloc}" 
     ;;
esac

