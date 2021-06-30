PROG = $(wildcard prog_*.f90)
MODS = $(wildcard mod_*.f90)
OBJS = $(patsubst %.f90,%.o,$(MODS))

PROGRAM = sy_n

FC      = gfortran
FCFLAGS = -fbacktrace -Wall -Wtabs -fcheck=all
LPFLAGS = -llapack
MPFLAGS = -fopenmp


# Actions for make
default: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) -o $@ $(PROG) $^ $(FCFLAGS) $(LPFLAGS) $(MPFLAGS)

$(OBJS): %.o : %.f90 
	$(FC) -c $< $(FCFLAGS) $(MPFLAGS)


# Actions for make debug
debug:
	@echo $(PROG)
	@echo $(MODS)
	@echo $(OBJS)
	@echo $(FCFLAGS)
	@echo $(LPFLAGS)


# Actions for make clean
clean:
	rm $(PROGRAM) $(OBJS) $(patsubst %.f90,%.mod,$(MODS))

.PHONY = default debug clean
