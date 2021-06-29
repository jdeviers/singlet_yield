PROG = $(wildcard prog_*.f90)
MODS = $(wildcard mod_*.f90)
OBJS = $(patsubst %.f90,%.o,$(MODS))

PROGRAM = sy_n

FC      = gfortran
FCFLAGS = -fbacktrace -Wall -Wtabs
LPFLAGS = -llapack


# Actions for make
default: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC) -o $@ $(PROG) $^ $(FCFLAGS) $(LPFLAGS)

$(OBJS): %.o : %.f90 
	$(FC) -c $< $(FCFLAGS)


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
