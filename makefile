#PROG = $(wildcard prog_*.f90)

PROG_comp = prog_sy.f90
PROG_time = prog_sy_timing.f90
PROG_para = prog_sy_parallel.f90

MODS = $(wildcard mod_*.f90)
OBJS = $(patsubst %.f90,%.o,$(MODS))

PROGRAM_comp = sy_c
PROGRAM_time = sy_t
PROGRAM_para = sy_p

FC      = gfortran
FCFLAGS = -fbacktrace -Wall -Wtabs -fcheck=all
MPFLAGS = -fopenmp


# Actions for make
default: $(PROGRAM_comp)

$(PROGRAM_comp): $(OBJS)
	$(FC) -o $@ $(PROG_comp) $^ $(FCFLAGS) $(MPFLAGS)

$(OBJS): %.o : %.f90 
	$(FC) -c $< $(FCFLAGS) $(MPFLAGS)


# Actions for make timing
timing: $(PROGRAM_time)

$(PROGRAM_time): $(OBJS)
	$(FC) -o $@ $(PROG_time) $^ $(FCFLAGS) $(MPFLAGS)


# Actions for make parallel
parallel: $(PROGRAM_para)

$(PROGRAM_para): $(OBJS)
	$(FC) -o $@ $(PROG_para) $^ $(FCFLAGS) $(MPFLAGS)


# Actions for make debug
debug:
	@echo $(PROG_comp)
	@echo $(PROG_time)
	@echo $(PROG_para)
	@echo $(MODS)
	@echo $(OBJS)
	@echo $(FCFLAGS)
	@echo $(LPFLAGS)


# Actions for make clean
clean:
	rm sy_* $(OBJS) $(patsubst %.f90,%.mod,$(MODS))

.PHONY = default timing parallel debug clean
