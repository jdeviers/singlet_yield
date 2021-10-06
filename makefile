PROG_comp = prog_sy.f90
PROG_time = prog_sy_timing.f90
PROG_para = prog_sy_parallel.f90
PROG_rndm = prog_sy_randomSampling.f90

D_SRC=./sources/
D_MOD=./modules/
D_BIN=./

MODS = $(wildcard mod_*.f90)
OBJS = $(patsubst %.f90,%.o,$(MODS))

PROGRAM_comp = sy_c
PROGRAM_time = sy_t
PROGRAM_para = sy_p
PROGRAM_rndm = sy_rs

FC      = gfortran
FCFLAGS = -fbacktrace -Wall -Wno-tabs -fcheck=all
MPFLAGS = -fopenmp


# Actions for make
default: $(PROGRAM_comp)

$(PROGRAM_comp): $(OBJS)
	$(FC) -o $@ $(PROG_comp) $^ $(FCFLAGS) $(MPFLAGS)

$(OBJS): %.o : %.f90 
	$(FC) -c $< $(FCFLAGS) $(MPFLAGS)


mod_sy_proc.o mod_random_sampling.o: mod_rwfile.o


# Actions for make timing
timing: $(PROGRAM_time)

$(PROGRAM_time): $(OBJS)
	$(FC) -o $@ $(PROG_time) $^ $(FCFLAGS) $(MPFLAGS)


# Actions for make parallel
parallel: $(PROGRAM_para)

$(PROGRAM_para): $(OBJS)
	$(FC) -o $@ $(PROG_para) $^ $(FCFLAGS) $(MPFLAGS)


# Actions for make random
random: $(PROGRAM_rndm)

$(PROGRAM_rndm): $(OBJS)
	$(FC) -o $@ $(PROG_rndm) $^ $(FCFLAGS) $(MPFLAGS)


# Actions for make debug
debug:
	@echo $(PROG_comp)
	@echo $(PROG_time)
	@echo $(PROG_para)
	@echo $(PROG_rndm)
	@echo $(MODS)
	@echo $(OBJS)
	@echo $(FCFLAGS)
	@echo $(LPFLAGS)


# Actions for make clean
clean:
	rm sy_* $(OBJS) $(patsubst %.f90,%.mod,$(MODS))

.PHONY = default timing parallel random debug clean
