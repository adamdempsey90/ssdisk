EXECUTABLES=../run_step_fixed ../run_step_bary
SOURCES=allocate.c average.c boundary.c cfl.c artificial_visc.c disk.c main.c input.c read_pars.c lambda.c output.c planets.c potential.c pressure.c source_step.c transport.c update.c vanleer.c viscosity.c stockholm.c advection.c ppa.c compute_vmed.c convolution.c
HEADER=evolve.h

CC=gcc

MATHLIB=-lm
WARNFLAG=-Wall -Wextra -std=c99
OPTFLAG=-O3

UNAME := $(shell echo $(USER))
OMPLIB=
OMPFLAG=
FFTWLIB=
FFTWINC=

ifeq ($(UNAME),amd616)

OMPLIB=-lgomp
#OMPLIB=-lrfftw_threads -lfftw_threads -lgomp
OMPFLAG=-fopenmp

#FFTWLIB=-L/software/FFTW/2.1.5-gcc/lib/  -lrfftw -lfftw
#FFTWINC=-I/software/FFTW/2.1.5-gcc/include/ -D_FFTW
endif

ifeq ($(UNAME),jupiter)
CC=gcc-6
OMPLIB=-lgomp
OMPFLAG=-fopenmp
endif

ifeq ($(UNAME),zeus)
CC=gcc-6
OMPLIB=-lgomp
OMPFLAG=-fopenmp
endif

LDFLAGS= $(FFTWLIB) $(OMPLIB) $(MATHLIB)
CFLAGS=-c $(OPTFLAG) $(WARNFLAG) $(OMPFLAG) $(FFTWINC)

FBIN=bin/fixed/
BBIN=bin/bary/

OBJECTS=$(SOURCES:.c=.o)
BOBJECTS=$(addprefix $(BBIN),$(OBJECTS))
FOBJECTS=$(addprefix $(FBIN),$(OBJECTS))

all: $(EXECUTABLES)

../run_step_fixed: $(FOBJECTS)
	$(CC) $(FOBJECTS) $(LDLIB) $(LDFLAGS) -o $@

../run_step_bary: $(BOBJECTS)
	$(CC) $(BOBJECTS) $(LDLIB) $(LDFLAGS) -o $@

$(FBIN)%.o: %.c $(HEADER)
	$(CC) $(INCLIB) $(CFLAGS) -DFIXEDPSYS $< -o $@

$(BBIN)%.o: %.c $(HEADER)
	$(CC) $(INCLIB) $(CFLAGS) $< -o $@

clean:
	rm $(BOBJECTS) $(FOBJECTS) $(EXECUTABLES)










