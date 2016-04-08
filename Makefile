# Edit the following variables as needed
#HDF_INSTALL = $(HOME)/hdf5
#
CC = mpicc
# define any compile-time flags
CFLAGS = -Werror -Wall -g -pedantic -std=gnu99 -Wno-long-long
OPTIMIZATION = -fopenmp -O2 -vec-report2 -xAVX -Wno-strict-aliasing -fomit-frame-pointer
# CFLAGS = -Wall -g -pedantic -std=gnu99 -Wno-long-long
# OPTIMIZATION = -fopenmp -O2 -ffast-math -Wno-strict-aliasing -fomit-frame-pointer

INCLUDES = -I$(HDF5_INCL)
LFLAGS =
HDF5LIB = -L$(HDF5_ROOT)/lib -lhdf5
# LIBS = $(HDF5LIB) -ldl -lm -lgsl -lgslcblas
LIBS = $(HDF5LIB) -ldl -lm

# define the C source files
# SRCS_CHAOTICB = cbmpi.c diagnostics.c force_free.c quick_sort.c StepperBS.c \
# 	   tracking.c wlcs.c domain.c
SRCS = domain.c wlcs.c particle_info.c diagnostics.c emfields.c \
	   quick_sort.c tracking.c force_free.c data_io.c velocity_field.c \
	   interpolation.c magnetic_field.c bessel.c
SRCS_CHAOTICB = main.c $(SRCS)
SRCS_MAGNETIC = magnetic_ene.c $(SRCS)

# define the C object files
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS_CHAOTICB = $(SRCS_CHAOTICB:.c=.o)

OBJS_MAGNETIC = $(SRCS_MAGNETIC:.c=.o)

# define the executable file
CHAOTICB = test_particle
MAGNETIC_ENE = eneb

#
.PHONY: depend clean

all:	$(CHAOTICB) $(MAGNETIC_ENE)
	@echo  Programs are successfully compiled!

chaoticb:	$(CHAOTICB)
	@echo  $(CHAOTICB) are successfully compiled!

magnetic_ene:	$(MAGNETIC_ENE)
	@echo  $(MAGNETIC_ENE) is successfully compiled!

$(CHAOTICB): $(OBJS_CHAOTICB)
	$(CC) $(CFLAGS) $(OPTIMIZATION) $(INCLUDES) -o $(CHAOTICB) \
		$(OBJS_CHAOTICB) $(LFLAGS) $(LIBS)

$(MAGNETIC_ENE): $(OBJS_MAGNETIC)
	$(CC) $(CFLAGS) $(OPTIMIZATION) $(INCLUDES) -o $(MAGNETIC_ENE) \
		$(OBJS_MAGNETIC) $(LFLAGS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(OPTIMIZATION) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN) $(TRAJ)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
