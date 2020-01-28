# project name (generate executable with this name)
TARGET   = test_particle
#
CC     = mpicc
LD     = mpicc
AR     = ar clq
RANLIB = ranlib
RM 	   = rm -f
CFLAGS = -Werror -Wall -g -pedantic -std=gnu99 -Wno-long-long -fopenmp
OPTIMIZATION = -O2 -qopt-report-phase=vec -fp-model precise -xAVX -Wno-strict-aliasing -fomit-frame-pointer
# OPTIMIZATION =

# Specify HDF5_ROOT
HDF5_ROOT = /usr/projects/hpcsoft/toss3/grizzly/hdf5/1.8.16_intel-18.0.2_openmpi-2.1.2/
HDF5_INC = $(HDF5_ROOT)/include
HDF5_LIB = -L$(HDF5_ROOT)/lib -lhdf5_hl -lhdf5
LIBS = $(HDF5_LIB) -ldl
LFLAGS =

CFLAGS += -I$(HDF5_INC)

# define the C source files
SRCDIR   = src
INCDIR	 = include
OBJDIR   = obj
BINDIR   = bin

CFLAGS += -I./$(INCDIR)

SOURCES  := $(wildcard $(SRCDIR)/*.c)
INCLUDES := $(wildcard $(INCDIR)/*.h)
OBJECTS  := $(addprefix $(OBJDIR)/,$(notdir $(SOURCES:.c=.o)))

all:	$(BINDIR)/$(TARGET)
	@echo  Programs are successfully compiled!

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(LD) $(CFLAGS) $(OPTIMIZATION) $(LFLAGS) $(LIBS) -o $@ $^
	@echo $(BINDIR)/$(TARGET) are successfully compiled!

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CFLAGS) $(OPTIMIZATION) -c $< -o $@

.PHONY: depend clean

clean:
	$(RM) $(OBJECTS) $(OBJDIR)/*.optrpt *~ $(BINDIR)/$(TARGET)

depend: $(SOURCES)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
