
BINARIES = testGauss testCrystalBall testDoubleCrystalBall testExponential testResos

BDIR = bin
HDIR = ../include
CDIR = src
ODIR = obj
LDIR = ../lib
DDIR = dep

#define MPI_ON
#1
#endef

HDRDIRS = $(HDIR)
LIBDIRS = $(LDIR)
LIBLIST = minuit cfit cfitmodels


#-------------------------------------------------------------------
# There should be no need to change much below this point.
#-------------------------------------------------------------------

CXX     = g++
CXXL    = g++

ifdef MPI_ON
# Directory of the MPI implementation.
MPIDIR	?= /afs/slac.stanford.edu/package/OpenMPI
CXXL     = $(MPIDIR)/bin/mpic++
HDRDIRS += $(MPIDIR)/include $(MPIDIR)/include/openmpi
LIBDIRS += $(MPIDIR)/lib
LIBLIST += mpi
endif

HDRSTR = $(foreach dir,$(HDRDIRS),-I $(dir))
LIBSTR = $(foreach dir,$(LIBDIRS),-L $(dir)) $(foreach lib,$(LIBLIST),-l$(lib))

CFLAGS  = -std=c++0x -g -O -Wall -fPIC $(HDRSTR)
DFLAGS  = -std=c++0x
LFLAGS  = -std=c++0x -g -O -Wall -fPIC $(LIBSTR)

ifdef MPI_ON
CFLAGS  += -DMPI_ON
LFLAGS  += -DMPI_ON
endif

RM      = rm -rf
LN      = ln -fs

ELVES  = $(foreach bin,$(BINARIES),$(BDIR)/$(bin))
OFILES = $(foreach bin,$(BINARIES),$(ODIR)/$(bin).o)
DFILES = $(foreach bin,$(BINARIES),$(DDIR)/$(bin).d)

all:    elves
elves:  $(ELVES)
.PHONY: tidy sweep clean

# Link rule.
$(ELVES): $(BDIR)/%: $(ODIR)/%.o $(MAKEFILE_LIST)
	@ mkdir -p $(dir $@)
	$(CXXL) -o $@ $< $(LFLAGS)
#	strip $@

# Rules for dependency files.
$(DFILES): $(DDIR)/%.d: $(CDIR)/%.cc $(MAKEFILE_LIST)
	@ mkdir -p $(dir $@)
	@ echo $(CXX) -MM -MF $@ -MT $$(echo $@ | sed "s/^$(DDIR)/$(ODIR)/g;s/\.d$$/\.o/g") $< -I $(HDIR)
	@ $(CXX) -MM -MF $@ -MT $$(echo $@ | sed "s/^$(DDIR)/$(ODIR)/g;s/\.d$$/\.o/g") $< $(DFLAGS) -I $(HDIR)

$(OFILES): $(ODIR)/%.o: $(CDIR)/%.cc $(DDIR)/%.d $(MAKEFILE_LIST)
	@ mkdir -p $(dir $@)
	$(CXX) -o $@ -c $< $(CFLAGS)

# Include all existent dependency files.
include $(wildcard $(DFILES))

tidy:

# This rule cleans up any backups or links to libraries only needed
#    at compilation time.
sweep:
	@ $(RM) $(ELF).out *~ `find -name '*~'` `find -name '\#*'` `find -name '\#*'`
	@ if [ $(ODIR) != "." ] && [ $(ODIR) != ".." ];	\
	then						\
		$(RM) $(ODIR);				\
	else						\
		$(RM) *.o;				\
	fi
	@ if [ $(DDIR) != "." ] && [ $(DDIR) != ".." ];	\
	then						\
		$(RM) $(DDIR);				\
	else						\
		$(RM) *.d;				\
	fi

clean: tidy sweep
