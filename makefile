
# define MPI_ON
# 1
# endef

all:
	$(MAKE) MPI_ON=$(MPI_ON) -f makefile.cfit
	$(MAKE) MPI_ON=$(MPI_ON) -f makefile.models

.PHONY: install tidy sweep clean

install:
	$(MAKE) -f makefile.cfit   install
	$(MAKE) -f makefile.models install

tidy:
	$(MAKE) -f makefile.cfit tidy

sweep:
	$(MAKE) -f makefile.cfit sweep

clean:
	$(MAKE) -f makefile.cfit clean
