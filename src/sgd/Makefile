# General makefile

MAKE=make
SHELL=/bin/sh

SUBDIRS=crf


world: all
	@echo "================================================"
	@echo "CONGRATULATIONS: The compilation was successful."
	@echo "To know what to do next, check the README file."
	@echo "================================================"


all clean:
	@for n in ${SUBDIRS} ; \
	  do ( cd $$n && ${MAKE} ${@}) || exit ; done

all: 

.PHONY: world all depend check
