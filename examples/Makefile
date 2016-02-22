
ifndef REB_DIR
ifneq ($(wildcard ../../rebound/.*),) # Check for REBOUND in default location
REB_DIR=../../rebound
endif

ifneq ($(wildcard ../../../rebound/.*),) # Check for REBOUNDx being inside REBOUND directory
REB_DIR=../../
endif
endif

ifndef REB_DIR # REBOUND is not in default location and REB_DIR is not set
    $(error REBOUNDx not in the same directory as REBOUND.  To use a custom location, you Must set the REB_DIR environment variable for the path to your rebound directory, e.g., export REB_DIR=/Users/dtamayo/rebound.  See reboundx.readthedocs.org)
endif

include $(REB_DIR)/src/Makefile.defs

REBOUNDXDIR=../src/
DEBUGOPT = -std=c99 -Wpointer-arith -D_GNU_SOURCE -O0 -march=native -I/usr/local/include -Wall -Wno-deprecated -g
SOURCES=$(REBOUNDXDIR)core.c $(REBOUNDXDIR)gr.c $(REBOUNDXDIR)gr_full.c $(REBOUNDXDIR)gr_potential.c $(REBOUNDXDIR)modify_mass.c $(REBOUNDXDIR)modify_orbits_direct.c $(REBOUNDXDIR)modify_orbits_forces.c $(REBOUNDXDIR)radiation_forces.c $(REBOUNDXDIR)rebxtools.c 
SUBDIRS := $(wildcard */.)

all: librebound.so libreboundx.so
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I../src/ -I$(REB_DIR)/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lreboundx -lrebound $(LIB) -o rebound
	@echo ""
	@echo "Problem file compiled successfully."

eclipse: librebound.so
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I../src/ -Wl,-rpath,./ $(DEBUGOPT) $(PREDEF) problem.c $(SOURCES) -L. -L ../ -lrebound $(LIB) -o rebound
	@echo ""
	@echo "Problem file compiled successfully."
 
allexamples: $(SUBDIRS)
$(SUBDIRS):
	@echo "Trying to compile example $(subst /.,,$@):"
	$(MAKE) -s -C $@ clean 
	$(MAKE) -s -C $@ 
	@echo "\033[1A\033[55CSuccess."

.PHONY: allexamples $(SUBDIRS)

librebound.so:
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C $(REB_DIR)/src/
	@echo "Creating link for shared library librebound.so ..."
	@-rm -f librebound.so
	@ln -s $(REB_DIR)/src/librebound.so .

libreboundx.so: 
	@echo "Compiling shared library libreboundx.so ..."
	$(MAKE) -C ../src/
	@-rm -f libreboundx.so
	@ln -s ../src/libreboundx.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C $(REB_DIR)/src/ clean
	@-rm -f libreboundx.so
	@-rm -f problem.c
	$(MAKE) -C ../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
