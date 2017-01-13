export OPENGL=1

ifndef REB_DIR
ifneq ($(wildcard ../../../rebound/.*),) # Check for REBOUND in default location
REB_DIR=../../../rebound
endif
ifneq ($(wildcard ../../../../rebound/.*),) # Check for REBOUNDx being inside REBOUND directory
REB_DIR=../../../
endif
endif
ifndef REB_DIR # REBOUND is not in default location and REB_DIR is not set
    $(error REBOUNDx not in the same directory as REBOUND.  To use a custom location, you Must set the REB_DIR environment variable for the path to your rebound directory, e.g., export REB_DIR=/Users/dtamayo/rebound.  See reboundx.readthedocs.org)
endif
PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)

include $(REB_DIR)/src/Makefile.defs

REBX_DIR=../../

all: librebound.so libreboundx.so
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I$(REBX_DIR)/src/ -I$(REB_DIR)/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lreboundx -lrebound $(LIB) -o rebound
	@echo ""
	@echo "Problem file compiled successfully."

librebound.so:
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C $(REB_DIR)/src/
	@echo "Creating link for shared library librebound.so ..."
	@-rm -f librebound.so
	@ln -s $(REB_DIR)/src/librebound.so .

libreboundx.so: 
	@echo "Compiling shared library libreboundx.so ..."
	$(MAKE) -C $(REBX_DIR)/src/
	@-rm -f libreboundx.so
	@ln -s $(REBX_DIR)/src/libreboundx.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C $(REB_DIR)/src/ clean
	@echo "Cleaning up shared library libreboundx.so ..."
	@-rm -f libreboundx.so
	$(MAKE) -C $(REBX_DIR)/src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
