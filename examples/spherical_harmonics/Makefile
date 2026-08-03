export OPENGL=0

REB_DIR = ../../../harmonics/Code/rebound

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

libreboundx.so: librebound.so 
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
