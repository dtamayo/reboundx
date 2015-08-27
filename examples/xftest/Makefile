export OPENGL=1

ifndef REB_DIR
$(error Must set environment variable for the path to your rebound directory, e.g., export REB_DIR=/Users/dtamayo/rebound.  See https://github.com/dtamayo/reboundx/blob/master/README.rst)
endif

include $(REB_DIR)/src/Makefile.defs

all: librebound libreboundx
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -g -I../../src/ -I$(REB_DIR)/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lreboundx -lrebound $(LIB) -o rebound
	@echo ""
	@echo "Problem file compiled successfully."

librebound:
	@echo "Compiling shared library libreboundx.so ..."
	$(MAKE) -C $(REB_DIR)/src/
	@echo "Creating link for shared library librebound.so ..."
	@-rm -f librebound.so
	@ln -s $(REB_DIR)/src/librebound.so .

libreboundx: 
	@echo "Compiling shared library libreboundx.so ..."
	$(MAKE) -C ../../src/
	@-rm -f libreboundx.so
	@ln -s ../../src/libreboundx.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C $(REB_DIR)/src/ clean
	@-rm -f libreboundx.so
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
