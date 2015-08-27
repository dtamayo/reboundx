export OPENGL=1

REB_DIR=/Users/rein/git/rebound/

include $(REB_DIR)src/Makefile.defs


all: librebound libreboundxf
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -g -I../../src/ -I$(REB_DIR)src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lreboundxf -lrebound $(LIB) -o rebound
	@echo ""
	@echo "Problem file compiled successfully."

librebound:
	@echo "Creating link for shared library librebound.so ..."
	@-rm -f librebound.so
	@ln -s $(REB_DIR)src/librebound.so .

libreboundxf: 
	@echo "Compiling shared library libreboundxf.so ..."
	$(MAKE) -C ../../../reboundxf/src/
	@-rm -f libreboundxf.so
	@ln -s ../../../reboundxf/src/libreboundxf.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	@-rm -f libreboundxf.so
	$(MAKE) -C ../../src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound


