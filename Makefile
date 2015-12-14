# This makefile compiles the shared libraries and updates the documentation.
#
ifndef REB_DIR
ifneq ($(wildcard ../rebound/.*),) # Check for REBOUND in default location
REB_DIR=../rebound
endif

ifneq ($(wildcard ../../rebound/.*),) # Check for REBOUNDx being inside REBOUND directory
REB_DIR=../
endif
endif

all: librebound.so libreboundx.so

librebound.so:
	@echo "Compiling shared library librebound.so ..."
	@echo $(REB_DIR)
	$(MAKE) -C $(REB_DIR)/src/
	@echo "Creating link for shared library librebound.so ..."
	@-rm -f librebound.so
	@ln -f -s $(REB_DIR)/src/librebound.so .

libreboundx.so: 
	@echo "Compiling shared library libreboundx.so ..."
	$(MAKE) -C src/
	@-rm -f libreboundx.so
	@ln -f -s src/libreboundx.so .

clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc clean

.PHONY: doc
doc: 
	cd doc/doxygen && doxygen
	$(MAKE) -C doc html
		
