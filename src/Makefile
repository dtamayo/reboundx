ifndef REB_DIR
	
ifneq ($(wildcard ../../rebound/.*),) # Check for REBOUND in default location
REB_DIR=../../rebound
endif

ifneq ($(wildcard ../../../rebound/.*),) # Check for REBOUNDx being inside the REBOUND directory
REB_DIR=../../
endif

endif

ifndef REB_DIR # REBOUND is not in default location and REB_DIR is not set
    $(error REBOUNDx not in the same directory as REBOUND.  To use a custom location, you Must set the REB_DIR environment variable for the path to your rebound directory, e.g., export REB_DIR=/Users/dtamayo/rebound.  See reboundx.readthedocs.org)
endif

include $(REB_DIR)/src/Makefile.defs
OPT+= -fPIC -DLIBREBOUNDX

SOURCES=core.c gr.c gr_full.c gr_potential.c modify_mass.c modify_orbits_direct.c modify_orbits_forces.c radiation_forces.c rebxtools.c 
OBJECTS=$(SOURCES:.c=.o)
HEADERS=rebxtools.h reboundx.h modify_orbits_forces.h modify_orbits_direct.h core.h

all: $(SOURCES) libreboundx.so
	
%.o: %.c $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CC) -c $(OPT) $(PREDEF) -I$(REB_DIR)/src -o $@ $<

libreboundx.so: $(OBJECTS)
	@echo ""        
	@echo "Linking shared library $@ ..."
	$(CC) $(OPT) -shared $(OBJECTS) $(LIB) -lrebound -L$(REB_DIR)/src -o $@ 
	
	@echo ""        
	@echo "The shared library $@ has been created successfully."
	
	
clean:
	@echo "Cleaning up shared library libreboundx.so ..."
	@-rm -f libreboundx.so
	@-rm -f librebound.so
	@-rm -f *.o
	
