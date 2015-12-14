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
OPT+= -fPIC -DLIBREBOUNDX -ferror-limit=50

SOURCES=gr.c modify_orbits.c radiation_forces.c reboundx.c rebxtools.c 
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) libreboundx.so
	
%.o: %.c $(HEADERS)
	@echo $(OPT)
	@echo "Compiling source file $< ..."
	$(CC) -g -c $(OPT) $(PREDEF) -o $@ $<

libreboundx.so: $(OBJECTS)
	@echo ""        
	@echo "Linking shared library $@ ..."
	$(CC) $(OPT) -shared  $(OBJECTS) $(LIB) -o $@ 
	
	@echo ""        
	@echo "The shared library $< has been created successfully."
	
	
clean:
	@echo "Cleaning up shared library libreboundx.so ..."
	@-rm -f libreboundx.so
	@-rm -f librebound.so
	@-rm -f *.o
	
