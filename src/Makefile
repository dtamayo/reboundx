ifndef REB_DIR
$(error Must set environment variable for the path to your rebound directory, e.g., export REB_DIR=/Users/dtamayo/rebound.  See https://github.com/dtamayo/reboundx/blob/master/README.rst)
endif

include $(REB_DIR)/src/Makefile.defs
OPT+= -fPIC -DLIBREBOUNDX

SOURCES=reboundx.c rebxtools.c modify_orbits_direct.c modify_orbits_forces.c gr.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) libreboundx.so
	
%.o: %.c $(HEADERS)
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
	@-rm -f *.o
	
