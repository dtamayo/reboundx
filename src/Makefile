ifndef REB_DIR
$(error Must set environment variable for the path to your rebound directory, e.g., export REB_DIR=/Users/dtamayo/rebound.  See https://github.com/dtamayo/reboundx/blob/master/README.rst)
endif

OPT+= -fPIC -DLIBREBOUNDX

SOURCES=reboundx.c rebxtools.c modify_orbits_direct.c modify_orbits_forces.c gr.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) librebound.so libreboundx.so

%.o: %.c $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CC) -g -c $(OPT) -I$(REB_DIR)/src $(PREDEF) -o $@ $<

librebound.so:
	@echo "Compiling shared library libreboundx.so ..."
	$(MAKE) -C $(REB_DIR)/src/
	@echo "Creating link for shared library librebound.so ..."
	@-rm -f librebound.so
	@ln -s $(REB_DIR)/src/librebound.so .

libreboundx.so: $(OBJECTS)
	@echo ""        
	@echo "Linking shared library $@ ..."
	$(CC) $(OPT) -shared  $(OBJECTS) -L$(REB_DIR)/src/ -lrebound $(LIB) -o $@ 
	
	@echo ""        
	@echo "The shared library $< has been created successfully."
	
	
clean:
	@-rm -f librebound.so
	$(MAKE) -C $(REB_DIR)/src/ clean
	@echo "Cleaning up shared library libreboundx.so ..."
	@-rm -f libreboundx.so
	@-rm -f *.o
	
