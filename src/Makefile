OPT+= -fPIC -DLIBREBOUNDX

SOURCES=reboundx.c rebxtools.c modify_orbits_direct.c modify_orbits_forces.c gr.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) libreboundx.so 

%.o: %.c $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CC) -g -c $(OPT) -I$(REB_DIR)/src $(PREDEF) -o $@ $<

libreboundx.so: $(OBJECTS)
	@echo ""        
	@echo "Linking shared library $@ ..."
	$(CC) $(OPT) -shared  $(OBJECTS) $(LIB) -o $@ 
	
	@echo ""        
	@echo "The shared library $< has been created successfully."
	
	
clean:
	@echo "Removing object files *.o ..."
	@-rm -f *.o
	@echo "Removing shared library libreboundx.so ..."
	@-rm -f *.so
	
