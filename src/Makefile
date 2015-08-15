OPT+= -fPIC -DLIBREBOUNDXF

SOURCES=libreboundxf.c xftools.c modify_elements_direct.c modify_elements_forces.c gr.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) libreboundxf.so 

%.o: %.c $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CC) -g -c $(OPT) $(PREDEF) -o $@ $<

libreboundxf.so: $(OBJECTS)
	@echo ""        
	@echo "Linking shared library $@ ..."
	$(CC) $(OPT) -shared  $(OBJECTS) $(LIB) -o $@ 
	
	@echo ""        
	@echo "The shared library $< has been created successfully."
	
	
clean:
	@echo "Removing object files *.o ..."
	@-rm -f *.o
	@echo "Removing shared library libreboundxf.so ..."
	@-rm -f *.so
	
