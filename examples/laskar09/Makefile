PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
export OPENGL=0

all: 
	@-ln -fs ../$(PROBLEMDIR)/problem.c ../problem.c
	$(MAKE) -C ../
	@-cp ../rebound .
	@-cp ../librebound.so .
	@-cp ../libreboundx.so .

eclipse:
	@-ln -fs ../$(PROBLEMDIR)/problem.c ../problem.c
	$(MAKE) -C ../ eclipse
	@cp ../rebound .
	@cp ../librebound.so .

clean:
	$(MAKE) -C ../ clean
	@-rm -f librebound.so
	@-rm -f libreboundx.so
	@-rm -f rebound
