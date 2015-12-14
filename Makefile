clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc clean

.PHONY: doc
doc: 
	cd doc/doxygen && doxygen
	$(MAKE) -C doc html
		
