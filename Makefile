clean:
	$(MAKE) -C src clean
	$(MAKE) -C doc clean
	@-rm -f *.so
	@pip uninstall -y reboundx
	@python setup.py clean --all
	@rm -rf reboundx.*

.PHONY: doc
doc: 
	cd doc/doxygen && doxygen
	$(MAKE) -C doc html
		
