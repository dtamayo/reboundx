SUBDIRS := $(wildcard */.)

all: $(SUBDIRS)
$(SUBDIRS):
	@echo "************* Trying to compile example $(subst /.,,$@):"
ifndef SKIP_CLEAN
	$(MAKE) -s -C $@ clean 
endif
	$(MAKE) OPENGL=0 -s -C $@ 
	@echo "\033[1A\033[55CSuccess."

.PHONY: all $(SUBDIRS)
