PREFIX = /usr/local

all:
	@echo "Nothing to build, type 'make install' to install."

install:
	[ -d $(PREFIX)/include ] || mkdir -p $(PREFIX)/include
	cp src/ufc/ufc.h $(PREFIX)/include
