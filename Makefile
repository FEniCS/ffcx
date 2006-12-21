PREFIX = /usr/local

all:
	@echo "Nothing to build, type 'make install' to install."

install:
	cp src/ufc/ufc.h $(PREFIX)/include
