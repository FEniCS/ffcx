PREFIX = /usr/local

default: all

pc:
	rm -f ufc-1.pc
	@echo "Name: UFC" >> ufc-1.pc
	@echo "Version: 1.0" >> ufc-1.pc
	@echo "Description: Unified Form-assembly Code" >> ufc-1.pc
	@echo "Cflags: -I$(PREFIX)/include" >> ufc-1.pc

all:
	@echo "Nothing to build, type 'make install' to install."

install: pc
	mkdir -p $(PREFIX)/include
	cp src/ufc/ufc.h $(PREFIX)/include
	mkdir -p $(PREFIX)/lib/pkgconfig
	cp ufc-1.pc $(PREFIX)/lib/pkgconfig
	cd src/utils/python && python setup.py install --prefix=$(PREFIX)
