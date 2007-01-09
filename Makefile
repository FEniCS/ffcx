PREFIX = /usr/local


default: all 

ufc.pc:
	@echo "Name: UFC ">> ufc.pc 
	@echo "Version: 1.0 ">> ufc.pc 
	@echo "Description: Unified Form-assembly Code" >> ufc.pc  
	@echo "Cflags: -I$(PREFIX)/include ">> ufc.pc  

	


all:
	@echo "Nothing to build, type 'make install' to install."

install: ufc.pc 
	mkdir -p $(PREFIX)/include
	cp src/ufc/ufc.h $(PREFIX)/include
	mkdir -p $(PREFIX)/lib/pkgconfig
	cp ufc.pc $(PREFIX)/lib/pkgconfig/. 
	cd src/utils/python &&  python setup.py install  --prefix=$(PREFIX)

