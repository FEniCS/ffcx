#
# An overly complicated makefile for compiling a swig module.
#

MODULENAME=ufc_benchmark

CXX=g++

ufcinclude=-I../../../ufc

# Python location and version
PYTHONROOT:=$(shell python -c 'import sys; print sys.prefix')
PYTHONVER:=$(shell python -c 'import sys; print sys.version[:3]')


all: _$(MODULENAME).so
	echo Done.

# link module
_$(MODULENAME).so: $(MODULENAME)_wrap.o $(MODULENAME).o
	$(CXX) -shared -o _$(MODULENAME).so $(MODULENAME)_wrap.o $(MODULENAME).o

# compile wrapper
$(MODULENAME)_wrap.o: $(MODULENAME)_wrap.cxx
	$(CXX) -I$(PYTHONROOT)/include/python$(PYTHONVER) -o $@ -c $<

# generate wrapper
$(MODULENAME)_wrap.cxx: $(MODULENAME).i $(MODULENAME).h
	swig -c++ -python $(ufcinclude) $(MODULENAME).i

# compile module code
$(MODULENAME).o: *.h *.cpp
	$(CXX) -c -o $(MODULENAME).o $(MODULENAME).cpp

clean:
	rm -f $(MODULENAME).o
	rm -f $(MODULENAME).py
	rm -f $(MODULENAME).pyc
	rm -f _$(MODULENAME).so
	rm -f $(MODULENAME)_wrap.cxx
	rm -f $(MODULENAME)_wrap.o
	rm -rf build
