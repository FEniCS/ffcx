UFLACS - UFL Analyser and Compiler System
=========================================

Description
-----------

Uflacs, the UFL Analyser and Compiler System, is a collection of
algorithms for processing symbolic UFL forms and expressions.
The main feature is efficient translation of tensor intensive
symbolic expressions into a low level expression representation and C++ code.


Licencing
---------

Uflacs is Copyright (2011-2014) by Martin Sandve Alnæs.

This version of uflacs is released under the LGPL v3 licence.


Installing
----------

Either install to default python location as root::

    sudo python setup.py install

Or install to your own python path directory::

    python setup.py install --prefix=/path/to/my/own/site-packages


Testing
-------

To run the Python tests you need to install the py.test framework.
Then run all tests simply by executing

    cd test && py.test

To run unittests of generated C++ code you also need the Google C++ Testing Framework:

    https://code.google.com/p/googletest/


Contact
-------

Bitbucket site:

    http://www.bitbucket.org/fenics-project/uflacs/

FEniCS Project site:

    http://www.fenicsproject.org/

Author:

    Martin Sandve Alnæs (martinal@simula.no)

