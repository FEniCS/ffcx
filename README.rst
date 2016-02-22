-----------------------------
FFC: The FEniCS Form Compiler
-----------------------------

FFC is a compiler for finite element variational forms. From a
high-level description of the form, it generates efficient low-level
C++ code that can be used to assemble the corresponding discrete
operator (tensor). In particular, a bilinear form may be assembled
into a matrix and a linear form may be assembled into a vector.

FFC may be used either from the command line (by invoking the ``ffc``
command) or as a Python module (``import ffc``).

FFC is part of the FEniCS project (http://www.fenicsproject.org) and
functions as a just-in-time (JIT) compiler for DOLFIN.

For further introduction to FFC, open the FFC user manual available in
the subdirectory ``doc/manual/`` of this source tree, or try out the
demos available in the subdirectory ``src/demo/`` of this source tree.


License
-------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this program. If not, see
<http://www.gnu.org/licenses/>.


Dependencies
------------

#. Python, version 2.7 or later

#. The latest version of FIAT, Instant and UFL

   You need to have FIAT, Instant and UFL installed. They are
   available from the web page: https://bitbucket.org/fenics-project/.

#. The Python NumPy module

#. The Python Six module

#. SWIG, version 3.0.3 or higher


Notes
-----

From February 2014, the code generation interface UFC is distributed
as part of FFC, and the UFC repository has been merged into the FFC
repository. From this point onwards, UFC version numbers are reset to
the same version numbers as for FFC.
