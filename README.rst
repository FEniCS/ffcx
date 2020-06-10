===============================
FFC-X: The FEniCS Form Compiler
===============================

.. image:: https://github.com/FEniCS/ffcx/workflows/FFCX%20CI/badge.svg
   :target: https://github.com/FEniCS/ffcx/workflows/FFCX%20CI
.. image:: https://readthedocs.org/projects/fenics-ffcx/badge/?version=latest
   :target: http://fenics-ffcx.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://coveralls.io/repos/github/FEniCS/ffcx/badge.svg?branch=master
   :target: https://coveralls.io/github/FEniCS/ffcx?branch=master

FFC-X is a new version of the FEniCS Form Compiler. It is being actively
developed and is compatible with DOLFIN-X.

FFC-X is a compiler for finite element variational forms. From a
high-level description of the form in the Unified Form Language (UFL),
it generates efficient low-level C code that can be used to assemble the
corresponding discrete operator (tensor). In particular, a bilinear form
may be assembled into a matrix and a linear form may be assembled into a
vector.  FFCX may be used either from the command line (by invoking the
``ffcx`` command) or as a Python module (``import ffcx``).

FFC-X is part of the FEniCS Project. For more information, visit
http://www.fenicsproject.org


Documentation
=============

Documentation can be viewed at http://fenics-ffcx.readthedocs.org/.


License
=======

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.
