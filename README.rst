=============================
FFC: The FEniCS Form Compiler
=============================

FFC is a compiler for finite element variational forms. From a
high-level description of the form, it generates efficient low-level
C++ code that can be used to assemble the corresponding discrete
operator (tensor). In particular, a bilinear form may be assembled
into a matrix and a linear form may be assembled into a vector.  FFC
may be used either from the command line (by invoking the ``ffc``
command) or as a Python module (``import ffc``).

FFC is part of the FEniCS Project.

For more information, visit http://www.fenicsproject.org


Documentation
=============

Documentation can be viewed at http://fenics-ffc.readthedocs.org/.

.. image:: https://readthedocs.org/projects/fenics-ffc/badge/?version=latest
   :target: http://fenics.readthedocs.io/projects/ffc/en/latest/?badge=latest
   :alt: Documentation Status


Automated Testing
=================

We use Bitbucket Pipelines and Atlassian Bamboo to perform automated
testing.

.. image:: https://bitbucket-badges.useast.atlassian.io/badge/fenics-project/ffc.svg
   :target: https://bitbucket.org/fenics-project/ffc/addon/pipelines/home
   :alt: Pipelines Build Status

.. image:: http://fenics-bamboo.simula.no:8085/plugins/servlet/wittified/build-status/FFC-FD
   :target: http://fenics-bamboo.simula.no:8085/browse/FFC-FD/latest
   :alt: Bamboo Build Status


Code Coverage
=============

Code coverage reports can be viewed at
https://coveralls.io/bitbucket/fenics-project/ffc.

.. image:: https://coveralls.io/repos/bitbucket/fenics-project/ffc/badge.svg?branch=master
   :target: https://coveralls.io/bitbucket/fenics-project/ffc?branch=master
   :alt: Coverage Status


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
