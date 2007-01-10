"""Code generation format strings for UFC (Unified Form-assembly Code) v. 1.0.

Three format strings are defined for each of the following UFC classes:

    function
    dof_map
    finite_element
    cell_integral
    exterior_facet_integral
    interior_facet_integral

The strings are named '<classname>_header', '<classname>_implementation',
and '<classname>_combined'. The header and implementation contain the
definition and declaration respectively, and are meant to be placed in
.h and .cpp files, while the combined version is for an implementation
within a single .h header.

Each string has the following format variables: 'classname',
'members', 'constructor', 'destructor', plus one for each interface
function with name equal to the function name.

For more information about UFC and the FEniCS project, visit

    http://www.fenics.org/ufc/

"""

# -*- coding: utf-8 -*-
__author__  = "Martin Sandve Alnaes, Hans Petter Langtangen, Anders Logg, Kent-Andre Mardal and Ola Skavhaug"
__date__    = "2007-01-10"
__version__ = "1.0-rc1"
__license__ = "This code is released into the public domain"

UFC_VERSION = __version__

from dof_map import *
from finite_element import *
from integrals import *
from form import *
