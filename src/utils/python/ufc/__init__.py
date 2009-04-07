"""Code generation format strings for UFC (Unified Form-assembly Code) v. 1.1.2.

Three format strings are defined for each of the following UFC classes:

    function
    finite_element
    dof_map
    cell_integral
    exterior_facet_integral
    interior_facet_integral
    form

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
__author__  = "Martin Sandve Alnaes, Anders Logg, Kent-Andre Mardal, Ola Skavhaug, and Hans Petter Langtangen"
__date__    = "2009-04-07"
__version__ = "1.1"
__license__ = "This code is released into the public domain"

UFC_VERSION_MAJOR = 1
UFC_VERSION_MINOR = 1
UFC_VERSION_MAINTENANCE = 1

UFC_VERSION = __version__

from function import *
from finite_element import *
from dof_map import *
from integrals import *
from form import *
from build import build_ufc_module
