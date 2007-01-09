"""Defines a set of strings to help in the generation of UFC compliant code.

Three format strings are defined for each of the classes:
function, dof_map, finite_element, cell_integral, exterior_facet_integral, and interior_facet_integral

The strings are named 'classname_header', 'classname_implementation', and 'classname_combined'.
The header and implementation contain the definition and declaration respectively, and are meant to
be placed in .h and .cpp files, while the combined version is for an implementation within the header.

Each string has the following format variables:
'classname', 'members', 'constructor', 'destructor',
plus one for each interface function with name equal to the function name.

In addition, a utility function generate_code(format_string, format_dict, num_indent_spaces=6) allows
the user to specify a format dictionary without the constructor, destructor or extra members.
"""
# -*- coding: utf-8 -*-
__author__  = "Martin Sandve Alnæs"
__date__    = "9th January 2007"
__version__ = "1.0-rc1"

UFC_VERSION = __version__


shapes = ["interval", "triangle", "quadrilateral", "tetrahedron", "hexahedron"]

from dof_map import *
from finite_element import *
from integrals import *
from form import *

import utility

