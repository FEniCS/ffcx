"""Code generation format strings for UFC (Unified Form-assembly Code) v. 2.3.0.

Three format strings are defined for each of the following UFC classes:

    function
    finite_element
    dofmap
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

For more information about UFC and the FEniCS Project, visit

    http://www.fenicsproject.org
    https://bitbucket.org/fenics-project/ufc

"""

# -*- coding: utf-8 -*-
__author__  = "Martin Sandve Alnaes, Anders Logg, Kent-Andre Mardal, Ola Skavhaug, and Hans Petter Langtangen"
__date__    = "2014-01-07"
__version__ = "2.3.0"
__license__ = "This code is released into the public domain"

UFC_VERSION_MAJOR = 2
UFC_VERSION_MINOR = 3
UFC_VERSION_MAINTENANCE = 0

UFC_VERSION = __version__

from .function import *
from .finite_element import *
from .dofmap import *
from .integrals import *
from .form import *
from .build import build_ufc_module

# Check if we should import development/experimental version
#import os
#if "UFC_DEV" in os.environ:
#    from finite_element_1_6 import *
#    from dofmap_1_6 import *

templates = {"function_header":                        function_header,
             "function_implementation":                function_implementation,
             "function_combined":                      function_combined,
             "finite_element_header":                  finite_element_header,
             "finite_element_implementation":          finite_element_implementation,
             "finite_element_combined":                finite_element_combined,
             "dofmap_header":                          dofmap_header,
             "dofmap_implementation":                  dofmap_implementation,
             "dofmap_combined":                        dofmap_combined,
             "cell_integral_header":                   cell_integral_header,
             "cell_integral_implementation":           cell_integral_implementation,
             "cell_integral_combined":                 cell_integral_combined,
             "exterior_facet_integral_header":         exterior_facet_integral_header,
             "exterior_facet_integral_implementation": exterior_facet_integral_implementation,
             "exterior_facet_integral_combined":       exterior_facet_integral_combined,
             "interior_facet_integral_header":         interior_facet_integral_header,
             "interior_facet_integral_implementation": interior_facet_integral_implementation,
             "interior_facet_integral_combined":       interior_facet_integral_combined,
             "point_integral_header":                  point_integral_header,
             "point_integral_implementation":          point_integral_implementation,
             "point_integral_combined":                point_integral_combined,
             "form_header":                            form_header,
             "form_implementation":                    form_implementation,
             "form_combined":                          form_combined}
