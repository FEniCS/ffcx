# -*- coding: utf-8 -*-
"""Code generation format strings for UFC (Unified Form-assembly Code) version 2016.2.0.dev0

Five format strings are defined for each of the following UFC classes:

    cell_integral
    exterior_facet_integral
    interior_facet_integral
    custom_integral
    cutcell_integral
    interface_integral
    overlap_integral
    function

    finite_element
    dofmap
    coordinate_mapping
    form

The strings are named:

    '<classname>_header'
    '<classname>_implementation'
    '<classname>_combined'
    '<classname>_jit_header'
    '<classname>_jit_implementation'

The header and implementation contain the definition and declaration
of the class respectively, and are meant to be placed in .h and .cpp files,
while the combined version is for an implementation within a single .h header.
The _jit_ versions are used in the jit compiler and contains some additional
factory functions exported as extern "C" to allow construction of compiled
objects through ctypes without dealing with C++ ABI and name mangling issues.

Each string has at least the following format variables: 'classname',
'members', 'constructor', 'destructor', plus one for each interface
function with name equal to the function name.

For more information about UFC and the FEniCS Project, visit

    http://www.fenicsproject.org
    https://bitbucket.org/fenics-project/ffc

"""

__author__ = "Martin Sandve Alnæs, Anders Logg, Kent-Andre Mardal, Ola Skavhaug, and Hans Petter Langtangen"
__date__ = "2016-09-27"
__version__ = "2016.2.0.dev0"
__license__ = "This code is released into the public domain"

from os.path import dirname, abspath

from .function import *
from .finite_element import *
from .dofmap import *
from .coordinate_mapping import *
from .integrals import *
from .form import *


def get_include_path():
    "Return location of UFC header files"
    return dirname(abspath(__file__))


# Platform specific snippets for controlling visilibity of exported symbols in generated shared libraries
visibility_snippet = """
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif
"""


# Generic factory function signature
factory_decl = """
extern "C" %(basename)s * create_%(publicname)s();
"""


# Generic factory function implementation. Set basename to the base class,
# and note that publicname and privatename does not need to match, allowing
# multiple factory functions to return the same object.
factory_impl = """
extern "C" DLL_EXPORT %(basename)s * create_%(publicname)s()
{
 return new %(privatename)s();
}
"""


def all_ufc_classnames():
    "Build list of all classnames."
    integral_names = ["cell", "exterior_facet", "interior_facet", "vertex", "custom", "cutcell", "interface", "overlap"]
    integral_classnames = [integral_name + "_integral" for integral_name in integral_names]
    jitable_classnames = ["finite_element", "dofmap", "coordinate_mapping", "form"]
    classnames = ["function"] + jitable_classnames + integral_classnames
    return classnames


def _build_templates():
    "Build collection of all templates to store in the templates dict."
    templates = {}
    classnames = all_ufc_classnames()

    for classname in classnames:
        # Expect all classes to have header, implementation, and combined versions
        header = globals()[classname + "_header"]
        implementation = globals()[classname + "_implementation"]
        combined = globals()[classname + "_combined"]

        # Construct jit header with just factory function signature
        jit_header = factory_decl % {
            "basename": "ufc::" + classname,
            "publicname": "%(classname)s",
            "privatename": "%(classname)s",
            }

        # Construct jit implementation template with class declaration,
        # factory function implementation, and class definition
        _fac_impl = factory_impl % {
            "basename": "ufc::" + classname,
            "publicname": "%(classname)s",
            "privatename": "%(classname)s",
            }
        jit_implementation = header + _fac_impl + implementation

        # Store all in templates dict
        templates[classname + "_header"] = header
        templates[classname + "_implementation"] = implementation
        templates[classname + "_combined"] = combined
        templates[classname + "_jit_header"] = jit_header
        templates[classname + "_jit_implementation"] = jit_implementation

    return templates

templates = _build_templates()
