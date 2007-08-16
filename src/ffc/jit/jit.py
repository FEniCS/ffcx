"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-07-20 -- 2007-08-16"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python module
from os import system
from commands import getoutput
from distutils import sysconfig

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *

# FFC compiler modules
from ffc.compiler.compiler import compile

# Global counter for numbering forms
counter = 0

def jit(form, representation=FFC_REPRESENTATION, language=FFC_LANGUAGE, options=FFC_OPTIONS):
    "Just-in-time compile the given form or element"

    # Choose prefix
    global counter
    prefix = "ffc_form_%d" % counter
    counter += 1

    # Check that we don't get a list
    if isinstance(form, list):
        raise RuntimeError, "Just-in-time compiler requires a single form (not a list of forms"

    # Compile form
    (form_data, form_representation) = compile(form, prefix, representation, language, options)

    # Filename of code to wrap
    filename = "../" + prefix + ".h"

    # FIXME: Move this to top when we have added dependence on Instant
    from instant import create_extension, header_and_libs_from_pkgconfig

    # Get include directory for ufc.h (might be better way to do this?)
    (path, dummy, dummy, dummy) = header_and_libs_from_pkgconfig("ufc-1")
    if len(path) == 0:
        path = [("/").join(sysconfig.get_python_inc().split("/")[:-2]) + "/include"]
    ufc_include = '%%include "%s/ufc.h"' % path[0]

    # Wrap code into a Python module using Instant
    module_name = prefix + "_module"
    create_extension(wrap_headers=[filename], module=module_name, additional_declarations=ufc_include)

    # Get name of form
    rank = form_data[0].rank
    if rank == 0:
        form_name = prefix + "Functional"
    elif rank == 1:
        form_name = prefix + "LinearForm"
    elif rank == 2:
        form_name = prefix + "BilinearForm"
    else:
        form_name = prefix

    # Return the form, module and form data
    exec("import %s as compiled_module" % module_name)
    exec("compiled_form = compiled_module.%s()" % form_name)
    return (compiled_form, compiled_module, form_data)
