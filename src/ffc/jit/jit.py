"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-07-20 -- 2007-12-20"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

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

# Options for JIT-compiler, evaluate_basis and evaluate_basis_derivatives turned off
FFC_OPTIONS_JIT = FFC_OPTIONS.copy()
FFC_OPTIONS_JIT["no-evaluate_basis"] = True
FFC_OPTIONS_JIT["no-evaluate_basis_derivatives"] = True

# Compiler options, don't optimize by default (could be added to options)
CPP_ARGS = "-O0"

def jit(form, representation=FFC_REPRESENTATION, language=FFC_LANGUAGE, options=FFC_OPTIONS_JIT):
    "Just-in-time compile the given form or element"

    # Choose prefix
    global counter
    prefix = "ffc_form_%d" % counter
    counter += 1

    # Check that we don't get a list
    if isinstance(form, list):
        raise RuntimeError, "Just-in-time compiler requires a single form (not a list of forms"

    # Compile form
    debug("Calling FFC just-in-time (JIT) compiler, this may take some time...", -1)
    (form_data, form_representation) = compile(form, prefix, representation, language, options)
    debug("done", -1)

    # Filename of code to wrap
    filename = prefix + ".h"

    # FIXME: Move this to top when we have added dependence on Instant
    import instant
    instant.USE_CACHE = 1

    # Get include directory for ufc.h (might be better way to do this?)
    (path, dummy, dummy, dummy) = instant.header_and_libs_from_pkgconfig("ufc-1")

    if len(path) == 0:
        path = [("/").join(sysconfig.get_python_inc().split("/")[:-2]) + "/include"]
    ufc_include = '%%include "%s/ufc.h"' % path[0]

    # Wrap code into a Python module using Instant
    debug("Creating Python extension (compiling and linking), this may take some time...", -1)
    module_name = prefix + "_module"
    instant.create_extension(wrap_headers=[filename], module=module_name, additional_declarations=ufc_include, include_dirs=path, cppargs=CPP_ARGS)
    debug("done", -1)

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
