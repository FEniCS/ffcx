"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-07-20 -- 2007-07-20"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *

# FFC compiler modules
from ffc.compiler.compiler import compile

# Global counter for numbering forms
counter = 0

def jit(form, representation=FFC_REPRESENTATION, output_language=FFC_LANGUAGE, options=FFC_OPTIONS, return_module=False):
    "Just-in-time compile the given form or element"

    # Choose prefix
    global counter
    prefix = "ffc_form_%d" % counter
    counter += 1

    # Check that we don't get a list
    if isinstance(form, list):
        raise RuntimeError, "Just-in-time compiler requires a single form (not a list of forms"

    # Compile form
    form_data = compile(form, prefix, representation, output_language, options)

    # Get code as a string
    file = open(prefix + ".h")
    code = file.read()
    file.close()

    # Wrap code into a Python module using Instant
    from instant import create_extension
    module_name = prefix + "_module"
    create_extension(code=code, module=module_name)

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

    # Return the module and the form
    exec("import %s as compiled_module" % module_name)
    exec("compiled_form = compiled_module.%s()" % form_name)
    if return_module:
        return (compiled_form, compiled_module)
    else:
        return compiled_form
