"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-07-20 -- 2008-09-16"
__copyright__ = "Copyright (C) 2007-2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Johan Hake, 2008
# Mofified by Ilmar Wilbers, 2008

# Python modules
import instant
import distutils
import os

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *

# FFC fem modules
from ffc.fem.finiteelement import FiniteElement
from ffc.fem.mixedelement import MixedElement

# FFC compiler modules
from ffc.compiler.compiler import compile
from ffc.compiler.language.algebra import Form, TestFunction
from ffc.compiler.language.builtins import dx

# FFC jit modules
from jitobject import JITObject

# Special Options for JIT-compilation
FFC_OPTIONS_JIT = FFC_OPTIONS.copy()
FFC_OPTIONS_JIT["no-evaluate_basis_derivatives"] = True

# Set debug level for Instant
instant.set_logging_level("warning")

def jit(object, options=None):
    """Just-in-time compile the given form or element
    
    Parameters:
    
      object  : The object to be compiled
      options : An option dictionary
    """

    # Check if we get an element or a form
    if isinstance(object, FiniteElement) or isinstance(object, MixedElement):
        return jit_element(object)
    else:
        return jit_form(object)

def jit_form(form, options=None):
    "Just-in-time compile the given form"

    # Make sure that we get a form
    if not isinstance(form, Form):
        form = Form(form)

    # Check options
    options = check_options(form, options)

    # Wrap input
    jit_object = JITObject(form, options)

    # Check cache
    module = instant.import_module(jit_object, cache_dir=options["cache_dir"])
    if module: return extract_form(form, module)

    # Generate code
    debug("Calling FFC just-in-time (JIT) compiler, this may take some time...", -1)
    signature = jit_object.signature()
    compile(form, signature, options)
    debug("done", -1)

    # Wrap code into a Python module using Instant
    debug("Creating Python extension (compiling and linking), this may take some time...", -1)
    filename = signature + ".h"
    (cppargs, path, ufc_include) = extract_instant_flags(options)
    module = instant.build_module(wrap_headers=[filename],
                                  additional_declarations=ufc_include,
                                  include_dirs=path,
                                  cppargs=cppargs,
                                  signature=signature,
                                  cache_dir=options["cache_dir"])
    os.unlink(filename)
    debug("done", -1)

    return extract_form(form, module)

def jit_element(element, options=None):
    "Just-in-time compile the given element"
    
    # Check that we get an element
    if not (isinstance(element, FiniteElement) or isinstance(element, MixedElement)):
        raise RuntimeError, "Expecting a finite element."

    # Create simplest possible dummy form
    if element.value_dimension(0) > 1:
        form = TestFunction(element)[0]*dx
    else:
        form = TestFunction(element)*dx

    # Compile form
    (compiled_form, module, form_data) = jit_form(form, options)

    # Extract element and dofmap
    return extract_element_and_dofmap(module)

def check_options(form, options):
    "Check options and add any missing options"

    # Form can not be a list
    if isinstance(form, list):
        raise RuntimeError, "JIT compiler requires a single form (not a list of forms)."

    # Copy options
    if options is None:
        options = {}
    else:
        options = options.copy()

    # Check for invalid options
    for key in options:
        if not key in FFC_OPTIONS_JIT:
            warning('Unknown option "%s" for JIT compiler, ignoring.' % key)

    # Add defaults for missing options
    for key in FFC_OPTIONS_JIT:
        if not key in options:
            options[key] = FFC_OPTIONS_JIT[key]

    # Don't postfix form names
    options["form_postfix"] = False

    return options

def extract_form(form, module):
    "Extract form from module"
    return (getattr(module, module.__name__)(), module, form.form_data)

def extract_element_and_dofmap(module):
    "Extract element and dofmap from module"
    name = module.__name__
    element = getattr(module, name + "_finite_element_0")()
    dofmap  = getattr(module, name + "_dof_map_0")()
    return (element, dofmap)

def extract_instant_flags(options):
    "Extract flags for Instant"

    # Get C++ compiler options
    if options["cpp optimize"]: cppargs = "-O2"
    else: cppargs = "-O0"

    # Get include directory for ufc.h (might be better way to do this?)
    (path, dummy, dummy, dummy) = instant.header_and_libs_from_pkgconfig("ufc-1")
    if len(path) == 0: path = [("/").join(distutils.sysconfig.get_python_inc().split("/")[:-2]) + "/include"]
    ufc_include = '%%include "%s/ufc.h"' % path[0]

    return (cppargs, path, ufc_include)
