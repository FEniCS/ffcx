"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-07-20"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Johan Hake, 2008-2009
# Modified by Ilmar Wilbers, 2008
# Modified by Kristian B. Oelgaard, 2009
# Last changed: 2009-12-16

# Python modules.
import os
import instant
import ufc_utils

# UFL modules.
from ufl.classes import Form
from ufl.classes import FiniteElementBase
from ufl.classes import TestFunction
from ufl.objects import dx
from ufl.algorithms import as_form
from ufl.algorithms import preprocess
from ufl.algorithms import FormData

# FFC modules.
from log import log
from log import info
from log import warning
from log import debug
from log import error
from log import set_level
from log import set_prefix
from log import INFO
from constants import FFC_OPTIONS
from mixedelement import MixedElement
from compiler import compile_form
from jitobject import JITObject

# Special Options for JIT-compilation
FFC_OPTIONS_JIT = FFC_OPTIONS.copy()
FFC_OPTIONS_JIT["no-evaluate_basis_derivatives"] = True

# Set debug level for Instant
instant.set_logging_level("warning")

# Form data cache
_form_data_cache = {}

def jit(object, options=None):
    """Just-in-time compile the given form or element

    Parameters:

      object  : The object to be compiled
      options : An option dictionary
    """

    # Check if we get an element or a form
    if isinstance(object, FiniteElementBase):
        return jit_element(object, options)
    else:
        return jit_form(object, options)

def jit_form(form, options=None):
    "Just-in-time compile the given form"

    # Check that we get a Form
    if not isinstance(form, Form):
        form = as_form(form)

    # Check options
    options = check_options(form, options)

    # Set log level
    set_level(options["log_level"])
    set_prefix(options["log_prefix"])

    # Preprocess form
    preprocessed_form = preprocess(form)

    # Wrap input
    jit_object = JITObject(form, preprocessed_form, options)

    # Check cache
    module = instant.import_module(jit_object, cache_dir=options["cache_dir"])
    if module:

        # Get compiled form from Instant cache
        compiled_form = getattr(module, module.__name__ + "_form_0")()

        # Get form data from in-memory cache or create it
        if form in _form_data_cache:
            form_data = _form_data_cache[form]
        else:
            form_data = FormData(preprocessed_form)
            _form_data_cache[form] = form_data
        return (compiled_form, module, form_data)

    # Write a message
    if options["representation"] == "auto":
        representation = ""
    else:
        representation = " (%s representation)" % options["representation"]
    log(INFO + 5, "Calling FFC just-in-time (JIT) compiler%s, this may take some time." % representation)

    # Generate code
    signature = jit_object.signature()
    preprocessed_form, form_data = compile_form(preprocessed_form, signature, options)[0]

    # Create python extension module using Instant (through UFC)
    debug("Creating Python extension (compiling and linking), this may take some time...")
    module = ufc_utils.build_ufc_module(
        signature + ".h",
        source_directory = os.curdir,
        signature = signature,
        sources = [signature + ".cpp"] if options["split"] else [],
        cppargs  = ["-O2"] if options["cpp optimize"] else ["-O0"] ,
        cache_dir = options["cache_dir"])

    # Remove code
    os.unlink(signature + ".h")
    if options["split"] :
        os.unlink(signature + ".cpp")

    # Extract compiled form
    compiled_form = getattr(module, module.__name__ + "_form_0")()

    # Store form data in cache
    _form_data_cache[form] = form_data

    return compiled_form, module, form_data

def jit_element(element, options=None):
    "Just-in-time compile the given element"

    # Check that we get an element
    if not isinstance(element, FiniteElementBase):
        error("Expecting a finite element.")

    # Create simplest possible dummy form
    v = TestFunction(element)
    for i in range(len(element.value_shape())):
        v = v[i]
    form = v*dx

    # Compile form
    (compiled_form, module, form_data) = jit_form(form, options)

    return extract_element_and_dofmap(module)

def check_options(form, options):
    "Check options and add any missing options"

    # Form can not be a list
    if isinstance(form, list):
        error("JIT compiler requires a single form (not a list of forms).")

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

def extract_element_and_dofmap(module):
    "Extract element and dofmap from module"
    return (getattr(module, module.__name__ + "_0_finite_element_0")(),
            getattr(module, module.__name__ + "_0_dof_map_0")())
