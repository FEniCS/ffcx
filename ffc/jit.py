"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-07-20"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Johan Hake, 2008-2009
# Modified by Ilmar Wilbers, 2008
# Modified by Kristian B. Oelgaard, 2009
# Last changed: 2010-02-01

# Python modules
import os
import instant
import ufc_utils

# UFL modules
from ufl.classes import Form, FiniteElementBase, TestFunction
from ufl.objects import dx
from ufl.algorithms import as_form, preprocess, FormData

# FFC modules
from log import log
from log import info
from log import warning
from log import debug
from log import error
from log import set_level
from log import set_prefix
from log import INFO
from parameters import default_parameters
from mixedelement import MixedElement
from compiler import compile_form
from jitobject import JITObject

# Special Options for JIT-compilation
FFC_PARAMETERS_JIT = default_parameters()
FFC_PARAMETERS_JIT["no-evaluate_basis_derivatives"] = True

# Set debug level for Instant
instant.set_logging_level("warning")

# Form data cache
_form_data_cache = {}

def jit(object, parameters=None):
    """Just-in-time compile the given form or element

    Parameters:

      object  : The object to be compiled
      parameters : An option dictionary
    """

    # Check if we get an element or a form
    if isinstance(object, FiniteElementBase):
        return jit_element(object, parameters)
    else:
        return jit_form(object, parameters)

def jit_form(form, parameters=None):
    "Just-in-time compile the given form"

    # Check that we get a Form
    if not isinstance(form, Form):
        form = as_form(form)

    # Check parameters
    parameters = _check_parameters(form, parameters)

    # Set log level
    set_level(parameters["log_level"])
    set_prefix(parameters["log_prefix"])

    # Preprocess form
    preprocessed_form = preprocess(form)

    # Wrap input
    jit_object = JITObject(form, preprocessed_form, parameters)

    # Get cache dir
    cache_dir = parameters["cache_dir"]
    if cache_dir == "": cache_dir = None

    # Check cache
    module = instant.import_module(jit_object, cache_dir=cache_dir)
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
    if parameters["representation"] == "auto":
        representation = ""
    else:
        representation = " (%s representation)" % parameters["representation"]
    log(INFO + 5, "Calling FFC just-in-time (JIT) compiler%s, this may take some time." % representation)

    # Generate code
    signature = jit_object.signature()
    analysis = compile_form(preprocessed_form, prefix=signature, parameters=parameters)
    # FIXME: Make this more sane, less cryptic
    form_data = analysis[0][0][1]

    # Create python extension module using Instant (through UFC)
    debug("Creating Python extension (compiling and linking), this may take some time...")
    module = ufc_utils.build_ufc_module(
        signature + ".h",
        source_directory = os.curdir,
        signature = signature,
        sources = [signature + ".cpp"] if parameters["split"] else [],
        cppargs  = ["-O2"] if parameters["cpp_optimize"] else ["-O0"] ,
        cache_dir = cache_dir)

    # Remove code
    os.unlink(signature + ".h")
    if parameters["split"] :
        os.unlink(signature + ".cpp")

    # Extract compiled form
    compiled_form = getattr(module, module.__name__ + "_form_0")()

    # Store form data in cache
    _form_data_cache[form] = form_data

    return compiled_form, module, form_data

def jit_element(element, parameters=None):
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
    (compiled_form, module, form_data) = jit_form(form, parameters)

    return _extract_element_and_dofmap(module, form_data)

def _check_parameters(form, parameters):
    "Check parameters and add any missing parameters"

    # Form can not be a list
    if isinstance(form, list):
        error("JIT compiler requires a single form (not a list of forms).")

    # Copy parameters
    if parameters is None:
        parameters = {}
    else:
        parameters = parameters.copy()

    # Check for invalid parameters
    for key in parameters:
        if not key in FFC_PARAMETERS_JIT:
            warning('Unknown option "%s" for JIT compiler, ignoring.' % key)

    # Add defaults for missing parameters
    for key in FFC_PARAMETERS_JIT:
        if not key in parameters:
            parameters[key] = FFC_PARAMETERS_JIT[key]

    # Don't postfix form names
    parameters["form_postfix"] = False

    return parameters

def _extract_element_and_dofmap(module, form_data):
    """
    Extract element and dofmap from module. Code will be generated for
    all unique elements (including sub elements) and to get the top
    level element we need to extract the last element.
    """
    i = len(form_data.unique_sub_elements) - 1
    return (getattr(module, module.__name__ + ("_finite_element_%d" % i))(),
            getattr(module, module.__name__ + ("_dof_map_%d" % i))())
