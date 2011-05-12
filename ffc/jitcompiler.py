"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

# Copyright (C) 2007-2009 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC.  If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Johan Hake, 2008-2009
# Modified by Ilmar Wilbers, 2008
# Modified by Kristian B. Oelgaard, 2009
# Modified by Joachim Haga, 2011.
#
# First added:  2007-07-20
# Last changed: 2011-04-26

# Python modules
import os, sys
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

def jit(object, parameters=None, common_cell=None):
    """Just-in-time compile the given form or element

    Parameters:

      object     : The object to be compiled
      parameters : A set of parameters
    """

    # Check if we get an element or a form
    if isinstance(object, FiniteElementBase):
        return jit_element(object, parameters)
    else:
        return jit_form(object, parameters, common_cell)

def jit_form(form, parameters=None, common_cell=None):
    "Just-in-time compile the given form."

    print

    # Check that we get a Form
    if not isinstance(form, Form):
        form = as_form(form)

    # Check parameters
    parameters = _check_parameters(form, parameters)

    # Set log level
    set_level(parameters["log_level"])
    set_prefix(parameters["log_prefix"])

    # Compute form metadata and extract preprocessed form
    form_data = form.compute_form_data(common_cell=common_cell)
    preprocessed_form = form_data.preprocessed_form

    # Wrap input
    jit_object = JITObject(form, preprocessed_form, parameters, common_cell)

    # Set prefix for generated code
    prefix = "ffc_form_" + jit_object.signature()

    # Use Instant cache if possible
    cache_dir = parameters["cache_dir"]
    if cache_dir == "": cache_dir = None
    module = instant.import_module(jit_object, cache_dir=cache_dir)
    if module:

        print "Reusing form from cache!"
        debug("Reusing form from cache.")
        compiled_form = getattr(module, prefix + "_form_0")()
        return (compiled_form, module, form_data)

    try:

        # Take lock to serialise code generation and compilation.
        lock = instant.locking.get_lock(instant.get_default_cache_dir(),
                                        'ffc_' + jit_object.signature())

        # Retry Instant cache. The module may have been created while we waited
        # for the lock, even if it didn't exist before.
        module = instant.import_module(jit_object, cache_dir=cache_dir)
        if module:
            compiled_form = getattr(module, module.__name__ + "_form_0")()
            return (compiled_form, module, form_data)

        # Write a message
        log(INFO + 5,
            "Calling FFC just-in-time (JIT) compiler, this may take some time.")

        # Generate code
        compile_form(preprocessed_form,
                     prefix=prefix,
                     parameters=parameters,
                     common_cell=common_cell)

        # Build module using Instant (through UFC)
        debug("Compiling and linking Python extension module, this may take some time.")
        hfile = prefix + ".h"
        cppfile = prefix + ".cpp"
        module = ufc_utils.build_ufc_module(
            hfile,
            swig_binary=parameters["swig_binary"],
            swig_path=parameters["swig_path"],
            source_directory = os.curdir,
            signature = jit_object.signature(),
            sources = [cppfile] if parameters["split"] else [],
            cppargs = parameters["cpp_optimize_flags"].split() \
                      if parameters["cpp_optimize"] else ["-O0"],
            cache_dir = cache_dir)

        # Remove code
        os.unlink(hfile)
        if parameters["split"] :
            os.unlink(cppfile)

        # Extract compiled form
        compiled_form = getattr(module, prefix + "_form_0")()

        return compiled_form, module, form_data

    finally:
        instant.locking.release_lock(lock)

def jit_element(element, parameters=None):
    "Just-in-time compile the given element"

    # FIXME: We need a new solution for this.

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
            getattr(module, module.__name__ + ("_dofmap_%d" % i))())
