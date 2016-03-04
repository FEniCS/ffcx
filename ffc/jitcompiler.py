"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

# Copyright (C) 2007-2015 Anders Logg
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Johan Hake, 2008-2009
# Modified by Ilmar Wilbers, 2008
# Modified by Kristian B. Oelgaard, 2009
# Modified by Joachim Haga, 2011.
# Modified by Martin Alnaes, 2013-2015

# Python modules
import os, sys
import instant

# UFL modules
from ufl import TestFunction, ds, dx
from ufl.classes import Form, FiniteElementBase
from ufl.algorithms import extract_elements, extract_sub_elements, compute_form_data

# FFC modules
from ffc.log import log
from ffc.log import info
from ffc.log import warning
from ffc.log import debug
from ffc.log import error
from ffc.log import set_level
from ffc.log import set_prefix
from ffc.log import INFO
from ffc.parameters import default_parameters
from ffc.mixedelement import MixedElement
from ffc.compiler import compile_form, compile_element
from ffc.jitobject import JITObject
from ffc.quadratureelement import default_quadrature_degree
from ffc.backends.ufc import build_ufc_module

# Special Options for JIT-compilation
FFC_PARAMETERS_JIT = default_parameters()
FFC_PARAMETERS_JIT["no-evaluate_basis_derivatives"] = True

# Set debug level for Instant
instant.set_log_level("warning")


def check_swig_version(compiled_module):
    # Check swig version of compiled module
    import ufc
    if compiled_module and compiled_module.swigversion != ufc.__swigversion__:
        error("Incompatible swig versions detected. UFC swig "\
              "version is not the same as extension module swig "\
              "version: '%s' != '%s' " % \
              (ufc.__swigversion__, compiled_module.swigversion))


def jit(ufl_object, parameters=None):
    """Just-in-time compile the given form or element

    Parameters:

      ufl_object : The UFL object to be compiled
      parameters : A set of parameters
    """
    form = ufl_object

    # Check that we get a Form
    if isinstance(form, Form):
        kind = "form"
    elif isinstance(form, FiniteElementBase):
        kind = "element"
    else:
        error("Expecting a UFL form or element: %s" % repr(form))

    # Check parameters
    parameters = _check_parameters(form, parameters)

    # Set log level
    set_level(parameters["log_level"])
    set_prefix(parameters["log_prefix"])

    # Wrap input
    jit_object = JITObject(form, parameters)

    # Set prefix for generated code
    module_name = "ffc_%s_%s" % (kind, jit_object.signature())

    # Use Instant cache if possible
    cache_dir = parameters["cache_dir"] or None
    module = instant.import_module(module_name, cache_dir=cache_dir)
    if module:
        debug("Reusing form from cache.")
    else:
        # Take lock to serialise file removal.
        # Need to add "_0" to lock as instant.import_module acquire
        # lock with name: module_name
        with instant.file_lock(instant.get_default_cache_dir(),
                               module_name + "_0") as lock:

            # Retry Instant cache. The module may have been created while we waited
            # for the lock, even if it didn't exist before.
            module = instant.import_module(module_name, cache_dir=cache_dir)
            if module:
                debug("Reusing form from cache.")
            else:
                # Write a message
                log(INFO + 5,
                    "Calling FFC just-in-time (JIT) compiler, this may take some time.")

                # Generate code
                if kind == "form":
                    compile_form(form,
                                 prefix=module_name,
                                 parameters=parameters,
                                 jit=True)
                elif kind == "element":
                    compile_element(form,
                                    prefix=module_name,
                                    parameters=parameters,
                                    jit=True)

                # Build module using Instant (through UFC)
                debug("Compiling and linking Python extension module, this may take some time.")
                hfile   = module_name + ".h"
                cppfile = module_name + ".cpp"

                if parameters["cpp_optimize"]:
                    cppargs = parameters["cpp_optimize_flags"].split()
                else:
                    cppargs = ["-O0"]

                module = build_ufc_module(
                    hfile,
                    source_directory = os.curdir,
                    signature = module_name,
                    sources = [cppfile],
                    cppargs = cppargs,
                    cache_dir = cache_dir)

                # Remove code
                if os.path.isfile(hfile):
                    os.unlink(hfile)
                if os.path.isfile(cppfile):
                    os.unlink(cppfile)

    # Construct instance of compiled form
    check_swig_version(module)
    prefix = module_name
    if kind == "form":
        compiled_form = _instantiate_form(module, prefix)
        return compiled_form, module, prefix
    elif kind == "element":
        return _instantiate_element_and_dofmap(module, prefix)


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

from ffc.cpp import make_classname
def _instantiate_form(module, prefix):
    "Extract the form from module with only one form."
    form_id = 0
    classname = make_classname(prefix, "form", form_id)
    return getattr(module, "create_" + classname)()

def _instantiate_element_and_dofmap(module, prefix):
    """Extract element and dofmap from module."""
    fe_classname = make_classname(prefix, "finite_element", "main")
    dm_classname = make_classname(prefix, "dofmap", "main")
    fe = getattr(module, "create_" + fe_classname)()
    dm = getattr(module, "create_" + dm_classname)()
    return (fe, dm)
