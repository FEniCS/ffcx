"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

# Copyright (C) 2007-2013 Anders Logg
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
# Modified by Martin Alnaes, 2013
#
# First added:  2007-07-20
# Last changed: 2013-01-25

# Python modules
import os, sys
import instant
import ufc_utils
import ufc

# UFL modules
from ufl.classes import Form, FiniteElementBase, TestFunction
from ufl.domains import as_domain
from ufl.objects import dx
from ufl.algorithms import as_form, extract_common_cell, extract_elements, extract_sub_elements
from ufl.common import istr, tstr

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
from ffc.quadratureelement import default_quadrature_degree

# Special Options for JIT-compilation
FFC_PARAMETERS_JIT = default_parameters()
FFC_PARAMETERS_JIT["no-evaluate_basis_derivatives"] = True

# Set debug level for Instant
instant.set_log_level("warning")

def jit(ufl_object, parameters=None, common_cell=None):
    """Just-in-time compile the given form or element

    Parameters:

      ufl_object : The UFL object to be compiled
      parameters : A set of parameters
    """

    # Check if we get an element or a form
    if isinstance(ufl_object, FiniteElementBase):
        return jit_element(ufl_object, parameters)
    else:
        return jit_form(ufl_object, parameters, common_cell)

def _auto_select_degree(elements):
    """
    Automatically select degree for all elements of the form in cases
    where this has not been specified by the user. This feature is
    used by DOLFIN to allow the specification of Expressions with
    undefined degrees.
    """

    # Extract common degree
    common_degree = max([e.degree() for e in elements] or [None])
    if common_degree is None:
        common_degree = default_quadrature_degree

    # Degree must be at least 1 (to work with Lagrange elements)
    common_degree = max(1, common_degree)

    return common_degree

def _compute_element_mapping(form, common_cell):
    "Compute element mapping for element replacement"

    # Extract all elements
    elements = extract_elements(form)
    elements = extract_sub_elements(elements)

    # Get cell and degree
    # FIXME: implement extract_common_top_domain(s) instead of this
    common_cell = extract_common_cell(form, common_cell)
    common_domain = as_domain(common_cell) # FIXME:
    common_degree = _auto_select_degree(elements)

    # Compute element map
    element_mapping = {}
    for element in elements:

        # Flag for whether element needs to be reconstructed
        reconstruct = False

        # Set cell
        domain = element.domain()
        if domain is None:
            info("Adjusting missing element domain to %s." % \
                     (common_domain,))
            domain = common_domain
            reconstruct = True

        # Set degree
        degree = element.degree()
        if degree is None:
            info("Adjusting element degree from %s to %d" % \
                     (istr(degree), common_degree))
            degree = common_degree
            reconstruct = True

        # Reconstruct element and add to map
        if reconstruct:
            element_mapping[element] = element.reconstruct(domain=domain,
                                                           degree=degree)

    return element_mapping

def check_swig_version(compiled_module):
    
    # Check swig version of compiled module
    if compiled_module and compiled_module.swigversion != ufc.__swigversion__:
        error("Incompatible swig versions detected. UFC swig "\
              "version is not the same as extension module swig "\
              "version: '%s' != '%s' " % \
              (ufc.__swigversion__, compiled_module.swigversion))

def jit_form(form, parameters=None, common_cell=None):
    "Just-in-time compile the given form."

    # Check that we get a Form
    if not isinstance(form, Form):
        form = as_form(form)

    # Check parameters
    parameters = _check_parameters(form, parameters)

    # Set log level
    set_level(parameters["log_level"])
    set_prefix(parameters["log_prefix"])

    # Compute element mapping for element replacement
    element_mapping = _compute_element_mapping(form, common_cell)

    # Compute form metadata and extract preprocessed form
    form_data = form.compute_form_data(common_cell=common_cell,
                                       element_mapping=element_mapping)

    # Wrap input
    jit_object = JITObject(form, parameters)

    # Set prefix for generated code
    module_name = "ffc_form_" + jit_object.signature()

    # Use Instant cache if possible
    cache_dir = parameters["cache_dir"] or None
    module = instant.import_module(module_name, cache_dir=cache_dir)
    if module:
        check_swig_version(module)
        debug("Reusing form from cache.")
        compiled_form = _extract_form(module, module_name)
        return (compiled_form, module, form_data, module_name)

    # Take lock to serialise file removal.
    # Need to add "_0" to lock as instant.import_module acquire
    # lock with name: module_name
    with instant.file_lock(instant.get_default_cache_dir(), \
                           module_name + "_0") as lock:

        # Retry Instant cache. The module may have been created while we waited
        # for the lock, even if it didn't exist before.
        module = instant.import_module(module_name, cache_dir=cache_dir)
        if module:
            check_swig_version(module)
            debug("Reusing form from cache.")
            compiled_form = _extract_form(module, module_name)
            return (compiled_form, module, form_data, module_name)

        # Write a message
        log(INFO + 5,
            "Calling FFC just-in-time (JIT) compiler, this may take some time.")

        # Generate code
        compile_form(form,
                     prefix=module_name,
                     parameters=parameters)

        # Build module using Instant (through UFC)
        debug("Compiling and linking Python extension module, this may take some time.")
        hfile   = module_name + ".h"
        cppfile = module_name + ".cpp"
        module = ufc_utils.build_ufc_module(
            hfile,
            source_directory = os.curdir,
            signature = module_name,
            sources = [cppfile] if parameters["split"] else [],
            cppargs = parameters["cpp_optimize_flags"].split() \
                      if parameters["cpp_optimize"] else ["-O0"],
            cache_dir = cache_dir)

        # Remove code
        if os.path.isfile(hfile):
            os.unlink(hfile)
        if parameters["split"] :
            if os.path.isfile(cppfile):
                os.unlink(cppfile)

        check_swig_version(module)

        # Extract compiled form
        compiled_form = _extract_form(module, module_name)

        return compiled_form, module, form_data, module_name

def jit_element(element, parameters=None):
    "Just-in-time compile the given element"

    # FIXME: We need a new solution for this.

    # Check that we get an element
    if not isinstance(element, FiniteElementBase):
        error("Expecting a finite element.")

    # Create simplest possible dummy form
    v = TestFunction(element)
    ii = (0,)*v.rank()
    form = v[ii]*dx

    # Compile form
    compiled_form, module, form_data, prefix = jit_form(form, parameters)

    return _extract_element_and_dofmap(module, prefix, form_data)

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

def _extract_form(module, prefix):
    "Extract form from module."
    return getattr(module, prefix + "_form_0")()

def _extract_element_and_dofmap(module, prefix, form_data):
    """
    Extract element and dofmap from module. Code will be generated for
    all unique elements (including sub elements) and to get the top
    level element we need to extract the last element.
    """
    i = len(form_data.unique_sub_elements) - 1
    return (getattr(module, prefix + ("_finite_element_%d" % i))(),
            getattr(module, prefix + ("_dofmap_%d" % i))())
