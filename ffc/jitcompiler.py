"""This module provides a just-in-time (JIT) form compiler.
It uses dijitso to wrap the generated code into a Python module."""

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
import os
import sys
import dijitso

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
from ffc.parameters import validate_jit_parameters
from ffc.mixedelement import MixedElement
from ffc.compiler import compile_form, compile_element
from ffc.formatting import write_code
from ffc.jitobject import JITObject
from ffc.quadratureelement import default_quadrature_degree
from ffc.ufc_config import get_ufc_include


def jit_generate(ufl_object, module_name, parameters):
    "Generate code and return as strings."
    if isinstance(ufl_object, Form):
        code_h, code_c = compile_form(ufl_object, prefix=module_name,
                                      parameters=parameters, jit=True)
    elif isinstance(ufl_object, FiniteElementBase):
        code_h, code_c = compile_element(ufl_object, prefix=module_name,
                                         parameters=parameters, jit=True)
    return code_h, code_c


def jit_build_with_dijitso(ufl_object, module_name, parameters):

    def _generate(ufl_object, module_name, signature, parameters):
        # Write a message
        log(INFO + 5,
            "Calling FFC just-in-time (JIT) compiler, this may take some time.")
        code_h, code_c = jit_generate(ufl_object, module_name, parameters)
        dependencies = ()
        return code_h, code_c, dependencies

    # Translating the C++ flags from ffc parameters to dijitso
    # to get equivalent behaviour to instant code
    if parameters["cpp_optimize"]:
        build_params = {
            "cxxflags_opt": tuple(parameters["cpp_optimize_flags"].split()),
            "debug": False
        }
    else:
        build_params = {
            "cxxflags_debug": ("-O0",),
            "debug": True
        }

    # Add path to UFC include dir
    build_params["include_dirs"] = get_ufc_include()

    # FFC default is "", use "." to point to curdir
    cache_dir = parameters.get("cache_dir") or None
    if cache_dir:
        cache_params = {"cache_dir": cache_dir}
    else:
        cache_params = {}

    params = dijitso.validate_params({
        "cache": cache_params,
        "build": build_params,
        "generator": parameters,
    })

    module, signature = dijitso.jit(ufl_object, module_name, params, _generate)
    return module


def jit(ufl_object, parameters=None):
    """Just-in-time compile the given form or element

    Parameters:

      ufl_object : The UFL object to be compiled
      parameters : A set of parameters
    """
    # Check that we get a form or element
    if isinstance(ufl_object, Form):
        kind = "form"
    elif isinstance(ufl_object, FiniteElementBase):
        kind = "element"
    else:
        error("Expecting a UFL form or element, got: %s" % repr(ufl_object))

    # Check parameters
    parameters = validate_jit_parameters(parameters)

    # FIXME: Setting the log level here becomes a permanent side effect...
    # Set log level
    set_level(parameters["log_level"])
    set_prefix(parameters["log_prefix"])

    # Wrap input
    jit_object = JITObject(ufl_object, parameters)

    # Set prefix for generated code
    module_name = "ffc_%s_%s" % (kind, jit_object.signature())

    # Inspect cache and generate+build if necessary
    module = jit_build_with_dijitso(ufl_object, module_name, parameters)

    # Construct instance of compiled form
    if isinstance(ufl_object, Form):
        compiled_form = _instantiate_form(module, module_name)
        return compiled_form, module, module_name
    elif isinstance(ufl_object, FiniteElementBase):
        return _instantiate_element_and_dofmap(module, module_name)


from ffc.cpp import make_classname


def _instantiate_form(module, prefix):
    "Extract the form from module with only one form."
    form_id = 0
    classname = make_classname(prefix, "form", form_id)
    form = dijitso.extract_factory_function(module, "create_" + classname)()
    return form


def _instantiate_element_and_dofmap(module, prefix):
    """Extract element and dofmap from module."""
    fe_classname = make_classname(prefix, "finite_element", "main")
    dm_classname = make_classname(prefix, "dofmap", "main")
    fe = dijitso.extract_factory_function(module, "create_" + fe_classname)()
    dm = dijitso.extract_factory_function(module, "create_" + dm_classname)()

    return (fe, dm)
