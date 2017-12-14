# -*- coding: utf-8 -*-
"""This module provides a just-in-time (JIT) form compiler.
It uses dijitso to wrap the generated code into a Python module."""

# Copyright (C) 2007-2017 Anders Logg
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
# Modified by Martin Sandve Alnæs, 2013-2017

# Python modules
import os
import sys
from hashlib import sha1

# FEniCS modules
import ufl

# Not importing globally to keep dijitso optional if jit is not used
#import dijitso

# FFC modules
from ffc.log import log
from ffc.log import error
from ffc.log import set_level
from ffc.log import set_prefix
from ffc.log import INFO
from ffc.parameters import validate_jit_parameters, compute_jit_parameters_signature
from ffc.compiler import compile_form, compile_element, compile_coordinate_mapping
from ffc.backends.ufc import get_include_path as get_ufc_include_path
from ffc.backends.ufc import get_ufc_signature, get_ufc_templates_signature
from ffc import __version__ as FFC_VERSION
from ffc.classname import make_classname


def jit_generate(ufl_object, module_name, signature, parameters):
    "Callback function passed to dijitso.jit: generate code and return as strings."
    log(INFO + 5, "Calling FFC just-in-time (JIT) compiler, this may take some time.")

    # Generate actual code for this object
    if isinstance(ufl_object, ufl.Form):
        compile_object = compile_form
    elif isinstance(ufl_object, ufl.FiniteElementBase):
        compile_object = compile_element
    elif isinstance(ufl_object, ufl.Mesh):
        compile_object = compile_coordinate_mapping

    code_h, code_c, dependent_ufl_objects = compile_object(ufl_object,
            prefix=module_name, parameters=parameters, jit=True)

    # Jit compile dependent objects separately,
    # but pass indirect=True to skip instantiating objects.
    # (this is done in here such that it's only triggered
    # if parent jit module is missing, and it's done after
    # compile_object because some misformed ufl objects may
    # require analysis to determine (looking at you Expression...))
    dependencies = []
    for dep in dependent_ufl_objects["element"]:
        dep_module_name = jit(dep, parameters, indirect=True)
        dependencies.append(dep_module_name)
    for dep in dependent_ufl_objects["coordinate_mapping"]:
        dep_module_name = jit(dep, parameters, indirect=True)
        dependencies.append(dep_module_name)
    return code_h, code_c, dependencies


def _string_tuple(param):
    "Split a : separated string or convert a list to a tuple."
    if isinstance(param, (tuple, list)):
        pass
    elif isinstance(param, str):
        param = param.split(":")
    else:
        param = ()
    param = tuple(p for p in param if p)
    assert all(isinstance(p, str) for p in param)
    return param


def jit_build(ufl_object, module_name, parameters):
    "Wraps dijitso jit with some parameter conversion etc."
    import dijitso

    # FIXME: Expose more dijitso parameters?
    # FIXME: dijitso build params are not part of module_name here.
    #        Currently dijitso doesn't add to the module signature.

    # Translating the C++ flags from ffc parameters to dijitso
    # to get equivalent behaviour to instant code
    build_params = {}
    build_params["debug"] = not parameters["cpp_optimize"]
    build_params["cxxflags_opt"] = tuple(parameters["cpp_optimize_flags"].split())
    build_params["cxxflags_debug"] = ("-O0",)
    build_params["include_dirs"] = (get_ufc_include_path(),) + _string_tuple(parameters.get("external_include_dirs"))
    build_params["lib_dirs"] = _string_tuple(parameters.get("external_library_dirs"))
    build_params["libs"] = _string_tuple(parameters.get("external_libraries"))

    # Interpreting FFC default "" as None, use "." if you want to point to curdir
    cache_dir = parameters.get("cache_dir") or None
    if cache_dir:
        cache_params = {"cache_dir": cache_dir}
    else:
        cache_params = {}

    # This will do some rudimenrary checking of the params and fill in dijitso defaults
    params = dijitso.validate_params({
        "cache": cache_params,
        "build": build_params,
        "generator": parameters,  # ffc parameters, just passed on to jit_generate
    })

    # Carry out jit compilation, calling jit_generate only if needed
    module, signature = dijitso.jit(jitable=ufl_object,
                                    name=module_name,
                                    params=params,
                                    generate=jit_generate)
    return module


def compute_jit_prefix(ufl_object, parameters, kind=None):
    "Compute the prefix (module name) for jit modules."

    # Get signature from ufl object
    if isinstance(ufl_object, ufl.Form):
        kind = "form"
        object_signature = ufl_object.signature()
    elif isinstance(ufl_object, ufl.Mesh):
        # When coordinate mapping is represented by a Mesh, just getting its coordinate element
        kind = "coordinate_mapping"
        ufl_object = ufl_object.ufl_coordinate_element()
        object_signature = repr(ufl_object)  # ** must match below
    elif kind == "coordinate_mapping" and isinstance(ufl_object, ufl.FiniteElementBase):
        # When coordinate mapping is represented by its coordinate element
        object_signature = repr(ufl_object)  # ** must match above
    elif isinstance(ufl_object, ufl.FiniteElementBase):
        kind = "element"
        object_signature = repr(ufl_object)
    else:
        error("Unknown ufl object type %s" % (ufl_object.__class__.__name__,))

    # Compute deterministic string of relevant parameters
    parameters_signature = compute_jit_parameters_signature(parameters)

    # Increase this number at any time to invalidate cache
    # signatures if code generation has changed in important
    # ways without the change being visible in regular signatures:
    jit_version_bump = 3

    # Build combined signature
    signatures = [
        object_signature,
        parameters_signature,
        str(FFC_VERSION),
        str(jit_version_bump),
        get_ufc_signature(),
        get_ufc_templates_signature(),
        kind,
        ]
    string = ";".join(signatures)
    signature = sha1(string.encode('utf-8')).hexdigest()

    # Optionally shorten signature
    max_signature_length = parameters["max_signature_length"]
    if max_signature_length:
        signature = signature[:max_signature_length]

    # Combine into prefix with some info including kind
    prefix = ("ffc_%s_%s" % (kind, signature)).lower()
    return kind, prefix


class FFCError(Exception):
    pass


class FFCJitError(FFCError):
    pass


def jit(ufl_object, parameters=None, indirect=False):
    """Just-in-time compile the given form or element

    Parameters:

      ufl_object : The UFL object to be compiled
      parameters : A set of parameters
    """
    # Check parameters
    parameters = validate_jit_parameters(parameters)

    # FIXME: Setting the log level here becomes a permanent side effect...
    # Set log level
    set_level(parameters["log_level"])
    set_prefix(parameters["log_prefix"])

    # Make unique module name for generated code
    kind, module_name = compute_jit_prefix(ufl_object, parameters)

    # Inspect cache and generate+build if necessary
    module = jit_build(ufl_object, module_name, parameters)

    # Raise exception on failure to build or import module
    if module is None:
        # TODO: To communicate directory name here, need dijitso params to call
        #fail_dir = dijitso.cache.create_fail_dir_path(signature, dijitso_cache_params)
        raise FFCJitError("A directory with files to reproduce the jit build failure has been created.")

    # Construct instance of object from compiled code unless indirect
    if indirect:
        return module_name
    else:
        # FIXME: Streamline number of return arguments here across kinds
        if kind == "form":
            if parameters.get("representation") == "quadrature" or \
                    any(itg.metadata().get("representation") == "quadrature" for itg in ufl_object.integrals()):
                from ffc.quadrature.deprecation import issue_deprecation_warning
                issue_deprecation_warning()

            compiled_form = _instantiate_form(module, module_name)
            return (compiled_form, module, module_name)
            # TODO: module, module_name are never used in dolfin, drop?
            #return _instantiate_form(module, module_name)
        elif kind == "element":
            fe, dm = _instantiate_element_and_dofmap(module, module_name)
            return fe, dm
        elif kind == "coordinate_mapping":
            cm = _instantiate_coordinate_mapping(module, module_name)
            return cm
        else:
            error("Unknown kind %s" % (kind,))


def _instantiate_form(module, prefix):
    "Instantiate an object of the jit-compiled form."
    import dijitso
    classname = make_classname(prefix, "form", "main")
    form = dijitso.extract_factory_function(module, "create_" + classname)()
    return form


def _instantiate_element_and_dofmap(module, prefix):
    "Instantiate objects of the jit-compiled finite_element and dofmap."
    import dijitso
    fe_classname = make_classname(prefix, "finite_element", "main")
    dm_classname = make_classname(prefix, "dofmap", "main")
    fe = dijitso.extract_factory_function(module, "create_" + fe_classname)()
    dm = dijitso.extract_factory_function(module, "create_" + dm_classname)()
    return (fe, dm)


def _instantiate_coordinate_mapping(module, prefix):
    "Instantiate an object of the jit-compiled coordinate_mapping."
    import dijitso
    classname = make_classname(prefix, "coordinate_mapping", "main")
    form = dijitso.extract_factory_function(module, "create_" + classname)()
    return form
