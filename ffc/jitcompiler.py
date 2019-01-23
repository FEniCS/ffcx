# -*- coding: utf-8 -*-
# Copyright (C) 2007-2017 Anders Logg
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Just-in-time (JIT) form compiler. It uses dijitso to wrap the generated
code into a Python module.

"""

import hashlib
import logging
import os

import dijitso
import ufl
from ffc import __version__ as FFC_VERSION
from ffc import classname
from ffc.codegeneration import get_include_path, get_signature
from ffc import compiler
from ffc.parameters import (compute_jit_parameters_signature, validate_jit_parameters)

logger = logging.getLogger(__name__)


def generate(ufl_object, module_name, signature, parameters):
    """Callback function passed to dijitso.jit: generate code and return as strings."""
    logger.info("Calling FFC just-in-time (JIT) compiler.")

    # Pick the generator for actual code for this object
    if isinstance(ufl_object, ufl.Form):
        compile_object = compiler.compile_form
    elif isinstance(ufl_object, ufl.FiniteElementBase):
        compile_object = compiler.compile_element
    elif isinstance(ufl_object, ufl.Mesh):
        compile_object = compiler.compile_coordinate_mapping

    # Return C code for requested ufl_object, and return any UFL objects
    # that ufl_objects needs, e.g. a form will require some elements.
    code_h, code_c, dependent_ufl_objects = compile_object(
        ufl_object, prefix=module_name, parameters=parameters, jit=True)

    # Jit compile dependent objects separately, but pass indirect=True
    # to skip instantiating objects. (This is done in here such that
    # it's only triggered if parent jit module is missing, and it's done
    # after compile_object because some misformed ufl objects may
    # require analysis to determine (looking at you Expression...)).
    # Dependencies is the name (string) of the compiled shared library.
    dependencies = []
    for dep in dependent_ufl_objects["element"]:
        dep_module_name = jit(dep, parameters, indirect=True)
        dependencies.append(dep_module_name)
    for dep in dependent_ufl_objects["coordinate_mapping"]:
        dep_module_name = jit(dep, parameters, indirect=True)
        dependencies.append(dep_module_name)

    return code_h, code_c, dependencies


def _string_tuple(param):
    """Split a : separated string or convert a list to a tuple."""
    if isinstance(param, (tuple, list)):
        pass
    elif isinstance(param, str):
        param = param.split(":")
    else:
        param = ()
    param = tuple(p for p in param if p)
    assert all(isinstance(p, str) for p in param)
    return param


def build(ufl_object, module_name, parameters):
    """Wraps dijitso jit with some parameter conversion etc."""
    # FIXME: Expose more dijitso parameters?
    # FIXME: dijitso build params are not part of module_name here.
    #        Currently dijitso doesn't add to the module signature.

    # Translating the C++ flags from ffc parameters to dijitso to get
    # equivalent behaviour to instant code
    build_params = {}
    build_params["debug"] = not parameters["cpp_optimize"]
    build_params["cxxflags_opt"] = tuple(parameters["cpp_optimize_flags"].split())
    build_params["cxxflags_debug"] = ("-O0", )
    build_params["include_dirs"] = (get_include_path(), ) + _string_tuple(
        parameters.get("external_include_dirs"))
    build_params["lib_dirs"] = _string_tuple(parameters.get("external_library_dirs"))

    # Use C compiler (not C++) and add some libs/options
    build_params["cxx"] = os.getenv("CC", "cc")
    build_params["libs"] = _string_tuple("m" + parameters.get("external_libraries"))
    build_params["cxxflags"] = ["-Wall", "-shared", "-fPIC"]

    # Interpreting FFC default "" as None, use "." if you want to point to curdir
    cache_dir = parameters.get("cache_dir") or None
    if cache_dir:
        cache_params = {"cache_dir": cache_dir}
    else:
        cache_params = {}

    # Use standard C file extension
    cache_params["src_postfix"] = ".c"

    # This will do some rudimentary checking of the params and fill in
    # dijitso defaults
    params = dijitso.validate_params({
        "cache": cache_params,
        "build": build_params,
        "generator": parameters,  # ffc parameters, just passed on to generate
    })

    # Carry out jit compilation, calling generate only if needed
    module, signature = dijitso.jit(
        jitable=ufl_object, name=module_name, params=params, generate=generate)

    return module


def compute_prefix(ufl_object, tag, parameters, kind=None):
    """Compute the prefix (module name) for jit modules."""

    # Get signature from ufl object
    if isinstance(ufl_object, ufl.Form):
        kind = "form"
        object_signature = ufl_object.signature()
    elif isinstance(ufl_object, ufl.Mesh):
        # When coordinate mapping is represented by a Mesh, just getting
        # its coordinate element
        kind = "coordinate_mapping"
        ufl_object = ufl_object.ufl_coordinate_element()
        object_signature = repr(ufl_object)  # ** must match below
    elif kind == "coordinate_mapping" and isinstance(ufl_object, ufl.FiniteElementBase):
        # When coordinate mapping is represented by its coordinate
        # element
        object_signature = repr(ufl_object)  # ** must match above
    elif isinstance(ufl_object, ufl.FiniteElementBase):
        kind = "element"
        object_signature = repr(ufl_object)
    else:
        raise RuntimeError("Unknown ufl object type {}".format(ufl_object.__class__.__name__))

    # Compute deterministic string of relevant parameters
    parameters_signature = compute_jit_parameters_signature(parameters)

    # Increase this number at any time to invalidate cache signatures if
    # code generation has changed in important ways without the change
    # being visible in regular signatures:
    jit_version_bump = 3

    # Build combined signature
    signatures = [
        object_signature,
        parameters_signature,
        str(FFC_VERSION),
        str(jit_version_bump),
        str(tag),
        get_signature(),
        kind,
    ]
    string = ";".join(signatures)
    signature = hashlib.sha1(string.encode('utf-8')).hexdigest()

    # Combine into prefix with some info including kind
    prefix = "ffc_{}_{}".format(kind, signature).lower()
    return kind, prefix


def jit(ufl_object, parameters=None, indirect=False):
    """Just-in-time compile the given form or element

    Parameters
    ----------
      ufl_object : The UFL object to be compiled
      parameters : A set of parameters

    """
    # Check parameters
    parameters = validate_jit_parameters(parameters)

    # Make unique module name for generated code
    kind, module_name = compute_prefix(ufl_object, None, parameters)

    # Get module (inspect cache and generate+build if necessary)
    module = build(ufl_object, module_name, parameters)

    # Raise exception on failure to build or import module
    if module is None:
        # TODO: To communicate directory name here, need dijitso params
        # to call
        # fail_dir = dijitso.cache.create_fail_dir_path(signature, dijitso_cache_params)
        raise RuntimeError(
            "A directory with files to reproduce the jit build failure has been created.")

    # Construct instance of object from compiled code, unless indirect
    # in which case return the name
    if indirect:
        return module_name
    else:
        # FIXME: Streamline number of return arguments here across kinds
        if kind == "form":
            compiled_form = _instantiate_form(module, module_name)
            return compiled_form, module, module_name
            # TODO: module, module_name are never used in dolfin, drop?
            # return _instantiate_form(module, module_name)
        elif kind == "element":
            fe, dm = _instantiate_element_and_dofmap(module, module_name)
            return fe, dm
        elif kind == "coordinate_mapping":
            cm = _instantiate_coordinate_mapping(module, module_name)
            return cm
        else:
            raise RuntimeError("Unknown kind {}".format(kind))


def _instantiate_form(module, prefix):
    """Instantiate an object of the jit-compiled form."""
    name = classname.make_name(prefix, "form", "main")
    form = dijitso.extract_factory_function(module, "create_" + name)()
    return form


def _instantiate_element_and_dofmap(module, prefix):
    """Instantiate objects of the jit-compiled finite_element and dofmap."""
    fe_classname = classname.make_name(prefix, "finite_element", "main")
    dm_classname = classname.make_name(prefix, "dofmap", "main")
    fe = dijitso.extract_factory_function(module, "create_" + fe_classname)()
    dm = dijitso.extract_factory_function(module, "create_" + dm_classname)()
    return (fe, dm)


def _instantiate_coordinate_mapping(module, prefix):
    """Instantiate an object of the jit-compiled coordinate_mapping."""
    name = classname.make_name(prefix, "coordinate_mapping", "main")
    form = dijitso.extract_factory_function(module, "create_" + name)()
    return form
