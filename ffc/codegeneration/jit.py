# Copyright (C) 2004-2019 Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import importlib
import logging
import os
import re
import tempfile
import time
from pathlib import Path

import cffi

import ffc

logger = logging.getLogger(__name__)

UFC_HEADER_DECL = """
typedef {} ufc_scalar_t;  /* Hack to deal with scalar type */

typedef struct ufc_coordinate_mapping ufc_coordinate_mapping;
typedef struct ufc_finite_element ufc_finite_element;
typedef struct ufc_dofmap ufc_dofmap;

typedef enum
{{
interval = 10,
triangle = 20,
quadrilateral = 30,
tetrahedron = 40,
hexahedron = 50,
vertex = 60,
}} ufc_shape;
"""

# Get declarations directly from ufc.h
file_dir = os.path.dirname(os.path.abspath(__file__))
with open(file_dir + "/ufc.h", "r") as f:
    ufc_h = ''.join(f.readlines())

UFC_ELEMENT_DECL = '\n'.join(re.findall('typedef struct ufc_finite_element.*?ufc_finite_element;', ufc_h, re.DOTALL))
UFC_DOFMAP_DECL = '\n'.join(re.findall('typedef struct ufc_dofmap.*?ufc_dofmap;', ufc_h, re.DOTALL))
UFC_COORDINATEMAPPING_DECL = '\n'.join(re.findall('typedef struct ufc_coordinate_mapping.*?ufc_coordinate_mapping;',
                                                  ufc_h, re.DOTALL))
UFC_FORM_DECL = '\n'.join(re.findall('typedef struct ufc_form.*?ufc_form;', ufc_h, re.DOTALL))

UFC_INTEGRAL_DECL = '\n'.join(re.findall(r'typedef void \(ufc_tabulate_tensor\).*?\);', ufc_h, re.DOTALL))
UFC_INTEGRAL_DECL += '\n'.join(re.findall(r'typedef void \(ufc_tabulate_tensor_custom\).*?\);', ufc_h, re.DOTALL))
UFC_INTEGRAL_DECL += '\n'.join(re.findall('typedef struct ufc_integral.*?ufc_integral;',
                                          ufc_h, re.DOTALL))
UFC_INTEGRAL_DECL += '\n'.join(re.findall('typedef struct ufc_custom_integral.*?ufc_custom_integral;',
                                          ufc_h, re.DOTALL))
UFC_EXPRESSION_DECL = '\n'.join(re.findall('typedef struct ufc_expression.*?ufc_expression;', ufc_h, re.DOTALL))


def _compute_parameter_signature(parameters):
    """Return parameters signature (some parameters should not affect signature)."""
    return str(sorted(parameters.items()))


def get_cached_module(module_name, object_names, cache_dir, timeout):
    """Look for an existing C file and wait for compilation, or if it does not exist, create it."""
    cache_dir = Path(cache_dir)
    c_filename = cache_dir.joinpath(module_name).with_suffix(".c")
    ready_name = c_filename.with_suffix(".c.cached")

    # Ensure cache dir exists
    cache_dir.mkdir(exist_ok=True, parents=True)

    try:
        # Create C file with exclusive access
        open(c_filename, "x")
        return None, None
    except FileExistsError:
        logger.info("Cached C file already exists: " + str(c_filename))
        finder = importlib.machinery.FileFinder(
            str(cache_dir), (importlib.machinery.ExtensionFileLoader, importlib.machinery.EXTENSION_SUFFIXES))
        finder.invalidate_caches()

        # Now, wait for ready
        for i in range(timeout):
            if os.path.exists(ready_name):
                spec = finder.find_spec(module_name)
                if spec is None:
                    raise ModuleNotFoundError("Unable to find JIT module.")
                compiled_module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(compiled_module)

                compiled_objects = [getattr(compiled_module.lib, "create_" + name)() for name in object_names]
                return compiled_objects, compiled_module

            logger.info("Waiting for {} to appear.".format(str(ready_name)))
            time.sleep(1)
        raise TimeoutError("""JIT compilation timed out, probably due to a failed previous compile.
        Try cleaning cache (e.g. remove {}) or increase timeout parameter.""".format(c_filename))


def compile_elements(elements, parameters=None, cache_dir=None, timeout=10, cffi_extra_compile_args=None,
                     cffi_verbose=False, cffi_debug=None):
    """Compile a list of UFL elements and dofmaps into Python objects."""
    p = ffc.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    logger.info('Compiling elements: ' + str(elements))

    # Get a signature for these elements
    module_name = 'libffc_elements_' + \
        ffc.classname.compute_signature(elements, _compute_parameter_signature(p)
                                        + str(cffi_extra_compile_args) + str(cffi_debug))

    names = []
    for e in elements:
        name = ffc.ir.representation.make_finite_element_classname(e, "JIT")
        names.append(name)
        name = ffc.ir.representation.make_dofmap_classname(e, "JIT")
        names.append(name)

    if cache_dir is not None:
        obj, mod = get_cached_module(module_name, names, cache_dir, timeout)
        if obj is not None:
            # Pair up elements with dofmaps
            obj = list(zip(obj[::2], obj[1::2]))
            return obj, mod

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL
    element_template = "ufc_finite_element * create_{name}(void);\n"
    dofmap_template = "ufc_dofmap * create_{name}(void);\n"
    for i in range(len(elements)):
        decl += element_template.format(name=names[i * 2])
        decl += dofmap_template.format(name=names[i * 2 + 1])

    objects, module = _compile_objects(decl, elements, names, module_name, p, cache_dir,
                                       cffi_extra_compile_args, cffi_verbose, cffi_debug)
    # Pair up elements with dofmaps
    objects = list(zip(objects[::2], objects[1::2]))
    return objects, module


def compile_forms(forms, parameters=None, cache_dir=None, timeout=10, cffi_extra_compile_args=None,
                  cffi_verbose=False, cffi_debug=None):
    """Compile a list of UFL forms into UFC Python objects."""
    p = ffc.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    logger.info('Compiling forms: ' + str(forms))

    # Get a signature for these forms
    module_name = 'libffc_forms_' + \
        ffc.classname.compute_signature(forms, _compute_parameter_signature(p)
                                        + str(cffi_extra_compile_args) + str(cffi_debug))

    form_names = [ffc.classname.make_name("JIT", "form", ffc.classname.compute_signature([form], str(i)))
                  for i, form in enumerate(forms)]

    if cache_dir is not None:
        obj, mod = get_cached_module(module_name, form_names, cache_dir, timeout)
        if obj is not None:
            return obj, mod

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL + \
        UFC_COORDINATEMAPPING_DECL + UFC_INTEGRAL_DECL + UFC_FORM_DECL

    form_template = "ufc_form * create_{name}(void);\n"
    for name in form_names:
        decl += form_template.format(name=name)

    return _compile_objects(decl, forms, form_names, module_name, p, cache_dir,
                            cffi_extra_compile_args, cffi_verbose, cffi_debug)


def compile_expressions(expressions, parameters=None, cache_dir=None, timeout=10, cffi_extra_compile_args=None,
                        cffi_verbose=False, cffi_debug=None):
    """Compile a list of UFL expressions into UFC Python objects.

    Parameters
    ----------
    expressions
        List of (UFL expression, evaluation points).

    """
    p = ffc.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    logger.info('Compiling expressions: ' + str(expressions))

    # Get a signature for these forms
    module_name = 'libffc_expressions_' + ffc.classname.compute_signature(expressions, '', p)

    expr_names = [ffc.classname.make_name("JIT", "expression", ffc.classname.compute_signature([expression], "", p))
                  for expression in expressions]

    if cache_dir is not None:
        obj, mod = get_cached_module(module_name, expr_names, cache_dir, timeout)
        if obj is not None:
            return obj, mod

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL + \
        UFC_COORDINATEMAPPING_DECL + UFC_INTEGRAL_DECL + UFC_FORM_DECL + UFC_EXPRESSION_DECL

    expression_template = "ufc_expression* create_{name}(void);\n"
    for name in expr_names:
        decl += expression_template.format(name=name)

    return _compile_objects(decl, expressions, expr_names, module_name, p, cache_dir,
                            cffi_extra_compile_args, cffi_verbose, cffi_debug)


def compile_coordinate_maps(meshes, parameters=None, cache_dir=None, timeout=10, cffi_extra_compile_args=None,
                            cffi_verbose=False, cffi_debug=None):
    """Compile a list of UFL coordinate mappings into UFC Python objects."""
    p = ffc.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    logger.info('Compiling cmaps: ' + str(meshes))

    # Get a signature for these cmaps
    module_name = 'libffc_cmaps_' + \
        ffc.classname.compute_signature(meshes, _compute_parameter_signature(
            p) + str(cffi_extra_compile_args) + str(cffi_debug), True)

    cmap_names = [ffc.ir.representation.make_coordinate_map_classname(
        mesh.ufl_coordinate_element(), "JIT") for mesh in meshes]

    if cache_dir is not None:
        obj, mod = get_cached_module(module_name, cmap_names, cache_dir, timeout)
        if obj is not None:
            return obj, mod

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_COORDINATEMAPPING_DECL
    cmap_template = "ufc_coordinate_mapping * create_{name}(void);\n"

    for name in cmap_names:
        decl += cmap_template.format(name=name)

    return _compile_objects(decl, meshes, cmap_names, module_name, p, cache_dir,
                            cffi_extra_compile_args, cffi_verbose, cffi_debug)


def _compile_objects(decl, ufl_objects, object_names, module_name, parameters, cache_dir,
                     cffi_extra_compile_args, cffi_verbose, cffi_debug):
    if cache_dir is None:
        cache_dir = Path(tempfile.mkdtemp())
    else:
        cache_dir = Path(cache_dir)
    _, code_body = ffc.compiler.compile_ufl_objects(ufl_objects, prefix="JIT", parameters=parameters)

    ffibuilder = cffi.FFI()
    ffibuilder.set_source(module_name, code_body, include_dirs=[ffc.codegeneration.get_include_path()],
                          extra_compile_args=cffi_extra_compile_args)
    ffibuilder.cdef(decl)

    c_filename = cache_dir.joinpath(module_name + ".c")
    ready_name = c_filename.with_suffix(".c.cached")

    # Compile (ensuring that compile dir exists)
    cache_dir.mkdir(exist_ok=True, parents=True)

    try:
        ffibuilder.compile(tmpdir=cache_dir, verbose=cffi_verbose, debug=cffi_debug)
    except Exception:
        os.replace(c_filename, c_filename.with_suffix(".c.failed"))
        raise

    # Create a "status ready" file. If this fails, it is an error,
    # because it should not exist yet.
    fd = open(ready_name, "x")
    fd.close()

    # Create module finder that searches the compile path
    finder = importlib.machinery.FileFinder(
        str(cache_dir), (importlib.machinery.ExtensionFileLoader, importlib.machinery.EXTENSION_SUFFIXES))

    # Find module. Clear search cache to be sure dynamically created
    # (new) modules are found
    finder.invalidate_caches()
    spec = finder.find_spec(module_name)
    if spec is None:
        raise ModuleNotFoundError("Unable to find JIT module.")

    # Load module
    compiled_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(compiled_module)

    compiled_objects = [getattr(compiled_module.lib, "create_" + name)() for name in object_names]

    return compiled_objects, compiled_module
