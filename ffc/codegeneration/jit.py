# -*- coding: utf-8 -*-
# Copyright (C) 2018 Garth N. Wells
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
import ffc.config

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


def get_cached_module(module_name, object_names, parameters):
    """Look for an existing C file and wait for compilation, or if it does not exist, create it."""
    cache_dir = ffc.config.get_cache_path(parameters)
    timeout = int(parameters.get("timeout", 10))
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


def compile_elements(elements, parameters=None):
    """Compile a list of UFL elements and dofmaps into UFC Python objects."""
    p = ffc.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    logger.info('Compiling elements: ' + str(elements))

    # Get a signature for these elements
    module_name = 'libffc_elements_' + ffc.classname.compute_signature(elements, '', p)

    names = []
    for e in elements:
        name = ffc.ir.representation.make_finite_element_jit_classname(e, "JIT", p)
        names.append(name)
        name = ffc.ir.representation.make_dofmap_jit_classname(e, "JIT", p)
        names.append(name)

    if p['use_cache']:
        obj, mod = get_cached_module(module_name, names, p)
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

    objects, module = _compile_objects(decl, elements, names, module_name, p)
    # Pair up elements with dofmaps
    objects = list(zip(objects[::2], objects[1::2]))
    return objects, module


def compile_forms(forms, parameters=None):
    """Compile a list of UFL forms into UFC Python objects."""
    p = ffc.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    logger.info('Compiling forms: ' + str(forms))

    # Get a signature for these forms
    module_name = 'libffc_forms_' + ffc.classname.compute_signature(forms, '', p)

    form_names = [ffc.classname.make_name("JIT", "form", ffc.classname.compute_signature([form], str(i), p))
                  for i, form in enumerate(forms)]

    if p['use_cache']:
        obj, mod = get_cached_module(module_name, form_names, p)
        if obj is not None:
            return obj, mod

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL + \
        UFC_COORDINATEMAPPING_DECL + UFC_INTEGRAL_DECL + UFC_FORM_DECL

    form_template = "ufc_form * create_{name}(void);\n"
    for name in form_names:
        decl += form_template.format(name=name)

    return _compile_objects(decl, forms, form_names, module_name, p)


def compile_coordinate_maps(meshes, parameters=None):
    """Compile a list of UFL coordinate mappings into UFC Python objects."""
    p = ffc.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    logger.info('Compiling cmaps: ' + str(meshes))

    # Get a signature for these cmaps
    module_name = 'libffc_cmaps_' + ffc.classname.compute_signature(meshes, '', p, True)

    cmap_names = [ffc.ir.representation.make_coordinate_mapping_jit_classname(
        mesh.ufl_coordinate_element(), "JIT", p) for mesh in meshes]

    if p['use_cache']:
        obj, mod = get_cached_module(module_name, cmap_names, p)
        if obj is not None:
            return obj, mod

    scalar_type = p["scalar_type"].replace("complex", "_Complex")
    decl = UFC_HEADER_DECL.format(scalar_type) + UFC_COORDINATEMAPPING_DECL
    cmap_template = "ufc_coordinate_mapping * create_{name}(void);\n"

    for name in cmap_names:
        decl += cmap_template.format(name=name)

    return _compile_objects(decl, meshes, cmap_names, module_name, p)


def _compile_objects(decl, ufl_objects, object_names, module_name, parameters):
    if (parameters['use_cache']):
        compile_dir = ffc.config.get_cache_path(parameters)
    else:
        compile_dir = Path(tempfile.mkdtemp())
    _, code_body = ffc.compiler.compile_ufl_objects(ufl_objects, prefix="JIT", parameters=parameters)

    ffibuilder = cffi.FFI()
    ffibuilder.set_source(module_name, code_body, include_dirs=[ffc.codegeneration.get_include_path()],
                          extra_compile_args=['-g0'])  # turn off -g
    ffibuilder.cdef(decl)

    c_filename = compile_dir.joinpath(module_name + ".c")
    ready_name = c_filename.with_suffix(".c.cached")

    # Compile (ensuring that compile dir exists)
    compile_dir.mkdir(exist_ok=True, parents=True)
    ffibuilder.compile(tmpdir=compile_dir, verbose=False)

    # Create a "status ready" file. If this fails, it is an error,
    # because it should not exist yet.
    fd = open(ready_name, "x")
    fd.close()

    # Create module finder that searches the compile path
    finder = importlib.machinery.FileFinder(
        str(compile_dir), (importlib.machinery.ExtensionFileLoader, importlib.machinery.EXTENSION_SUFFIXES))

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
