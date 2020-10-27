# Copyright (C) 2004-2019 Garth N. Wells
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from contextlib import redirect_stdout
import importlib
import io
import logging
import os
import re
import tempfile
import time
from pathlib import Path

import cffi
import ffcx
import ffcx.naming

logger = logging.getLogger("ffcx")

# Get declarations directly from ufc.h
file_dir = os.path.dirname(os.path.abspath(__file__))
with open(file_dir + "/ufc.h", "r") as f:
    ufc_h = ''.join(f.readlines())

UFC_HEADER_DECL = "typedef {} ufc_scalar_t;  /* Hack to deal with scalar type */\n"
header = ufc_h.split("<HEADER_DECL>")[1].split("</HEADER_DECL>")[0].strip(" /\n")
header = header.replace("{", "{{").replace("}", "}}")
UFC_HEADER_DECL += header + "\n"
UFC_HEADER_DECL += "void free(void *); \n"

UFC_ELEMENT_DECL = '\n'.join(re.findall('typedef struct ufc_finite_element.*?ufc_finite_element;', ufc_h, re.DOTALL))
UFC_DOFMAP_DECL = '\n'.join(re.findall('typedef struct ufc_dofmap.*?ufc_dofmap;', ufc_h, re.DOTALL))
UFC_COORDINATEMAPPING_DECL = '\n'.join(re.findall('typedef struct ufc_coordinate_mapping.*?ufc_coordinate_mapping;',
                                                  ufc_h, re.DOTALL))
UFC_FORM_DECL = '\n'.join(re.findall('typedef struct ufc_form.*?ufc_form;', ufc_h, re.DOTALL))

UFC_INTEGRAL_DECL = '\n'.join(re.findall(r'typedef void ?\(ufc_tabulate_tensor\).*?\);', ufc_h, re.DOTALL))
UFC_INTEGRAL_DECL += '\n'.join(re.findall(r'typedef void ?\(ufc_tabulate_tensor_custom\).*?\);', ufc_h, re.DOTALL))
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
                     cffi_verbose=False, cffi_debug=None, cffi_libraries=None):
    """Compile a list of UFL elements and dofmaps into Python objects."""
    p = ffcx.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    p.update(ffcx.parameters.env_parameters())

    # Get a signature for these elements
    module_name = 'libffcx_elements_' + \
        ffcx.naming.compute_signature(elements, _compute_parameter_signature(p)
                                      + str(cffi_extra_compile_args) + str(cffi_debug))

    names = []
    for e in elements:
        name = ffcx.naming.finite_element_name(e, "JIT")
        names.append(name)
        name = ffcx.naming.dofmap_name(e, "JIT")
        names.append(name)

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        obj, mod = get_cached_module(module_name, names, cache_dir, timeout)
        if obj is not None:
            # Pair up elements with dofmaps
            obj = list(zip(obj[::2], obj[1::2]))
            return obj, mod
    else:
        cache_dir = Path(tempfile.mkdtemp())

    try:
        scalar_type = p["scalar_type"].replace("complex", "_Complex")
        decl = UFC_HEADER_DECL.format(scalar_type) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL
        element_template = "ufc_finite_element * create_{name}(void);\n"
        dofmap_template = "ufc_dofmap * create_{name}(void);\n"
        for i in range(len(elements)):
            decl += element_template.format(name=names[i * 2])
            decl += dofmap_template.format(name=names[i * 2 + 1])

        _compile_objects(decl, elements, names, module_name, p, cache_dir,
                         cffi_extra_compile_args, cffi_verbose, cffi_debug, cffi_libraries)
    except Exception:
        # remove c file so that it will not timeout next time
        c_filename = cache_dir.joinpath(module_name + ".c")
        os.replace(c_filename, c_filename.with_suffix(".c.failed"))
        raise

    objects, module = _load_objects(cache_dir, module_name, names)
    # Pair up elements with dofmaps
    objects = list(zip(objects[::2], objects[1::2]))
    return objects, module


def compile_forms(forms, parameters=None, cache_dir=None, timeout=10, cffi_extra_compile_args=None,
                  cffi_verbose=False, cffi_debug=None, cffi_libraries=None):
    """Compile a list of UFL forms into UFC Python objects."""
    p = ffcx.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    p.update(ffcx.parameters.env_parameters())

    # Get a signature for these forms
    module_name = 'libffcx_forms_' + \
        ffcx.naming.compute_signature(forms, _compute_parameter_signature(p)
                                      + str(cffi_extra_compile_args) + str(cffi_debug))

    form_names = [ffcx.naming.form_name(form, i) for i, form in enumerate(forms)]

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        obj, mod = get_cached_module(module_name, form_names, cache_dir, timeout)
        if obj is not None:
            return obj, mod
    else:
        cache_dir = Path(tempfile.mkdtemp())

    try:
        scalar_type = p["scalar_type"].replace("complex", "_Complex")
        decl = UFC_HEADER_DECL.format(scalar_type) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL + \
            UFC_COORDINATEMAPPING_DECL + UFC_INTEGRAL_DECL + UFC_FORM_DECL

        form_template = "ufc_form * create_{name}(void);\n"
        for name in form_names:
            decl += form_template.format(name=name)

        _compile_objects(decl, forms, form_names, module_name, p, cache_dir,
                         cffi_extra_compile_args, cffi_verbose, cffi_debug, cffi_libraries)
    except Exception:
        # remove c file so that it will not timeout next time
        c_filename = cache_dir.joinpath(module_name + ".c")
        os.replace(c_filename, c_filename.with_suffix(".c.failed"))
        raise

    obj, module = _load_objects(cache_dir, module_name, form_names)
    return obj, module


def compile_expressions(expressions, parameters=None, cache_dir=None, timeout=10, cffi_extra_compile_args=None,
                        cffi_verbose=False, cffi_debug=None, cffi_libraries=None):
    """Compile a list of UFL expressions into UFC Python objects.

    Parameters
    ----------
    expressions
        List of (UFL expression, evaluation points).

    """
    p = ffcx.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    p.update(ffcx.parameters.env_parameters())

    # Get a signature for these forms
    module_name = 'libffcx_expressions_' + ffcx.naming.compute_signature(expressions, '', p)

    expr_names = ["expression_{!s}".format(ffcx.naming.compute_signature([expression], "", p))
                  for expression in expressions]

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        obj, mod = get_cached_module(module_name, expr_names, cache_dir, timeout)
        if obj is not None:
            return obj, mod
    else:
        cache_dir = Path(tempfile.mkdtemp())

    try:
        scalar_type = p["scalar_type"].replace("complex", "_Complex")
        decl = UFC_HEADER_DECL.format(scalar_type) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL + \
            UFC_COORDINATEMAPPING_DECL + UFC_INTEGRAL_DECL + UFC_FORM_DECL + UFC_EXPRESSION_DECL

        expression_template = "ufc_expression* create_{name}(void);\n"
        for name in expr_names:
            decl += expression_template.format(name=name)

        _compile_objects(decl, expressions, expr_names, module_name, p, cache_dir,
                         cffi_extra_compile_args, cffi_verbose, cffi_debug, cffi_libraries)
    except Exception:
        # remove c file so that it will not timeout next time
        c_filename = cache_dir.joinpath(module_name + ".c")
        os.replace(c_filename, c_filename.with_suffix(".c.failed"))
        raise

    obj, module = _load_objects(cache_dir, module_name, expr_names)
    return obj, module


def compile_coordinate_maps(meshes, parameters=None, cache_dir=None, timeout=10, cffi_extra_compile_args=None,
                            cffi_verbose=False, cffi_debug=None, cffi_libraries=None):
    """Compile a list of UFL coordinate mappings into UFC Python objects."""
    p = ffcx.parameters.default_parameters()
    if parameters is not None:
        p.update(parameters)

    p.update(ffcx.parameters.env_parameters())

    # Get a signature for these cmaps
    module_name = 'libffcx_cmaps_' + \
        ffcx.naming.compute_signature(meshes, _compute_parameter_signature(
            p) + str(cffi_extra_compile_args) + str(cffi_debug), True)

    cmap_names = [ffcx.naming.coordinate_map_name(
        mesh.ufl_coordinate_element(), "JIT") for mesh in meshes]

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        obj, mod = get_cached_module(module_name, cmap_names, cache_dir, timeout)
        if obj is not None:
            return obj, mod
    else:
        cache_dir = Path(tempfile.mkdtemp())

    try:
        scalar_type = p["scalar_type"].replace("complex", "_Complex")
        decl = UFC_HEADER_DECL.format(scalar_type) + UFC_COORDINATEMAPPING_DECL + UFC_DOFMAP_DECL
        cmap_template = "ufc_coordinate_mapping * create_{name}(void);\n"

        for name in cmap_names:
            decl += cmap_template.format(name=name)

        _compile_objects(decl, meshes, cmap_names, module_name, p, cache_dir,
                         cffi_extra_compile_args, cffi_verbose, cffi_debug, cffi_libraries)
    except Exception:
        # remove c file so that it will not timeout next time
        c_filename = cache_dir.joinpath(module_name + ".c")
        os.replace(c_filename, c_filename.with_suffix(".c.failed"))
        raise

    obj, module = _load_objects(cache_dir, module_name, cmap_names)
    return obj, module


def _compile_objects(decl, ufl_objects, object_names, module_name, parameters, cache_dir,
                     cffi_extra_compile_args, cffi_verbose, cffi_debug, cffi_libraries):

    import ffcx.compiler

    _, code_body = ffcx.compiler.compile_ufl_objects(ufl_objects, prefix="JIT", parameters=parameters)

    ffibuilder = cffi.FFI()
    ffibuilder.set_source(module_name, code_body, include_dirs=[ffcx.codegeneration.get_include_path()],
                          extra_compile_args=cffi_extra_compile_args, libraries=cffi_libraries)
    ffibuilder.cdef(decl)

    c_filename = cache_dir.joinpath(module_name + ".c")
    ready_name = c_filename.with_suffix(".c.cached")

    # Compile (ensuring that compile dir exists)
    cache_dir.mkdir(exist_ok=True, parents=True)

    logger.info(79 * "#")
    logger.info("Calling JIT C compiler")
    logger.info(79 * "#")

    t0 = time.time()
    f = io.StringIO()
    with redirect_stdout(f):
        ffibuilder.compile(tmpdir=cache_dir, verbose=True, debug=cffi_debug)
    s = f.getvalue()
    if (cffi_verbose):
        print(s)

    logger.info("JIT C compiler finished in {:.4f}".format(time.time() - t0))

    # Create a "status ready" file. If this fails, it is an error,
    # because it should not exist yet.
    # Copy the stdout verbose output of the build into the ready file
    fd = open(ready_name, "x")
    fd.write(s)
    fd.close()


def _load_objects(cache_dir, module_name, object_names):

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

    compiled_objects = []
    for name in object_names:
        # Call UFC factory to create object data struct (calls malloc)
        obj = getattr(compiled_module.lib, "create_" + name)()

        # Set garbage collector to use C free()
        compiled_objects.append(compiled_module.ffi.gc(obj, compiled_module.lib.free))

    return compiled_objects, compiled_module
