# Copyright (C) 2004-2019 Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import importlib
import io
import logging
import os
import re
import tempfile
import time
from contextlib import redirect_stdout
from pathlib import Path

import cffi
import ffcx
import ffcx.naming

logger = logging.getLogger("ffcx")

# Get declarations directly from ufcx.h
file_dir = os.path.dirname(os.path.abspath(__file__))
with open(file_dir + "/ufcx.h", "r") as f:
    ufcx_h = ''.join(f.readlines())

header = ufcx_h.split("<HEADER_DECL>")[1].split("</HEADER_DECL>")[0].strip(" /\n")
header = header.replace("{", "{{").replace("}", "}}")
UFC_HEADER_DECL = header + "\n"

UFC_ELEMENT_DECL = '\n'.join(re.findall('typedef struct ufcx_finite_element.*?ufcx_finite_element;', ufcx_h, re.DOTALL))
UFC_DOFMAP_DECL = '\n'.join(re.findall('typedef struct ufcx_dofmap.*?ufcx_dofmap;', ufcx_h, re.DOTALL))
UFC_FORM_DECL = '\n'.join(re.findall('typedef struct ufcx_form.*?ufcx_form;', ufcx_h, re.DOTALL))

UFC_INTEGRAL_DECL = '\n'.join(re.findall(r'typedef void ?\(ufcx_tabulate_tensor_float32\).*?\);', ufcx_h, re.DOTALL))
UFC_INTEGRAL_DECL += '\n'.join(re.findall(r'typedef void ?\(ufcx_tabulate_tensor_float64\).*?\);', ufcx_h, re.DOTALL))
UFC_INTEGRAL_DECL += '\n'.join(re.findall(r'typedef void ?\(ufcx_tabulate_tensor_complex64\).*?\);', ufcx_h, re.DOTALL))
UFC_INTEGRAL_DECL += '\n'.join(re.findall(r'typedef void ?\(ufcx_tabulate_tensor_complex128\).*?\);',
                               ufcx_h, re.DOTALL))
UFC_INTEGRAL_DECL += '\n'.join(re.findall(r'typedef void ?\(ufcx_tabulate_tensor_longdouble\).*?\);',
                               ufcx_h, re.DOTALL))

UFC_INTEGRAL_DECL += '\n'.join(re.findall('typedef struct ufcx_integral.*?ufcx_integral;',
                                          ufcx_h, re.DOTALL))
UFC_EXPRESSION_DECL = '\n'.join(re.findall('typedef struct ufcx_expression.*?ufcx_expression;', ufcx_h, re.DOTALL))


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

                compiled_objects = [getattr(compiled_module.lib, name) for name in object_names]
                return compiled_objects, compiled_module

            logger.info(f"Waiting for {ready_name} to appear.")
            time.sleep(1)
        raise TimeoutError(f"""JIT compilation timed out, probably due to a failed previous compile.
        Try cleaning cache (e.g. remove {c_filename}) or increase timeout parameter.""")


def compile_elements(elements, parameters=None, cache_dir=None, timeout=10, cffi_extra_compile_args=None,
                     cffi_verbose=False, cffi_debug=None, cffi_libraries=None):
    """Compile a list of UFL elements and dofmaps into Python objects."""
    p = ffcx.parameters.get_parameters(parameters)

    # Get a signature for these elements
    module_name = 'libffcx_elements_' + \
        ffcx.naming.compute_signature(elements, _compute_parameter_signature(p)
                                      + str(cffi_extra_compile_args) + str(cffi_debug))

    names = []
    for e in elements:
        name = ffcx.naming.finite_element_name(e, module_name)
        names.append(name)
        name = ffcx.naming.dofmap_name(e, module_name)
        names.append(name)

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        obj, mod = get_cached_module(module_name, names, cache_dir, timeout)
        if obj is not None:
            # Pair up elements with dofmaps
            obj = list(zip(obj[::2], obj[1::2]))
            return obj, mod, (None, None)
    else:
        cache_dir = Path(tempfile.mkdtemp())

    try:
        decl = UFC_HEADER_DECL.format(p["scalar_type"]) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL
        element_template = "extern ufcx_finite_element {name};\n"
        dofmap_template = "extern ufcx_dofmap {name};\n"
        for i in range(len(elements)):
            decl += element_template.format(name=names[i * 2])
            decl += dofmap_template.format(name=names[i * 2 + 1])

        impl = _compile_objects(decl, elements, names, module_name, p, cache_dir,
                                cffi_extra_compile_args, cffi_verbose, cffi_debug, cffi_libraries)
    except Exception:
        # remove c file so that it will not timeout next time
        c_filename = cache_dir.joinpath(module_name + ".c")
        os.replace(c_filename, c_filename.with_suffix(".c.failed"))
        raise

    objects, module = _load_objects(cache_dir, module_name, names)
    # Pair up elements with dofmaps
    objects = list(zip(objects[::2], objects[1::2]))
    return objects, module, (decl, impl)


def compile_forms(forms, parameters=None, cache_dir=None, timeout=10, cffi_extra_compile_args=None,
                  cffi_verbose=False, cffi_debug=None, cffi_libraries=None):
    """Compile a list of UFL forms into UFC Python objects."""
    p = ffcx.parameters.get_parameters(parameters)

    # Get a signature for these forms
    module_name = 'libffcx_forms_' + \
        ffcx.naming.compute_signature(forms, _compute_parameter_signature(p)
                                      + str(cffi_extra_compile_args) + str(cffi_debug))

    form_names = [ffcx.naming.form_name(form, i, module_name) for i, form in enumerate(forms)]

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        obj, mod = get_cached_module(module_name, form_names, cache_dir, timeout)
        if obj is not None:
            return obj, mod, (None, None)
    else:
        cache_dir = Path(tempfile.mkdtemp())

    try:
        decl = UFC_HEADER_DECL.format(p["scalar_type"]) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL + \
            UFC_INTEGRAL_DECL + UFC_FORM_DECL

        form_template = "extern ufcx_form {name};\n"
        for name in form_names:
            decl += form_template.format(name=name)

        impl = _compile_objects(decl, forms, form_names, module_name, p, cache_dir,
                                cffi_extra_compile_args, cffi_verbose, cffi_debug, cffi_libraries)
    except Exception:
        # remove c file so that it will not timeout next time
        c_filename = cache_dir.joinpath(module_name + ".c")
        os.replace(c_filename, c_filename.with_suffix(".c.failed"))
        raise

    obj, module = _load_objects(cache_dir, module_name, form_names)
    return obj, module, (decl, impl)


def compile_expressions(expressions, parameters=None, cache_dir=None, timeout=10, cffi_extra_compile_args=None,
                        cffi_verbose=False, cffi_debug=None, cffi_libraries=None):
    """Compile a list of UFL expressions into UFC Python objects.

    Parameters
    ----------
    expressions
        List of (UFL expression, evaluation points).

    """
    p = ffcx.parameters.get_parameters(parameters)

    module_name = 'libffcx_expressions_' + \
        ffcx.naming.compute_signature(expressions, _compute_parameter_signature(p)
                                      + str(cffi_extra_compile_args) + str(cffi_debug))
    expr_names = [ffcx.naming.expression_name(expression, module_name) for expression in expressions]

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        obj, mod = get_cached_module(module_name, expr_names, cache_dir, timeout)
        if obj is not None:
            return obj, mod, (None, None)
    else:
        cache_dir = Path(tempfile.mkdtemp())

    try:
        decl = UFC_HEADER_DECL.format(p["scalar_type"]) + UFC_ELEMENT_DECL + UFC_DOFMAP_DECL + \
            UFC_INTEGRAL_DECL + UFC_FORM_DECL + UFC_EXPRESSION_DECL

        expression_template = "extern ufcx_expression {name};\n"
        for name in expr_names:
            decl += expression_template.format(name=name)

        impl = _compile_objects(decl, expressions, expr_names, module_name, p, cache_dir,
                                cffi_extra_compile_args, cffi_verbose, cffi_debug, cffi_libraries)
    except Exception:
        # remove c file so that it will not timeout next time
        c_filename = cache_dir.joinpath(module_name + ".c")
        os.replace(c_filename, c_filename.with_suffix(".c.failed"))
        raise

    obj, module = _load_objects(cache_dir, module_name, expr_names)
    return obj, module, (decl, impl)


def _compile_objects(decl, ufl_objects, object_names, module_name, parameters, cache_dir,
                     cffi_extra_compile_args, cffi_verbose, cffi_debug, cffi_libraries):

    import ffcx.compiler

    # JIT uses module_name as prefix, which is needed to make names of all struct/function
    # unique across modules
    _, code_body = ffcx.compiler.compile_ufl_objects(ufl_objects, prefix=module_name, parameters=parameters)

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

    return code_body


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
        obj = getattr(compiled_module.lib, name)
        compiled_objects.append(obj)

    return compiled_objects, compiled_module
