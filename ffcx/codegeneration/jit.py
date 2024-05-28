# Copyright (C) 2004-2019 Garth N. Wells
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Just-in-time compilation."""

from __future__ import annotations

import importlib
import io
import logging
import os
import re
import sys
import sysconfig
import tempfile
import time
from contextlib import redirect_stdout
from pathlib import Path

import cffi
import numpy as np
import numpy.typing as npt
import ufl

import ffcx
import ffcx.naming
from ffcx.codegeneration.C.file_template import libraries as _libraries

logger = logging.getLogger("ffcx")
root_logger = logging.getLogger()

# Get declarations directly from ufcx.h
file_dir = os.path.dirname(os.path.abspath(__file__))
with open(file_dir + "/ufcx.h") as f:
    ufcx_h = "".join(f.readlines())

# Emulate C preprocessor on __STDC_NO_COMPLEX__
if sys.platform.startswith("win32"):
    # Remove macro statements and content
    ufcx_h = re.sub(
        r"\#ifndef __STDC_NO_COMPLEX__.*?\#endif // __STDC_NO_COMPLEX__",
        "",
        ufcx_h,
        flags=re.DOTALL,
    )
else:
    # Remove only macros keeping content
    ufcx_h = ufcx_h.replace("#ifndef __STDC_NO_COMPLEX__", "")
    ufcx_h = ufcx_h.replace("#endif // __STDC_NO_COMPLEX__", "")

header = ufcx_h.split("<HEADER_DECL>")[1].split("</HEADER_DECL>")[0].strip(" /\n")
header = header.replace("{", "{{").replace("}", "}}")
UFC_HEADER_DECL = header + "\n"

UFC_FORM_DECL = "\n".join(re.findall("typedef struct ufcx_form.*?ufcx_form;", ufcx_h, re.DOTALL))

UFC_INTEGRAL_DECL = "\n".join(
    re.findall(r"typedef void ?\(ufcx_tabulate_tensor_float32\).*?\);", ufcx_h, re.DOTALL)
)
UFC_INTEGRAL_DECL += "\n".join(
    re.findall(r"typedef void ?\(ufcx_tabulate_tensor_float64\).*?\);", ufcx_h, re.DOTALL)
)
UFC_INTEGRAL_DECL += "\n".join(
    re.findall(r"typedef void ?\(ufcx_tabulate_tensor_complex64\).*?\);", ufcx_h, re.DOTALL)
)
UFC_INTEGRAL_DECL += "\n".join(
    re.findall(r"typedef void ?\(ufcx_tabulate_tensor_complex128\).*?\);", ufcx_h, re.DOTALL)
)

UFC_INTEGRAL_DECL += "\n".join(
    re.findall("typedef struct ufcx_integral.*?ufcx_integral;", ufcx_h, re.DOTALL)
)

UFC_EXPRESSION_DECL = "\n".join(
    re.findall("typedef struct ufcx_expression.*?ufcx_expression;", ufcx_h, re.DOTALL)
)


def _compute_option_signature(options):
    """Return options signature (some options should not affect signature)."""
    return str(sorted(options.items()))


def get_cached_module(module_name, object_names, cache_dir, timeout):
    """Look for an existing C file and wait for compilation, or if it does not exist, create it."""
    cache_dir = Path(cache_dir)
    c_filename = cache_dir.joinpath(module_name).with_suffix(".c")
    ready_name = c_filename.with_suffix(".c.cached")

    # Ensure cache dir exists
    cache_dir.mkdir(exist_ok=True, parents=True)

    try:
        # Create C file with exclusive access
        with open(c_filename, "x"):
            pass
        return None, None
    except FileExistsError:
        logger.info("Cached C file already exists: " + str(c_filename))
        finder = importlib.machinery.FileFinder(
            str(cache_dir),
            (importlib.machinery.ExtensionFileLoader, importlib.machinery.EXTENSION_SUFFIXES),
        )
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
        raise TimeoutError(
            "JIT compilation timed out, probably due to a failed previous compile. "
            f"Try cleaning cache (e.g. remove {c_filename}) or increase timeout option."
        )


def _compilation_signature(cffi_extra_compile_args, cffi_debug):
    """Compute the compilation-inputs part of the signature.

    Used to avoid cache conflicts across Python versions, architectures, installs.

    - SOABI includes platform, Python version, debug flags
    - CFLAGS includes prefixes, arch targets
    """
    if sys.platform.startswith("win32"):
        # NOTE: SOABI not defined on win32, EXT_SUFFIX contains e.g. '.cp312-win_amd64.pyd'
        return (
            str(cffi_extra_compile_args)
            + str(cffi_debug)
            + str(sysconfig.get_config_var("EXT_SUFFIX"))
        )
    else:
        return (
            str(cffi_extra_compile_args)
            + str(cffi_debug)
            + str(sysconfig.get_config_var("CFLAGS"))
            + str(sysconfig.get_config_var("SOABI"))
        )


def compile_forms(
    forms: list[ufl.Form],
    options: dict = {},
    cache_dir: Path | None = None,
    timeout: int = 10,
    cffi_extra_compile_args: list[str] = [],
    cffi_verbose: bool = False,
    cffi_debug: bool = False,
    cffi_libraries: list[str] = [],
    visualise: bool = False,
):
    """Compile a list of UFL forms into UFC Python objects.

    Args:
        forms: List of ufl.form to compile.
        options: Options
        cache_dir: Cache directory
        timeout: Timeout
        cffi_extra_compile_args: Extra compilation args for CFFI
        cffi_verbose: Use verbose compile
        cffi_debug: Use compiler debug mode
        cffi_libraries: libraries to use with compiler
        visualise: Toggle visualisation
    """
    p = ffcx.options.get_options(options)

    # Get a signature for these forms
    module_name = "libffcx_forms_" + ffcx.naming.compute_signature(
        forms,
        _compute_option_signature(p) + _compilation_signature(cffi_extra_compile_args, cffi_debug),
    )

    form_names = [ffcx.naming.form_name(form, i, module_name) for i, form in enumerate(forms)]

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        obj, mod = get_cached_module(module_name, form_names, cache_dir, timeout)
        if obj is not None:
            return obj, mod, (None, None)
    else:
        cache_dir = Path(tempfile.mkdtemp())

    try:
        decl = (
            UFC_HEADER_DECL.format(np.dtype(p["scalar_type"]).name)  # type: ignore
            + UFC_INTEGRAL_DECL
            + UFC_FORM_DECL
        )

        form_template = "extern ufcx_form {name};\n"
        for name in form_names:
            decl += form_template.format(name=name)

        impl = _compile_objects(
            decl,
            forms,
            form_names,
            module_name,
            p,
            cache_dir,
            cffi_extra_compile_args,
            cffi_verbose,
            cffi_debug,
            cffi_libraries,
            visualise=visualise,
        )
    except Exception as e:
        try:
            # remove c file so that it will not timeout next time
            c_filename = cache_dir.joinpath(module_name + ".c")
            os.replace(c_filename, c_filename.with_suffix(".c.failed"))
        except Exception:
            pass
        raise e

    obj, module = _load_objects(cache_dir, module_name, form_names)
    return obj, module, (decl, impl)


def compile_expressions(
    expressions: list[tuple[ufl.Expr, npt.NDArray[np.floating]]],
    options: dict = {},
    cache_dir: Path | None = None,
    timeout: int = 10,
    cffi_extra_compile_args: list[str] = [],
    cffi_verbose: bool = False,
    cffi_debug: bool = False,
    cffi_libraries: list[str] = [],
    visualise: bool = False,
):
    """Compile a list of UFL expressions into UFC Python objects.

    Args:
        expressions: List of (UFL expression, evaluation points).
        options: Options
        cache_dir: Cache directory
        timeout: Timeout
        cffi_extra_compile_args: Extra compilation args for CFFI
        cffi_verbose: Use verbose compile
        cffi_debug: Use compiler debug mode
        cffi_libraries: libraries to use with compiler
        visualise: Toggle visualisation
    """
    p = ffcx.options.get_options(options)

    module_name = "libffcx_expressions_" + ffcx.naming.compute_signature(
        expressions,
        _compute_option_signature(p) + _compilation_signature(cffi_extra_compile_args, cffi_debug),
    )
    expr_names = [
        ffcx.naming.expression_name(expression, module_name) for expression in expressions
    ]

    if cache_dir is not None:
        cache_dir = Path(cache_dir)
        obj, mod = get_cached_module(module_name, expr_names, cache_dir, timeout)
        if obj is not None:
            return obj, mod, (None, None)
    else:
        cache_dir = Path(tempfile.mkdtemp())

    try:
        decl = (
            UFC_HEADER_DECL.format(np.dtype(p["scalar_type"]).name)  # type: ignore
            + UFC_INTEGRAL_DECL
            + UFC_FORM_DECL
            + UFC_EXPRESSION_DECL
        )

        expression_template = "extern ufcx_expression {name};\n"
        for name in expr_names:
            decl += expression_template.format(name=name)

        impl = _compile_objects(
            decl,
            expressions,
            expr_names,
            module_name,
            p,
            cache_dir,
            cffi_extra_compile_args,
            cffi_verbose,
            cffi_debug,
            cffi_libraries,
            visualise=visualise,
        )
    except Exception as e:
        try:
            # remove c file so that it will not timeout next time
            c_filename = cache_dir.joinpath(module_name + ".c")
            os.replace(c_filename, c_filename.with_suffix(".c.failed"))
        except Exception:
            pass
        raise e

    obj, module = _load_objects(cache_dir, module_name, expr_names)
    return obj, module, (decl, impl)


def _compile_objects(
    decl,
    ufl_objects,
    object_names,
    module_name,
    options,
    cache_dir,
    cffi_extra_compile_args,
    cffi_verbose,
    cffi_debug,
    cffi_libraries,
    visualise: bool = False,
):
    import ffcx.compiler

    libraries = _libraries + cffi_libraries if cffi_libraries is not None else _libraries

    # JIT uses module_name as prefix, which is needed to make names of all struct/function
    # unique across modules
    _, code_body, _ = ffcx.compiler.compile_ufl_objects(
        ufl_objects, prefix=module_name, options=options, visualise=visualise
    )

    # Raise error immediately prior to compilation if no support for C99
    # _Complex. Doing this here allows FFCx to be used for complex codegen on
    # Windows.
    if sys.platform.startswith("win32"):
        if np.issubdtype(options["scalar_type"], np.complexfloating):
            raise NotImplementedError("win32 platform does not support C99 _Complex numbers")
        elif isinstance(options["scalar_type"], str) and "complex" in options["scalar_type"]:
            raise NotImplementedError("win32 platform does not support C99 _Complex numbers")

    # Compile in C17 mode
    if sys.platform.startswith("win32"):
        cffi_base_compile_args = ["-std:c17"]
    else:
        cffi_base_compile_args = ["-std=c17"]

    cffi_final_compile_args = cffi_base_compile_args + cffi_extra_compile_args

    ffibuilder = cffi.FFI()

    ffibuilder.set_source(
        module_name,
        code_body,
        include_dirs=[ffcx.codegeneration.get_include_path()],
        extra_compile_args=cffi_final_compile_args,
        libraries=libraries,
    )

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
    # Temporarily set root logger handlers to string buffer only
    # since CFFI logs into root logger
    old_handlers = root_logger.handlers.copy()
    root_logger.handlers = [logging.StreamHandler(f)]
    with redirect_stdout(f):
        ffibuilder.compile(tmpdir=cache_dir, verbose=True, debug=cffi_debug)
    s = f.getvalue()
    if cffi_verbose:
        print(s)

    logger.info(f"JIT C compiler finished in {time.time() - t0:.4f}")

    # Create a "status ready" file. If this fails, it is an error,
    # because it should not exist yet.
    # Copy the stdout verbose output of the build into the ready file
    fd = open(ready_name, "x")
    fd.write(s)
    fd.close()

    # Copy back the original handlers (in case someone is logging into
    # root logger and has custom handlers)
    root_logger.handlers = old_handlers

    return code_body


def _load_objects(cache_dir, module_name, object_names):
    # Create module finder that searches the compile path
    finder = importlib.machinery.FileFinder(
        str(cache_dir),
        (importlib.machinery.ExtensionFileLoader, importlib.machinery.EXTENSION_SUFFIXES),
    )

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
