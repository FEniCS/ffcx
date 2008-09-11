"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-07-20 -- 2008-09-06"
__copyright__ = "Copyright (C) 2007-2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Johan Hake, 2008
# Mofified by Ilmar Wilbers, 2008

# Python modules
from os import system
from commands import getoutput
from distutils import sysconfig
import md5, os, sys, shutil

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *

# FFC compiler modules
from ffc.compiler.compiler import compile
from ffc.compiler.language import algebra
from ffc.compiler.analysis import simplify, analyze

# FFC jit modules
from jitobject import JITObject

# Import Instant
import instant

# In-memory form cache
form_cache = {}

# FIXME: Change here for testing
use_ffc_cache = True

# Options for JIT-compiler, evaluate_basis and evaluate_basis_derivatives turned off
FFC_OPTIONS_JIT = FFC_OPTIONS.copy()
FFC_OPTIONS_JIT["no-evaluate_basis_derivatives"] = True

def jit(form, options=None):
    """Just-in-time compile the given form or element
    
    Parameters:
    
      form    : The form
      options : An option dictionary
    """

    print ""
    print "--- Calling FFC JIT compiler ---"

    # Check options
    options = check_options(form, options)

    # Wrap input
    jit_object = JITObject(form, options)

    # Check in-memory form cache
    if use_ffc_cache and form in form_cache: return form_cache[form]
    
    # Check cache
    module = instant.import_module(jit_object, cache_dir=options["cache_dir"])
    if module: return extract_form(form, module, jit_object)

    # Generate code
    debug("Calling FFC just-in-time (JIT) compiler, this may take some time...", -1)
    signature = jit_object.signature()
    compile(form, signature, options)
    debug("done", -1)

    # Wrap code into a Python module using Instant
    debug("Creating Python extension (compiling and linking), this may take some time...", -1)
    filename = signature + ".h"
    (cppargs, path, ufc_include) = extract_instant_flags(options)
    module = instant.build_module(wrap_headers=[filename],
                                  additional_declarations=ufc_include,
                                  include_dirs=path,
                                  cppargs=cppargs,
                                  signature=signature,
                                  cache_dir=options["cache_dir"])
    debug("done", -1)

    return extract_form(form, module, jit_object)

def check_options(form, options):
    "Check options and add any missing options"

    # Form can not be a list
    if isinstance(form, list):
        raise RuntimeError, "JIT compiler requires a single form (not a list of forms)."

    # Copy options
    options = options.copy()

    # Check for invalid options
    for key in options:
        if not key in FFC_OPTIONS_JIT:
            warning('Unknown option "%s" for JIT compiler, ignoring.' % key)

    # Add defaults for missing options
    for key in FFC_OPTIONS_JIT:
        if not key in options:
            options[key] = FFC_OPTIONS_JIT[key]

    # Don't postfix form names
    if "form_postfix" in options and options["form_postfix"]:
        warning("Forms cannot be postfixed when the JIT compiler is used.")
    options["form_postfix"] = False

    return options

def extract_form(form, module, jit_object):
    "Extract form from module"
    signature = jit_object.signature()
    compiled_form = getattr(module, signature)()
    form_data = jit_object.form_data
    if use_ffc_cache: form_cache[form] = (compiled_form, module, form_data)
    return (compiled_form, module, form_data)

def extract_instant_flags(options):
    "Extract flags for Instant"

    # Get C++ compiler options
    if options["optimize"]: cppargs = "-O2"
    else: cppargs = "-O0"

    # Get include directory for ufc.h (might be better way to do this?)
    (path, dummy, dummy, dummy) = instant.header_and_libs_from_pkgconfig("ufc-1")
    if len(path) == 0: path = [("/").join(sysconfig.get_python_inc().split("/")[:-2]) + "/include"]
    ufc_include = '%%include "%s/ufc.h"' % path[0]

    return (cppargs, path, ufc_include)
