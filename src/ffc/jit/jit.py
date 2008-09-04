"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-07-20 -- 2008-09-04"
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
from jitobject import wrap

# Import Instant
import instant

# Global counter for numbering forms
counter = 0

# In-memory form cache
form_cache = {}

# JIT objects (used for caching)
_jit_objects = {}

# Options for JIT-compiler, evaluate_basis and evaluate_basis_derivatives turned off
FFC_OPTIONS_JIT = FFC_OPTIONS.copy()
#FFC_OPTIONS_JIT["no-evaluate_basis"] = True
FFC_OPTIONS_JIT["no-evaluate_basis_derivatives"] = True

def jit(input_form, options=None):
    """Just-in-time compile the given form or element
    
    Parameters:
    
      input_form : The form
      options    : An option dictionary
    """

    # Check options
    _options = check_options(options)
     
    # Set C++ compiler options
    if _options["optimize"]:
        cppargs = "-O2"
    else:
        cppargs = "-O0"

    # Wrap input
    jit_object = wrap(input_form, _options)

    # Check in-memory form cache
    if input_form in form_cache:
        return form_cache[input_form]

    # Check that we don't get a list
    if isinstance(input_form, list):
        raise RuntimeError, "Just-in-time compiler requires a single form (not a list of forms"

    # Analyze and simplify form (to get checksum for simplified form and to get form_data)
    form = algebra.Form(input_form)
    form_data = analyze.analyze(form, simplify_form=False)

    # Compute md5 checksum of form signature
    signature = " ".join([str(form),
                          ", ".join([element.signature() for element in form_data.elements]),
                          _options["representation"], _options["language"], str(_options), cppargs])
    md5sum = "form_" + md5.new(signature).hexdigest()

    # Get name of form
    prefix = md5sum
    rank = form_data.rank
    if rank == 0:
        form_name = prefix + "Functional"
    elif rank == 1:
        form_name = prefix + "LinearForm"
    elif rank == 2:
        form_name = prefix + "BilinearForm"
    else:
        form_name = prefix

    # Check if we can reuse form from cache
    compiled_form = None
    compiled_module = instant.import_module(md5sum)
    if compiled_module:
        debug("Found form in cache, reusing previously built module (checksum %s)" % md5sum[5:], -1)
        try:
            exec("compiled_form = compiled_module.%s()" %form_name)
        except:
            debug("Form module in cache seems to be broken, need to rebuild module", -1)
            compiled_module = False

    # Need to rebuild and import module
    if compiled_module is None:
        
        # Build form module
        compiled_module = build_module(form, _options, md5sum, prefix, cppargs)
        try: 
            exec("compiled_form = compiled_module.%s()" %form_name)
        except:
            debug("Cannot find function %s after loading module, should never happen" %form_name, 1)

    # Add to form cache
    if not input_form in form_cache:
        form_cache[input_form] = (compiled_form, compiled_module, form_data)
    
    return (compiled_form, compiled_module, form_data)

def check_options(options):
    "Check options and add any missing options"
    _options = options.copy()
    for key in options:
        if not key in FFC_OPTIONS:
            warning('Unknown option "%s" for JIT compiler, ignoring.' % key)
    for key in FFC_OPTIONS:
        if not key in options:
            _options[key] = FFC_OPTIONS[key]
    return _options

def build_module(form, options, md5sum, prefix, cppargs):
    "Build module"

    # Compile form
    debug("Calling FFC just-in-time (JIT) compiler, this may take some time...", -1)
    compile(form, prefix, options)
    debug("done", -1)

    # Header file
    filename = prefix + ".h"

    # Get include directory for ufc.h (might be better way to do this?)
    (path, dummy, dummy, dummy) = instant.header_and_libs_from_pkgconfig("ufc-1")
    if len(path) == 0:
        path = [("/").join(sysconfig.get_python_inc().split("/")[:-2]) + "/include"]
    ufc_include = '%%include "%s/ufc.h"' % path[0]

    # Wrap code into a Python module using Instant
    debug("Creating Python extension (compiling and linking), this may take some time...", -1)
    module = instant.build_module(wrap_headers=[filename], additional_declarations=ufc_include, include_dirs=path, cppargs=cppargs, signature=md5sum)
    debug("done", -1)
    return module
