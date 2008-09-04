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
    options = check_options(input_form, options)

    # Testing new design
    #return new_jit(input_form, options)
    return old_jit(input_form, options)

def check_options(form, options):
    "Check options and add any missing options"

    # Form can not be a list
    if isinstance(form, list):
        raise RuntimeError, "JIT compiler requires a single form (not a list of forms)."

    # Copy options
    options = options.copy()

    # Check for invalid options
    for key in options:
        if not key in FFC_OPTIONS:
            warning('Unknown option "%s" for JIT compiler, ignoring.' % key)

    # Add defaults for missing options
    for key in FFC_OPTIONS:
        if not key in options:
            options[key] = FFC_OPTIONS[key]

    # Don't postfix form names
    if "form_postfix" in options and options["form_postfix"]:
        warning("Forms cannot be postfixed when the JIT compiler is used.")
    options["form_postfix"] = False

    return options

def extract_instant_flags(options):
    "Extract flags for Instant"

    # Get C++ compiler options
    if options["optimize"]:
        cppargs = "-O2"
    else:
        cppargs = "-O0"

    # Get include directory for ufc.h (might be better way to do this?)
    (path, dummy, dummy, dummy) = instant.header_and_libs_from_pkgconfig("ufc-1")
    if len(path) == 0: path = [("/").join(sysconfig.get_python_inc().split("/")[:-2]) + "/include"]
    ufc_include = '%%include "%s/ufc.h"' % path[0]

    return (cppargs, path, ufc_include)

def old_jit(input_form, options):

    # Analyze and simplify form (to get checksum for simplified form and to get form_data)
    form = algebra.Form(input_form)
    form_data = analyze.analyze(form, simplify_form=False)

    (cppargs, path, ufc_include) = extract_instant_flags(options)

    # Compute md5 checksum of form signature
    signature = " ".join([str(form),
                          ", ".join([element.signature() for element in form_data.elements]),
                          options["representation"], options["language"], str(options), cppargs])
    md5sum = "form_" + md5.new(signature).hexdigest()

    # Get name of form
    prefix = md5sum

    # Check if we can reuse form from cache
    compiled_form = None
    compiled_module = instant.import_module(md5sum)

    # Need to rebuild and import module
    if compiled_module is None:
        
        # Build form module
        compile(form, prefix, options)

   
        # Wrap code into a Python module using Instant
        filename = prefix + ".h"
        compiled_module = instant.build_module(wrap_headers=[filename], additional_declarations=ufc_include, include_dirs=path, cppargs=cppargs, signature=md5sum)
        
        try: 
            exec("compiled_form = compiled_module.%s()" % prefix)
        except:
            debug("Cannot find function %s after loading module, should never happen" % prefix, 1)

    # Add to form cache
    if not input_form in form_cache:
        form_cache[input_form] = (compiled_form, compiled_module, form_data)
    
    return (compiled_form, compiled_module, form_data)

def new_jit(input_form, options):

    # Check in-memory form cache
    if input_form in form_cache:
        return form_cache[input_form]

    # Wrap input
    jit_object = wrap(input_form, options)

    # Check cache
    module = instant.import_module(jit_object)
    print module

    # Compile form
    debug("Calling FFC just-in-time (JIT) compiler, this may take some time...", -1)
    compile(input_form, jit_object.signature(), options)
    debug("done", -1)

    # Wrap code into a Python module using Instant
    debug("Creating Python extension (compiling and linking), this may take some time...", -1)
    signature = jit_object.signature()
    filename = signature + ".h"
    (cppargs, path, ufc_include) = extract_instant_flags(options)
    module = instant.build_module(wrap_headers=[filename],
                                  additional_declarations=ufc_include,
                                  include_dirs=path,
                                  cppargs=cppargs,
                                  signature=signature)
    debug("done", -1)

    # Extract form
    exec("compiled_form = module.%s()" % signature)

    # Add to form cache
    if not input_form in form_cache:
        form_cache[input_form] = (compiled_form, module, jit_object.form_data)

    return (compiled_form, module, jit_object.form_data)
