"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-07-20 -- 2008-03-31"
__copyright__ = "Copyright (C) 2007-2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

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

# Global counter for numbering forms
counter = 0

# Options for JIT-compiler, evaluate_basis and evaluate_basis_derivatives turned off
FFC_OPTIONS_JIT = FFC_OPTIONS.copy()
#FFC_OPTIONS_JIT["no-evaluate_basis"] = True
FFC_OPTIONS_JIT["no-evaluate_basis_derivatives"] = True

# Compiler options, don't optimize by default (could be added to options)
CPP_ARGS = "-O0"

def jit(form, representation=FFC_REPRESENTATION, language=FFC_LANGUAGE, options=FFC_OPTIONS_JIT):
    "Just-in-time compile the given form or element"

    # Check that we don't get a list
    if isinstance(form, list):
        raise RuntimeError, "Just-in-time compiler requires a single form (not a list of forms"

    # Analyze and simplify form (to get checksum for simplified form and to get form_data)
    form = algebra.Form(form)
    form_data = analyze.analyze(form, simplify_form=False)

    # Compute md5 checksum of form signature
    signature = " ".join([str(form),
                          ", ".join([element.signature() for element in form_data.elements]),
                          representation, language, str(options)])
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

    # Make sure cache directory exists
    cache_dir = os.path.join((os.environ['HOME']), ".ffc", "cache")
    if not os.path.isdir(cache_dir):
        debug("Creating FFC form cache %s" % cache_dir, -1)
        os.makedirs(cache_dir)

    # Make sure any previous module from current director (to not confuse Instant)
    local_dir = prefix + "_module"
    if os.path.isdir(local_dir):
        shutil.rmtree(local_dir)

    # Check if we can reuse form from cache
    form_dir = os.path.join(cache_dir, md5sum)
    module_dir = os.path.join(form_dir, md5sum + "_module")
    if os.path.isdir(form_dir):
        debug("Found form in cache, reusing previously built module (checksum %s)" % md5sum[5:], 1)
        sys.path.append(form_dir)
        try:
            exec("import %s as compiled_module" % (md5sum + "_module"))
            exec("compiled_form = compiled_module.%s()" % form_name)
            return (compiled_form, compiled_module, form_data)
        except:
            debug("Form module in cache seems to be broken, need to rebuild module", -1)
            shutil.rmtree(form_dir, form_dir + "_broken")

    # Build form module
    build_module(form, representation, language, options, md5sum, form_dir, module_dir, prefix)

    # Return the form, module and form data
    sys.path.append(form_dir)
    exec("import %s as compiled_module" % (md5sum + "_module"))
    exec("compiled_form = compiled_module.%s()" % form_name)
    
    return (compiled_form, compiled_module, form_data)

def build_module(form, representation, language, options, md5sum, form_dir, module_dir, prefix):
    "Build module"

    # Make sure form directory exists
    os.makedirs(form_dir)

    # Compile form
    debug("Calling FFC just-in-time (JIT) compiler, this may take some time...", -1)
    compile(form, prefix, representation, language, options)
    debug("done", -1)

    # Move code to cache
    filename = os.path.join(form_dir, prefix + ".h")
    shutil.move(prefix + ".h", filename)

    # FIXME: Move this to top when we have added dependence on Instant
    import instant
    instant.USE_CACHE = 0
    instant.VERBOSE = 0
    instant.COPY_LOCAL_FILES = 0

    # Get include directory for ufc.h (might be better way to do this?)
    (path, dummy, dummy, dummy) = instant.header_and_libs_from_pkgconfig("ufc-1")
    if len(path) == 0:
        path = [("/").join(sysconfig.get_python_inc().split("/")[:-2]) + "/include"]
    ufc_include = '%%include "%s/ufc.h"' % path[0]

    # Wrap code into a Python module using Instant
    debug("Creating Python extension (compiling and linking), this may take some time...", -1)
    module_name = prefix + "_module"
    instant.create_extension(wrap_headers=[filename], module=module_name, additional_declarations=ufc_include, include_dirs=path, cppargs=CPP_ARGS)
    debug("done", -1)

    # Move module to cache
    shutil.move(module_name, module_dir)
