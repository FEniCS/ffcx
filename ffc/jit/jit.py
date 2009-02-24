"""This module provides a just-in-time (JIT) form compiler.
It uses Instant to wrap the generated code into a Python module."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-07-20 -- 2009-01-12"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Johan Hake, 2008 - 2009
# Modified by Ilmar Wilbers, 2008

# Python modules
import instant
from distutils import sysconfig
import os
import re
# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *

# FFC fem modules
from ffc.fem.finiteelement import FiniteElement
from ffc.fem.mixedelement import MixedElement

# FFC compiler modules
from ffc.compiler.compiler import compile

try:
    from ffc.compiler.uflcompiler import compile as uflcompile
    from ufl.classes import Form as UFLForm
    from ufl.classes import FiniteElementBase
except:
    pass

from ffc.compiler.language.algebra import Form, TestFunction
from ffc.compiler.language.builtins import dx

# FFC jit modules
from jitobject import JITObject

# Special Options for JIT-compilation
FFC_OPTIONS_JIT = FFC_OPTIONS.copy()
FFC_OPTIONS_JIT["no-evaluate_basis_derivatives"] = True

# Set debug level for Instant
instant.set_logging_level("warning")

def jit(object, options=None):
    """Just-in-time compile the given form or element
    
    Parameters:
    
      object  : The object to be compiled
      options : An option dictionary
    """

    # Check if we get an element or a form
    if isinstance(object, FiniteElement) or isinstance(object, MixedElement):
        return jit_element(object, options)
    else:
        return jit_form(object, options)


def jit_form(form, options=None):
    "Just-in-time compile the given form"

    # FIXME: Remove this test later
    if not options is None and "compiler" in options and options["compiler"] == "ufl":
        if not isinstance(form, UFLForm):
            form = UFLForm(form)
    else:
        if not isinstance(form, Form):
            form = Form(form)

    # Check options
    options = check_options(form, options)
    
    # Wrap input
    jit_object = JITObject(form, options)
    
    # Check cache
    signature = jit_object.signature()
    module = instant.import_module(jit_object, cache_dir=options["cache_dir"])
    if module: return extract_form(jit_object, module)
    
    # Generate code
    debug("Calling FFC just-in-time (JIT) compiler, this may take some time...", -1)
#    signature = jit_object.signature()

    if options["compiler"] == "ffc":
        compile(form, signature, options)
    elif options["compiler"] == "ufl":
        uflcompile(form, signature, options)
    else:
        raise RuntimeError(options["compiler"], "Unknown compiler.")

    debug("done", -1)
    
    # Wrap code into a Python module using Instant
    debug("Creating Python extension (compiling and linking), this may take some time...", -1)
    filename = signature + ".h"
    
    # Check system requirements
    (cppargs, cpp_path, swig_path) = configure_instant(options)

    # Add shared_ptr declarations
    declarations = extract_declarations(filename)
    
    # Set system headers
    system_headers = ["boost/shared_ptr.hpp"]
    
    # Call instant
    module = instant.build_module(wrap_headers=[filename],
                                  additional_declarations = declarations,
                                  system_headers=system_headers,
                                  include_dirs=cpp_path,
                                  swigargs=['-c++','-I.'],
                                  swig_include_dirs=swig_path,
                                  cppargs=cppargs,
                                  signature=signature,
                                  cache_dir=options["cache_dir"])
    os.unlink(filename)
    debug("done", -1)

    return extract_form(jit_object, module)

def jit_element(element, options=None):
    "Just-in-time compile the given element"
    
    # Check that we get an element
    if not (isinstance(element, FiniteElement) or isinstance(element, MixedElement)):
        raise RuntimeError, "Expecting a finite element."

    # Create simplest possible dummy form
    if element.value_dimension(0) > 1:
        form = TestFunction(element)[0]*dx
    else:
        form = TestFunction(element)*dx

    # Compile form
    (compiled_form, module, form_data) = jit_form(form, options)

    return extract_element_and_dofmap(module)

def check_options(form, options):
    "Check options and add any missing options"

    # Form can not be a list
    if isinstance(form, list):
        raise RuntimeError, "JIT compiler requires a single form (not a list of forms)."

    # Copy options
    if options is None:
        options = {}
    else:
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
    options["form_postfix"] = False

    return options

def extract_form(jit_object, module):
    "Extract form from module"
    return (getattr(module, module.__name__)(), module, jit_object.form_data)

def extract_element_and_dofmap(module):
    "Extract element and dofmap from module"
    name = module.__name__
    return (getattr(module, name + "_finite_element_0")(),
            getattr(module, name + "_dof_map_0")())

def configure_instant(options):
    "Check system requirements"

    # Get C++ compiler options
    if options["cpp optimize"]: cppargs = "-O2"
    else: cppargs = "-O0"

    # Get include directory for ufc.h (might be better way to do this?)
    (path, dummy, dummy, dummy) = instant.header_and_libs_from_pkgconfig("ufc-1")
    if len(path) == 0: path = [(os.sep).join(sysconfig.get_python_inc().split(os.sep)[:-2]) + os.sep + "include"]
    
    # Register the paths
    cpp_path, swig_path = [path[0]], [path[0]]
    
    # Check for swig installation
    result, output = instant.get_status_output("swig -version")
    if result == 1:
        raise OSError, "Could not find swig installation. Please install swig version 1.3.35 or higher.\n"

    # Check swig version for shared_ptr
    if not instant.check_swig_version("1.3.35"):
        raise OSError, "Your current swig version is %s, it needs to be 1.3.35 or higher.\n"%swig_version

    # Check if UFC is importable and what version of swig was used to
    # create the UFC extension module
    try: import ufc
    except: raise OSError, "Please install the python extenstion module of UFC on your system.\n"
    
    # Check that the form compiler will use the same swig version 
    # that UFC was compiled with
    if not instant.check_swig_version(ufc.__swigversion__,same=True):
        raise OSError, """The python extension module of UFC was not compiled with the present version of swig.
Install swig version %s or recompiled UFC with present swig
"""%ufc.__swigversion__
    
    # Set a default swig command and boost include directory
    boost_include_dir = []
    
    # Check for boost installation
    # Set a default directory for the boost installation
    if sys.platform == "darwin":
        # use fink as default
        default = os.path.join(os.path.sep,"sw")
    else:
        default = os.path.join(os.path.sep,"usr")
        
    # If BOOST_DIR is not set use default directory
    boost_dir = os.getenv("BOOST_DIR", default)
    boost_is_found = False
    for inc_dir in ["", "include"]:
        if os.path.isfile(os.path.join(boost_dir, inc_dir, "boost", "version.hpp")):
            boost_include_dir = [os.path.join(boost_dir, inc_dir)]
            boost_is_found = True
            break
        
    if not boost_is_found:
        raise OSError, """The Boost library was not found.
If Boost is installed in a nonstandard location,
set the environment variable BOOST_DIR.
"""

    # Add the boost_include_dir
    cpp_path += boost_include_dir
    return cppargs, cpp_path, swig_path

def extract_declarations(filename):
    "Extract information for shared_ptr"

    # Read the code
    code = open(filename).read()
    
    # Extract the class names
    derived_classes   = re.findall(r"class[ ]+([\w]+)[ ]*: public",code)
    ufc_classes       = re.findall(r"public[ ]+(ufc::[\w]+).*",code)
    ufc_proxy_classes = [s.replace("ufc::","") for s in ufc_classes]
    
    # Swig declarations
    declarations =r"""
%pythoncode %{
import ufc
'''
A hack to get passed a bug in swig.
This is fixed in swig version 1.3.37
%}
%import "swig/ufc.i"
%pythoncode %{
'''
%}
//Uncomment these to produce code for std::tr1::shared_ptr
//#define SWIG_SHARED_PTR_NAMESPACE std
//#define SWIG_SHARED_PTR_SUBNAMESPACE tr1
%include <boost_shared_ptr.i>
"""

    shared_ptr_format = "SWIG_SHARED_PTR_DERIVED(%(der_class)s,%(ufc_class)s,%(der_class)s)"
    declarations += "\n".join(\
        shared_ptr_format%{"ufc_proxy_class":c[0],"ufc_class":c[1],"der_class":c[2]}\
        for c in zip(ufc_proxy_classes, ufc_classes, derived_classes)\
        )

    # Register new objects created by a .create_foo function
    # FIXME: Add any created sub_elements and sub_dofmaps too.
    new_objects = []
    for ufc_class in ufc_classes:
        stripped = ufc_class.replace("ufc::","")
        if stripped not in new_objects:
            new_objects.append(stripped)
    
    new_object_str = "%%newobject %s::create_"%filename.replace(".h","")
    declarations += "\n".join(new_object_str + new_object + ";" \
                                         for new_object in new_objects)
    
    return declarations

