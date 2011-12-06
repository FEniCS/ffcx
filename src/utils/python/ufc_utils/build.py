__author__ = "Johan Hake (hake@simula.no)"
__date__ = "2009-03-06 -- 2011-12-06"
__copyright__ = "Copyright (C) 2009 Johan Hake"
__license__  = "GNU LGPL Version 2.1"

__all__ = ['build_ufc_module']

# Modified by Martin Alnes, 2009

import instant
import os, sys, re, glob

from distutils import sysconfig

def build_ufc_module(h_files, source_directory="", system_headers=None, \
                     swig_binary="swig", swig_path="", **kwargs):
    """Build a python extension module from ufc compliant source code.

    The compiled module will be imported and returned by the function.

    @param h_files:
       The name(s) of the header files that should be compiled and included in
       the python extension module.
    @param source_directory:
       The directory where the source files reside.
    @param system_headers:
       Extra headers that will be #included in the generated wrapper file.
    @param swig_binary:
       Name of the swig binary instant need to look for
    @param swig_path:
       Path to the swig binary

    Any additional keyword arguments are passed on to instant.build_module.
    """

    # Check h_files argument
    if isinstance(h_files, str):
        h_files = [h_files]
    assert isinstance(h_files, list) , "Provide a 'list' or a 'str' as 'h_files'."
    assert all(isinstance(f, str) for f in h_files), \
           "Elements of 'h_files' must be 'str'."

    h_files2 = [os.path.join(source_directory, fn) for fn in h_files]
    for f in h_files2:
        if not os.path.isfile(f):
            raise IOError, "The file '%s' does not exist." % f

    # Check system_headers argument
    system_headers = system_headers or []
    assert isinstance(system_headers, list), "Provide a 'list' as 'system_headers'"
    assert all(isinstance(header, str) for header in system_headers), \
           "Elements of 'system_headers' must be 'str'."

    system_headers.append("boost/shared_ptr.hpp")

    # Get the swig interface file declarations
    declarations = extract_declarations(h_files2)

    # Check system requirements
    (cpp_path, swig_include_dirs, library_dirs, libraries) = \
               configure_instant(swig_binary, swig_path)
    
    # Call instant and return module
    return instant.build_module(wrap_headers            = h_files,
                                source_directory        = source_directory,
                                additional_declarations = declarations,
                                system_headers          = system_headers,
                                include_dirs            = cpp_path,
                                library_dirs            = library_dirs,
                                libraries               = libraries,
                                swigargs                = ['-c++', '-I.','-O'],
                                swig_include_dirs       = swig_include_dirs,
                                **kwargs)

def configure_instant(swig_binary="swig", swig_path=""):
    "Check system requirements"

    # Get include directory for ufc.h (might be better way to do this?)
    (path, dummy, dummy, dummy) = instant.header_and_libs_from_pkgconfig("ufc-1")
    if len(path) == 0: path = [(os.sep).join(sysconfig.get_python_inc().\
                                    split(os.sep)[:-2]) + os.sep + "include"]

    # Register the paths
    cpp_path, swig_include_dirs = [path[0]], [path[0]]

    # Check for swig installation
    if not instant.check_and_set_swig_binary(swig_binary, swig_path):
        raise OSError, "Could not find swig installation. Pass an existing "\
              "swig binary or install SWIG version 1.3.35 or higher.\n"

    # Check swig version for shared_ptr
    if not instant.check_swig_version("1.3.35"):
        raise OSError, "Your current swig version is %s, it needs to be "\
              "1.3.35 or higher.\n" % instant.get_swig_version()

    # Check if UFC is importable and what version of swig was used to
    # create the UFC extension module
    try:
        import ufc
    except:
        raise OSError, "Please install the python extenstion module of UFC on your system.\n"

    # Check that the form compiler will use the same swig version
    # that UFC was compiled with
    if not instant.check_swig_version(ufc.__swigversion__, same=True):
        raise OSError, """The python extension module of UFC was not compiled with the present version of swig.
Install swig version %s or recompiled UFC with present swig
""" % ufc.__swigversion__

    # Set a default swig command and boost include directory
    boost_include_dir = []

    # Check for boost installation
    # Set a default directory for the boost installation
    if sys.platform == "darwin":
        # Use MacPorts as default
        default = os.path.join(os.path.sep, "opt", "local")
    else:
        default = os.path.join(os.path.sep, "usr")

    # If BOOST_DIR is not set use default directory
    boost_dir = os.getenv("BOOST_DIR", default)
    boost_is_found = False
    for inc_dir in ["", "include"]:
        if os.path.isfile(os.path.join(boost_dir, inc_dir, "boost", "version.hpp")):
            boost_include_dir = [os.path.join(boost_dir, inc_dir)]
            boost_is_found = True
            break

    if not boost_is_found:
        raise OSError, """The Boost headers was not found.
If Boost is installed in a nonstandard location,
set the environment variable BOOST_DIR.
"""
    
    # Add the boost_include_dir
    cpp_path += boost_include_dir
    
    # Check for boost_math library
    # FIXME: This is a hack and should be done properly, probably using
    # FIXME: cmake --find-packages
    # If BOOST_DIR is not set use default directory
    boost_dir = os.getenv("BOOST_DIR", default)
    boost_math_is_found = False
    for lib_dir in ["", "lib64", "lib"]:
        for math_lib in ["-mt", ""]:
            if glob.glob(os.path.join(boost_dir, lib_dir, \
                                      "libboost_math_tr1%s*"%math_lib)):
                boost_library_dir = [os.path.join(boost_dir, lib_dir)]
                boost_math_is_found = True
                boost_library = ["boost_math_tr1%s"%math_lib]
                break
    
    if not boost_math_is_found:
        raise OSError, """The Boost math library was not found.
If Boost math library is installed in a nonstandard location,
set the environment variable BOOST_DIR.
"""
    return cpp_path, swig_include_dirs, boost_library_dir, boost_library

def extract_declarations(h_files):
    "Extract information for shared_ptr"

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

// Swig shared_ptr macro declarations
"""

    for h_file in h_files:
        # Read the code
        code = open(h_file).read()

        # Extract the class names
        derived_classes   = re.findall(r"class[ ]+([\w]+)[ ]*: public", code)
        ufc_classes       = re.findall(r"public[ ]+(ufc::[\w]+).*", code)
        ufc_proxy_classes = [s.replace("ufc::", "") for s in ufc_classes]

        shared_ptr_format = "SWIG_SHARED_PTR_DERIVED(%(der_class)s,%(ufc_class)s,%(der_class)s)"
        new_share_ptr_format = "%%shared_ptr(%s)"

        # Write shared_ptr code for swig 2.0.0 or higher
        declarations += "\n"
        declarations += "#if SWIG_VERSION >= 0x020000\n"
        declarations += "\n".join(new_share_ptr_format%c for c in derived_classes)

        declarations += "\n"
        declarations += "#else\n"
        declarations += "\n"
        declarations += "\n".join(\
            shared_ptr_format % { "ufc_proxy_class": c[0], "ufc_class": c[1], "der_class": c[2] }\
            for c in zip(ufc_proxy_classes, ufc_classes, derived_classes)\
            )
        declarations += "\n"
        declarations += "\n"
        declarations += "#endif\n"
        declarations += "\n"
    return declarations
