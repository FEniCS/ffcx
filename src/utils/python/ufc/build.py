__author__ = "Johan Hake (hake@simula.no)"
__date__ = "2009-03-06 -- 2009-03-20"
__copyright__ = "Copyright (C) 2009 Johan Hake"
__license__  = "GNU LGPL Version 2.1"

__all__ = ['build_ufc_module']

import instant
import os, sys, re

from distutils import sysconfig

def old_build_ufc_module(signature, h_files,
                     cpp_files=None,
                     system_headers=None,
                     cpp_args=None,
                     cache_dir=None):

    """ Build a python extension module from ufc complient source code

    The compiled module will be imported an returned by the function.
    
    @param signature:
       The signature that identifies the extenstion module
    @param h_files:
       The name(s) of the header files that should be compiled and included in
       the python extension module
    @param cpp_files:
       Optional ufc-cpp files that will be compiled seperately and linked to the
       ufc module.
    @param system_headers:
       Extra headers that will be #included in the generated wrapper file
    @param cache_dir:
       An optional cache dir.
    """
    print "FIXME: old_build_ufc_module is deprecated, use build_ufc_module instead!"
    # Check signature argument
    assert isinstance(signature,str), "Provide a 'str' as 'signature'"
    
    # Check h_files argument
    assert isinstance(h_files,(list,str)) , "Provide a 'list' or a 'str' as 'h_files'"
    if isinstance(h_files,str):
        h_files = [h_files]
    for f in h_files:
        assert isinstance(f,str), "The provided h_file must be a 'str'"
        if not os.path.isfile(f):
            raise IOError, "The file '%s' does not excist"%f

    # Check cpp_files argument
    if not cpp_files is None:
        assert isinstance(cpp_files,(list,str)) , "Provide a 'list' or a 'str' as 'cpp_files'"
        if isinstance(cpp_files,str):
            cpp_files = [cpp_files]
        for f in cpp_files:
            assert isinstance(f,str), "The provided cpp_file must be a 'str'"
            if not os.path.isfile(f):
                raise IOError, "The file '%s' does not excist"%f
        if len(cpp_files) != len(h_files):
            raise TypeError, "The provided number of h_files and cpp_files must be the same."
    else:
        cpp_files = []

    # Check system_headers argument
    system_headers = system_headers or []
    assert isinstance(system_headers,list), "Provide a 'list' as 'system_headers'"
    assert all(isinstance(header,str) for header in system_headers), "Elements of 'system_headers' must be 'str'"

    system_headers.append("boost/shared_ptr.hpp")

    # Check cpp_args
    cpp_args = cpp_args or []
    
    # Get the swig interface file declarations
    declarations = extract_declarations(h_files)
    
    # Check system requirements
    (cpp_path, swig_path) = configure_instant()
    
    # Call instant and return module
    return instant.build_module(wrap_headers            = h_files,
                                sources                 = cpp_files,
                                additional_declarations = declarations,
                                system_headers          = system_headers,
                                include_dirs            = cpp_path,
                                swigargs                = ['-c++','-I.'],
                                swig_include_dirs       = swig_path,
                                cppargs                 = cpp_args,
                                signature               = signature,
                                cache_dir               = cache_dir)

def build_ufc_module(h_files, source_directory="", system_headers=None, **kwargs):
    """Build a python extension module from ufc compliant source code.
    
    The compiled module will be imported and returned by the function.
    
    @param h_files:
       The name(s) of the header files that should be compiled and included in
       the python extension module.
    @param source_directory:
       The directory where the source files reside.
    @param system_headers:
       Extra headers that will be #included in the generated wrapper file.
    
    Any additional keyword arguments are passed on to instant.build_module.
    """

    # Check h_files argument
    if isinstance(h_files, str):
        h_files = [h_files]
    assert isinstance(h_files, list) , "Provide a 'list' or a 'str' as 'h_files'."
    assert all(isinstance(f, str) for f in h_files), "Elements of 'h_files' must be 'str'."

    h_files2 = [os.path.join(source_directory, fn) for fn in h_files]
    for f in h_files2:
        if not os.path.isfile(f):
            raise IOError, "The file '%s' does not exist." % f
    
    # Check system_headers argument
    system_headers = system_headers or []
    assert isinstance(system_headers, list), "Provide a 'list' as 'system_headers'"
    assert all(isinstance(header, str) for header in system_headers), "Elements of 'system_headers' must be 'str'."
    
    system_headers.append("boost/shared_ptr.hpp")
    
    # Get the swig interface file declarations
    declarations = extract_declarations(h_files2)
    
    # Check system requirements
    (cpp_path, swig_path) = configure_instant()

    # Call instant and return module
    return instant.build_module(wrap_headers            = h_files,
                                source_directory        = source_directory,
                                additional_declarations = declarations,
                                system_headers          = system_headers,
                                include_dirs            = cpp_path,
                                swigargs                = ['-c++', '-I.'],
                                swig_include_dirs       = swig_path,
                                **kwargs)

def configure_instant():
    "Check system requirements"

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
        raise OSError, "Your current swig version is %s, it needs to be 1.3.35 or higher.\n" % swig_version

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
        # use fink as default
        default = os.path.join(os.path.sep, "sw")
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
        raise OSError, """The Boost library was not found.
If Boost is installed in a nonstandard location,
set the environment variable BOOST_DIR.
"""

    # Add the boost_include_dir
    cpp_path += boost_include_dir
    
    return cpp_path, swig_path

def extract_declarations(h_files):
    "Extract information for shared_ptr"

    # Swig declarations
    declarations =r"""
%pythoncode %{
import ufc
'''
A hack to get past a bug in swig.
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
        declarations += "\n"
        declarations += "\n".join(\
            shared_ptr_format % { "ufc_proxy_class": c[0], "ufc_class": c[1], "der_class": c[2] }\
            for c in zip(ufc_proxy_classes, ufc_classes, derived_classes)\
            )

        declarations += "\n\n// Swig %ignore declations\n"
        
        # Extract any provided form, elementa and dofmap classes
        form_classes    = []
        element_classes = []
        dof_map_classes = []
        for i, ufc_class in enumerate(ufc_classes):
            if ufc_class == "ufc::form" :
                form_classes.append(derived_classes[i])
            if ufc_class == "ufc::finite_element" :
                element_classes.append(derived_classes[i])
            if ufc_class == "ufc::dof_map" :
                dof_map_classes.append(derived_classes[i])
        
        # %ignore all foo::create_foo methods
        for form_class in form_classes:
            declarations += """
%%ignore %s::create_finite_element;
%%ignore %s::create_dof_map;
%%ignore %s::create_cell_integral;
%%ignore %s::create_exterior_facet_integral;
%%ignore %s::create_interior_facet_integral;
""" % ((form_class,)*5)
        for element_class in element_classes:
            declarations += """
%%ignore %s::create_sub_element;
""" % (element_class,)
        for dof_map_class in dof_map_classes:
            declarations += """
%%ignore %s::create_sub_dof_map;
""" % (dof_map_class,)
        
        declarations += "\n"
    return declarations

