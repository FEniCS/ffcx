__author__ = "Johan Hake (hake@simula.no)"
__date__ = "2009-03-06 -- 2014-05-20"
__license__  = "This code is released into the public domain"

__all__ = ['build_ufc_module']

# Modified by Martin Alnes, 2009

import instant
import os, sys, re, glob

from distutils import sysconfig

def build_ufc_module(h_files, source_directory="", system_headers=None, \
                     **kwargs):
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
    assert all(isinstance(f, str) for f in h_files), \
           "Elements of 'h_files' must be 'str'."

    h_files2 = [os.path.join(source_directory, fn) for fn in h_files]
    for f in h_files2:
        if not os.path.isfile(f):
            raise IOError("The file '%s' does not exist." % f)

    # Check system_headers argument
    system_headers = system_headers or []
    assert isinstance(system_headers, list), "Provide a 'list' as 'system_headers'"
    assert all(isinstance(header, str) for header in system_headers), \
           "Elements of 'system_headers' must be 'str'."

    system_headers.append("memory")

    # Get the swig interface file declarations
    declarations = extract_declarations(h_files2)
    declarations += """

// SWIG version
%inline %{
int get_swigversion() { return  SWIGVERSION; }
%}

%pythoncode %{
tmp = hex(get_swigversion())
swigversion = "%d.%d.%d"%(tuple(map(int, [tmp[-5], tmp[-3], tmp[-2:]])))
del tmp, get_swigversion
%}
"""

    # Call instant and return module
    return instant.build_module(wrap_headers            = h_files,
                                source_directory        = source_directory,
                                additional_declarations = declarations,
                                system_headers          = system_headers,
                                cmake_packages          = ["UFC"],
                                **kwargs)

def extract_declarations(h_files):
    "Extract information for shared_ptr"

    # Swig declarations
    declarations =r"""
// Use std::shared_ptr
#define SWIG_SHARED_PTR_NAMESPACE std
%include <std_shared_ptr.i>

// Swig shared_ptr macro declarations
"""

    for h_file in h_files:
        # Read the code
        with open(h_file) as file:
            code = file.read()

        # Extract the class names
        derived_classes   = re.findall(r"class[ ]+([\w]+)[ ]*: public", code)
        ufc_classes       = re.findall(r"public[ ]+(ufc::[\w]+).*", code)
        ufc_proxy_classes = [s.replace("ufc::", "") for s in ufc_classes]

        new_share_ptr_format = "%%shared_ptr(%s)"

        # Write shared_ptr code for swig 2.0.0 or higher
        declarations += "\n".join(new_share_ptr_format%c for c in derived_classes)
        declarations += "\n"
    return declarations
