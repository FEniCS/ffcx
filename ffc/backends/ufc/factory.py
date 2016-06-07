# Code generation format strings for UFC (Unified Form-assembly Code) v. 2016.1.0dev.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2015.

visibility_snippet = """\
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif
"""

factory_header = """\
class %(namespace)s%(classname)s;

extern "C" %(namespace)s%(classname)s * create_%(classname)s();
"""

factory_implementation = """\
extern "C" DLL_EXPORT %(namespace)s%(classname)s * create_%(classname)s()
{
  return new %(namespace)s%(classname)s();
}
"""
