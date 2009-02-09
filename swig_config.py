# Check for swig, swig version and boost installation
__all__ = ['configure_swig_env','get_status_output']
import os, sys, re

def configure_swig_env(swig_env):
    " Configure and check swig dependencies" 
    conf = swig_env.Configure(custom_tests={"checkSwig"  : checkSwig,
                                            "checkSwigVersion"  : checkSwigVersion,
                                            "checkBoost" : checkBoost,
                                            "checkPython": checkPython})
    if not conf.checkSwig():
        swig_env["enablePyUFC"] = False
        swig_env = conf.Finish()
        return swig_env, """*** Warning: swig not found.
Python extension module will not be installed.
"""
    
    if not conf.checkSwigVersion():
        swig_env["enablePyUFC"] = False
        swig_env = conf.Finish()
        return swig_env, """*** Warning: swig version must be > 1.3.35.
Python extension module will not be installed.
"""
    
    if not conf.checkBoost():
        swig_env["enablePyUFC"] = False
        swig_env = conf.Finish()
        return swig_env, """*** Warning: Boost installation not found.
Set boostDir option, or swig verison > 1.3.35 not found.
Python extension module will not be installed.
"""

    if not conf.checkPython():
        swig_env["enablePyUFC"] = False
        swig_env = conf.Finish()
        return swig_env, """
*** Warning: Unable to find Python development files on your system.
Perhaps you need to install the package python%s-dev?
Python extension module will not be installed.
""" % swig_env["PYTHON_VERSION"]

    # If the dependecies are met, return a prepared swig environment
    return prepareSwigEnv(swig_env), ""

def checkSwig(context):
    context.Message("Checking for SWIG...")
    result, output = get_status_output("swig -version")
    context.Result(not result)
    if result != 0:
        return False
    return True

def checkSwigVersion(context):
    context.Message("Checking for SWIG version > 1.3.35...")
    result, output = get_status_output("swig -version")
    swig_version = re.findall(r"SWIG Version ([0-9.]+)",output)[0]
    swig_enough = True
    swig_correct_version = [1,3,35]
    for i, v in enumerate([int(v) for v in swig_version.split(".")]):
        if swig_correct_version[i] < v:
            break
        elif swig_correct_version[i] == v:
            continue
        else:
            swig_enough = False
    context.Result(swig_enough)
    return swig_enough

def checkPython(context):
    context.Message("Checking for Python...")
    # Check for python development files
    python_version = "%s.%s" % (sys.version_info[0],sys.version_info[1])
    context.env["PYTHON_VERSION"] = python_version
    if context.env["PLATFORM"].startswith("win"):
        python_inc_dir = os.path.join(sys.prefix,"include")
        python_lib_dir = os.path.join(sys.prefix,"libs")
        python_lib = "python%s%s" % (sys.version_info[0],sys.version_info[1])
    else:
        python_inc_dir = os.path.join(sys.prefix,"include","python"+python_version)

    python_is_found = os.path.isfile(os.path.join(python_inc_dir, "Python.h"))
    
    # Add path to python if found
    if python_is_found:
        context.env.Append(CPPPATH=[python_inc_dir])
        if context.env["PLATFORM"].startswith("win"):
            context.env.Append(LIBPATH=[python_lib_dir])
            context.env.Append(LIBS=[python_lib])
    context.Result(python_is_found)
    return python_is_found


def checkBoost(context):
    context.Message("Checking for Boost...")
    # Set a default directory for the boost installation
    if context.env["PLATFORM"] == "darwin":
        # use fink as default
        default = os.path.join(os.path.sep,"sw")
    elif context.env["PLATFORM"].startswith("win"):
        default = r'C:\local'
    else:
        default = os.path.join(os.path.sep,"usr")
            
    boost_dir = context.env.get("boostDir",os.getenv("BOOST_DIR",default))
    boost_is_found = False
    
    for inc_dir in ["", "include"]:
        if os.path.isfile(os.path.join(boost_dir, inc_dir, "boost", "version.hpp")):
            boost_dir = os.path.join(boost_dir, inc_dir)
            boost_is_found = True
            break
    
    # Add path if boost is found
    if boost_is_found:
        context.env.Append(CPPPATH = [boost_dir])
    
    context.Result(boost_is_found)
    return boost_is_found

def prepareSwigEnv(swig_env):
    # Determine which compiler to be used:
    cxx_compilers = ["c++", "g++", "CC"]
    
    # Use CXX from os.environ if available:
    swig_env["CXX"] = os.environ.get("CXX", swig_env.Detect(cxx_compilers))
    
    if not swig_env["CXX"]:
        print "Unable to find any valid C++ compiler."
        # try to use g++ as default:
        swig_env["CXX"] = "g++"

    if "-Werror" in swig_env["CXXFLAGS"]:
            swig_env["CXXFLAGS"].remove("-Werror")
    
    swig_env["SWIGFLAGS"] = "-python -c++ -shadow".split()
    swig_env["SHLIBPREFIX"] = "_"
    if swig_env["PLATFORM"] == "darwin":
        swig_env.Append(CXXFLAGS=" -undefined dynamic_lookup")
        swig_env.Append(SHLINKFLAGS=" -undefined dynamic_lookup")
        swig_env.Append(LDMODULEFLAGS=" -undefined dynamic_lookup")
        if not swig_env.GetOption("clean"):
            # remove /usr/lib if present in LIBPATH.
            # there should maybe be some testing here to make sure 
            # it doesn't blow in our face.
            try: swig_env["LIBPATH"].remove('/usr/lib')
            except: pass
        swig_env['SHLINKFLAGS'] = swig_env['LDMODULEFLAGS']
        swig_env['SHLINKCOM']   = swig_env['LDMODULECOM']
        swig_env["SHLIBSUFFIX"] = ".so"
    elif swig_env["PLATFORM"].startswith("win"):
        swig_env["SHLIBSUFFIX"] = ".pyd"

    swig_env.Append(CPPPATH = [os.path.join("src","ufc")])
    return swig_env    

# Taken from http://ivory.idyll.org/blog/mar-07/replacing-commands-with-subprocess
def get_status_output(cmd, input=None, cwd=None, env=None):
    from subprocess import Popen, PIPE, STDOUT
    pipe = Popen(cmd, shell=True, cwd=cwd, env=env, stdout=PIPE, stderr=STDOUT)
    (output, errout) = pipe.communicate(input=input)
    assert not errout
    status = pipe.returncode
    return (status, output)

