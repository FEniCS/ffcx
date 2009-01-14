#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils import sysconfig
from os.path import join as pjoin
from os.path import pardir
import sys, os, re

# Taken from http://ivory.idyll.org/blog/mar-07/replacing-commands-with-subprocess
from subprocess import Popen, PIPE, STDOUT
def get_status_output(cmd, input=None, cwd=None, env=None):
    pipe = Popen(cmd, shell=True, cwd=cwd, env=env, stdout=PIPE, stderr=STDOUT)
    (output, errout) = pipe.communicate(input=input)
    assert not errout
    status = pipe.returncode
    return (status, output)

# Version number
major = 1
minor = 1

# Set prefix
try:
    prefix = [item for item in sys.argv[1:] \
              if "--prefix=" in item][0].split("=")[1]
except:
    try:
        prefix = sys.argv[sys.argv.index('--prefix')+1]
    except:
        prefix = sys.prefix

print "Installing UFC under %s..." % prefix

# Generate pkgconfig file
file = open("ufc-%d.pc" % major, "w")
file.write("Name: UFC\n")
file.write("Version: %d.%d\n" % (major, minor))
file.write("Description: Unified Form-assembly Code\n")
file.write("Cflags: -I%s\n" % repr(pjoin(prefix,"include"))[1:-1])
# FIXME: better way for this? ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
file.close()

# Check for swig installation
result, output = get_status_output("swig -version")
if result == 1:
    raise OSError, "Could not find swig installation. Please install swig version 1.3.35 or higher"

# Check swig version
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
        print """*** Warning: Your swig version need to be 1.3.35 or higher.
Python extension module will be compiled without shared_ptr support"""
    
# Set a default swig command and boost include directory
swig_cmd = "swig -python -c++ -O -I. -DNO_SHARED_PTR ufc.i"
boost_include_dir = []

# If swig version 1.3.35 or higher we check for boost installation
if swig_enough:
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
        if os.path.isfile(pjoin(boost_dir, inc_dir, "boost", "version.hpp")):
            boost_include_dir = [pjoin(boost_dir, inc_dir)]
            boost_is_found = True
            break
    
    if boost_is_found:
        print "Using Boost installation in: %s"%boost_include_dir[0]
        swig_cmd = "swig -python -c++ -O -I. ufc.i"
    else:
        # If no boost installation is found compile the ufc extension module
        # without shared_ptr support
        print """*** Warning: Boost was not found.
--------------------------------------------------------------    
Python extension module compiled without shared_ptr support

If boost is installed on your system you can specify the path
by setting the environment variable BOOST_DIR.
--------------------------------------------------------------"""
    
# Run swig on the ufc.i file
os.chdir(pjoin("src","ufc"))
print "running swig command"
print swig_cmd
os.system(swig_cmd)
os.chdir(os.pardir); os.chdir(os.pardir)
sources = [pjoin("src","ufc","ufc_wrap.cxx")]

setup(name = "UFC",
      version = "%d.%d" % (major, minor),
      description = "Unified Form-assembly Code",
      author = "Martin Sandve Alnaes, Hans Petter Langtangen, Anders Logg, Kent-Andre Mardal and Ola Skavhaug",
      author_email = "ufc@fenics.org",
      url = "http://www.fenics.org/ufc/",
      packages = ["","ufc_utils"],
      package_dir = {"":pjoin("src","ufc"),
                     "ufc_utils": pjoin("src","utils","python","ufc")},
      data_files = [("include", [pjoin("src","ufc","ufc.h")]),
                    (pjoin("lib","pkgconfig"), ["ufc-%d.pc" % major]),
                    (pjoin("include","swig"),[pjoin("src","ufc","ufc.i")])],
      ext_modules = [Extension("_ufc",
                     sources,
                     include_dirs=[pjoin("src","ufc")] + boost_include_dir
                               )])
