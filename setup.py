#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils import sysconfig
from os.path import join as pjoin
from os.path import pardir
import sys, os

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

# Checking for boost installation uses pkg-config
# If no boost installation is found compiled the ufc extension module
# without shared_ptr support

# Set pkg-config to always include system flags
env = os.environ.copy()
try:
    assert env["PKG_CONFIG_ALLOW_SYSTEM_CFLAGS"] == "0"
except:
    env["PKG_CONFIG_ALLOW_SYSTEM_CFLAGS"] = "1"

# Check if pkg-config is installed
result, output = get_status_output("pkg-config --version ")
if result != 0: 
    raise OSError("The pkg-config package is not installed on the system.")

# Check if boost is installed
result, output = get_status_output("pkg-config --exists boost")
if result == 0:
    boost_include_dir = [get_status_output("pkg-config --cflags-only-I boost",env=env)[1][2:]]
    swig_cmd = "swig -python -c++ -fcompact -O -I. -small ufc.i"
else:
    print "*** Warning: Boost was not found. Python extension module compiled without shared_ptr support"
    boost_include_dir = []
    swig_cmd = "swig -python -c++ -fcompact -O -I. -DNO_SHARED_PTR -small ufc.i"
    
# Run swig on the ufc.i file
os.chdir(pjoin("src","ufc"))
print "Running swig command"
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
