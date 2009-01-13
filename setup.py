#!/usr/bin/env python

from distutils.core import setup
from distutils import sysconfig
from os.path import join as pjoin
import sys

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

setup(name = "UFC",
      version = "%d.%d" % (major, minor),
      description = "Unified Form-assembly Code",
      author = "Martin Sandve Alnaes, Hans Petter Langtangen, Anders Logg, Kent-Andre Mardal and Ola Skavhaug",
      author_email = "ufc@fenics.org",
      url = "http://www.fenics.org/ufc/",
      packages = ["ufc"],
      package_dir = {"ufc": pjoin("src","utils","python","ufc")},
      data_files = [("include", [pjoin("src","ufc","ufc.h")]),
                    (pjoin("lib","pkgconfig"), ["ufc-%d.pc" % major]),
                    (pjoin("include","swig"),[pjoin("src","ufc","ufc.i")])])
