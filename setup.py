#!/usr/bin/env python

from distutils.core import setup
from distutils import sysconfig
import sys

# Version number
major = 1
minor = 0

# Set prefix
try:    prefix = [item for item in sys.argv[1:] if "--prefix=" in item][0].split("=")[1]
except: prefix = ("/").join(sysconfig.get_python_inc().split("/")[:-2])
print "Installing UFC under %s..." % prefix

# Generate pkgconfig file
file = open("ufc-%d.pc" % major, "w")
file.write("Name: UFC\n")
file.write("Version: %d.%d\n" % (major, minor))
file.write("Description: Unified Form-assembly Code\n")
file.write("Cflags: -I%s/include\n" % prefix)
file.close()

setup(name = "UFC",
      version = "%d.%d" % (major, minor),
      description = "Unified Form-assembly Code",
      author = "Martin Sandve Alnaes, Hans Petter Langtangen, Anders Logg, Kent-Andre Mardal and Ola Skavhaug",
      author_email = "ufc@fenics.org",
      url = "http://www.fenics.org/ufc/",
      packages = ["ufc"],
      package_dir = {"ufc": "src/utils/python/ufc/"},
      data_files = [("%s/include" % prefix, ["src/ufc/ufc.h"]),
                    ("%s/lib/pkgconfig" % prefix, ["ufc-%d.pc" % major])])
