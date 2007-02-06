#!/usr/bin/env python

from distutils.core import setup, Extension
import sys

try:    prefix = [ item for item in sys.argv[1:] if "--prefix=" in item ][0].split("=")[1]
except: prefix = "/usr/local"

setup(name = "UFC",
      version = "1.0",
      description = "Unified Form-assembly Code",
      author = "Martin Sandve Alnaes, Hans Petter Langtangen, Anders Logg, Kent-Andre Mardal and Ola Skavhaug",
      url = "http://www.fenics.org/ufc/",
      packages = ["ufc"],
      package_dir= {"ufc": "src/utils/python/ufc"},
      data_files=[("%s/include" % prefix, ["src/ufc/ufc.h"])]
    )
