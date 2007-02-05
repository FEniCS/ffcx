#!/usr/bin/env python

from distutils.core import setup, Extension

setup(name = "UFC",
      version = "1.0",
      description = "Unified Form-assembly Code",
      author = "Martin Sandve Alnaes, Hans Petter Langtangen, Anders Logg, Kent-Andre Mardal and Ola Skavhaug",
      url = "http://www.fenics.org/ufc/",
      packages = ["ufc"],
      package_dir= {"ufc": "src/utils/python/ufc"},
      ext_modules = [Extension("ufc", ["src/ufc/ufc.c"])])
