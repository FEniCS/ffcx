#!/usr/bin/env python

from distutils.core import setup
from os import chdir

chdir("src")

setup(name="FFC",
      version="0.1.4",
      description="The FEniCS Form Compiler",
      author="Anders Logg",
      author_email="logg@tti-c.org",
      url="http://www.fenics.org/ffc/",
      packages=["ffc/",
                "ffc/common",
                "ffc/compiler",
                "ffc/format",
                "ffc/parser"],
      scripts=["bin/ffc"])
