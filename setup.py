#!/usr/bin/env python

from distutils.core import setup
from os import chdir

chdir("src")

setup(name="FFC",
      version="0.1.3",
      description="The FEniCS Form Compiler",
      author="Anders Logg",
      author_email="logg@tti-c.org",
      url="http://www.fenics.org/ffc/",
      packages=["ffc/",
                "ffc/parser",
                "ffc/compiler",
                "ffc/format"],
      scripts=["bin/ffc"])
