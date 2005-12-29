#!/usr/bin/env python

from distutils.core import setup
from os import chdir

chdir("src")

setup(name="FFC",
      version="0.2.5",
      description="The FEniCS Form Compiler",
      author="Anders Logg",
      author_email="logg@tti-c.org",
      url="http://www.fenics.org/ffc/",
      packages=["ffc/",
                "ffc/common",
                "ffc/compiler",
                "ffc/format",
                "ffc/parser"],
      scripts=["bin/ffc"],
      data_files=[("share/man/man1", ["../doc/man/man1/ffc.1.gz"])])
