#!/usr/bin/env python

from distutils.core import setup
from os import chdir

chdir("src")

setup(name = "FFC",
      version = "0.3.5-dev",
      description = "The FEniCS Form Compiler",
      author = "Anders Logg",
      author_email = "logg@simula.no",
      url = "http://www.fenics.org/ffc/",
      packages = ["ffc/",
                  "ffc/common",
                  "ffc/fem",
                  "ffc/compiler",
                  "ffc/compiler/optimization",
                  "ffc/compiler/format",
                  "ffc/compiler/codegeneration",
                  "ffc/compiler/codegeneration/common",
                  "ffc/compiler/codegeneration/quadrature",
                  "ffc/compiler/codegeneration/tensor",
                  "ffc/compiler/representation",
                  "ffc/compiler/representation/quadrature",
                  "ffc/compiler/representation/tensor",
                  "ffc/compiler/language",
                  "ffc/compiler/analysis",
                  "ffc/old"],
      scripts = ["bin/ffc"],
      data_files = [("share/man/man1", ["../doc/man/man1/ffc.1.gz"])])
