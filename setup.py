#!/usr/bin/env python

from distutils.core import setup
from os import chdir
from os.path import join

setup(name = "FFC",
      version = "0.4.1",
      description = "The FEniCS Form Compiler",
      author = "Anders Logg",
      author_email = "logg@simula.no",
      url = "http://www.fenics.org/ffc/",
      packages = ["ffc",
                  "ffc.common",
                  "ffc.fem",
                  "ffc.compiler",
                  "ffc.compiler.optimization",
                  "ffc.compiler.format",
                  "ffc.compiler.codegeneration",
                  "ffc.compiler.codegeneration.common",
                  "ffc.compiler.codegeneration.quadrature",
                  "ffc.compiler.codegeneration.tensor",
                  "ffc.compiler.representation",
                  "ffc.compiler.representation.quadrature",
                  "ffc.compiler.representation.tensor",
                  "ffc.compiler.language",
                  "ffc.compiler.analysis"],
      package_dir={"ffc": join("src", "ffc")},
      scripts = [join("src","bin","ffc")],
      data_files = [(join("share", "man", "man1"), [join("doc", "man", "man1", "ffc.1.gz")])])
