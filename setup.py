#!/usr/bin/env python

import sys, platform
from distutils.core import setup
from os import chdir
from os.path import join, splitext

scripts = [join("scripts", "ffc"), join("scripts", "ffc-clean")]

if platform.system() == "Windows" or "bdist_wininst" in sys.argv:

    # In the Windows command prompt we can't execute Python scripts 
    # without the .py extension. A solution is to create batch files
    # that runs the different scripts.

    # Try to determine the installation prefix
    if platform.system() == "Windows":
        prefix = sys.prefix
    else:
        # We are running bdist_wininst on a non-Windows platform
        pymajor, pyminor = sysconfig.get_python_version().split(".")
        prefix = "C:\\Python%s%s" % (pymajor, pyminor)

    # If --prefix is specified we use this instead of the default:
    for arg in sys.argv:
        if "--prefix" in arg:
            prefix = arg.split("=")[1]
            break

    # Create batch files for Windows
    for batch_file in ["ffc.bat", "ffc-clean.bat"]:
        f = open(batch_file, "w")
        f.write("@python %s %%*" % join(prefix, "Scripts", splitext(batch_file)[0]))
        f.close()
        scripts.append(batch_file)

setup(name = "FFC",
      version = "0.6.1",
      description = "The FEniCS Form Compiler",
      author = "Anders Logg et al.",
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
                  "ffc.compiler.analysis",
                  "ffc.jit"],
      package_dir={"ffc": "ffc"},
      scripts = scripts,
      data_files = [(join("share", "man", "man1"),
                     [join("doc", "man", "man1", "ffc.1.gz"),
                      join("doc", "man", "man1", "ffc-clean.1.gz")])])
