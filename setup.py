#!/usr/bin/env python

import sys, platform
from distutils.core import setup, Extension
from os import chdir
from os.path import join, split
import numpy

scripts = [join("scripts", "ffc")]

if platform.system() == "Windows" or "bdist_wininst" in sys.argv:
    # In the Windows command prompt we can't execute Python scripts
    # without a .py extension. A solution is to create batch files
    # that runs the different scripts.
    batch_files = []
    for script in scripts:
        batch_file = script + ".bat"
        f = open(batch_file, "w")
        f.write('python "%%~dp0\%s" %%*\n' % split(script)[1])
        f.close()
        batch_files.append(batch_file)
    scripts.extend(batch_files)

ext = Extension("ffc.time_elements_ext",
                ["ffc/ext/time_elements_interface.cpp",
                 "ffc/ext/time_elements.cpp",
                 "ffc/ext/LobattoQuadrature.cpp",
                 "ffc/ext/RadauQuadrature.cpp",
                 "ffc/ext/Legendre.cpp"],
                include_dirs=[numpy.get_include()])

setup(name = "FFC",
      version = "1.2.0",
      description = "The FEniCS Form Compiler",
      author = "Anders Logg, Kristian Oelgaard, Marie Rognes et al.",
      author_email = "ffc@lists.launchpad.net",
      url = "http://www.fenicsproject.org",
      packages = ["ffc",
                  "ffc.quadrature", "ffc.tensor", "ffc.uflacsrepr",
                  "ffc.errorcontrol",
                  "ffc.dolfin"],
      package_dir={"ffc": "ffc"},
      scripts = scripts,
      ext_modules = [ext],
      data_files = [(join("share", "man", "man1"),
                     [join("doc", "man", "man1", "ffc.1.gz")])])
