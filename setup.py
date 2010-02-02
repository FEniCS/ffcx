#!/usr/bin/env python

import sys, platform
from distutils.core import setup
from os import chdir
from os.path import join, split

scripts = [join("scripts", "ffc"), join("scripts", "ffc-clean")]

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

setup(name = "FFC",
      version = "0.9.0+",
      description = "The FEniCS Form Compiler",
      author = "Anders Logg, Kristian Oelgaard, Marie Rognes et al.",
      author_email = "ffc@lists.launchpad.net",
      url = "http://www.fenics.org/ffc/",
      packages = ["ffc", "ffc.quadrature", "ffc.tensor"],
      package_dir={"ffc": "ffc"},
      scripts = scripts,
      data_files = [(join("share", "man", "man1"),
                     [join("doc", "man", "man1", "ffc.1.gz"),
                      join("doc", "man", "man1", "ffc-clean.1.gz")])])
