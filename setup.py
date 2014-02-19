#!/usr/bin/env python

from distutils.core import setup
from distutils import sysconfig
import sys
import os
import platform

# Version number
major = 1
minor = 3
maintenance = 0
#development = "" # Select this for release
development = "+" # Select this otherwise


packages = [
    "uflacs",
    "uflacs.commands",
    "uflacs.utils",
    "uflacs.codeutils",
    "uflacs.geometry",
    "uflacs.analysis",
    "uflacs.generation",
    "uflacs.backends",
    "uflacs.backends.toy",
    "uflacs.backends.latex",
    "uflacs.backends.ffc",
    "uflacs.backends.dolfin",
    ]


scripts = [os.path.join("scripts", "uflacs")]
if platform.system() == "Windows" or "bdist_wininst" in sys.argv:
    # In the Windows command prompt we can't execute Python scripts
    # without a .py extension. A solution is to create batch files
    # that runs the different scripts.
    batch_files = []
    for script in scripts:
        batch_file = script + ".bat"
        f = open(batch_file, "w")
        f.write('python "%%~dp0\%s" %%*' % os.path.split(script)[1])
        f.close()
        batch_files.append(batch_file)
    scripts.extend(batch_files)


CLASSIFIERS = """
Development Status :: 2 - Pre-Alpha
Environment :: Console
Intended Audience :: Developers
Intended Audience :: Science/Research
Programming Language :: Python :: 2.5
License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)
Topic :: Scientific/Engineering :: Mathematics
Topic :: Software Development :: Compilers
Topic :: Software Development :: Libraries :: Python Modules
Topic :: Utilities
"""

setup(name = "uflacs",
      version = "%d.%d.%d%s" % (major, minor, maintenance, development),
      description = "UFL Analyser and Compiler System",
      author = "Martin Sandve Alnaes",
      author_email = "martinal@simula.no",
      url = 'http://bitbucket.com/fenics-project/uflacs',
      classifiers = CLASSIFIERS.split('\n')[1:-1],
      scripts = scripts,
      packages = packages,
      package_dir = {"uflacs": "uflacs"},
#     data_files = [(os.path.join("share", "man", "man1"),
#                    [os.path.join("doc", "man", "man1", "uflacs.1.gz")])]
    )

