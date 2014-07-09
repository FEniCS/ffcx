#!/usr/bin/env python

from distutils.core import setup
from distutils import sysconfig
import sys
import os
import platform

# Version number
major = 1
minor = 4
maintenance = 0
#development = "" # Select this for release
development = "+" # Select this otherwise
version = "%d.%d.%d%s" % (major, minor, maintenance, development),

packages = [
    "uflacs",
    "uflacs.commands",
    "uflacs.utils",
    "uflacs.codeutils",
    "uflacs.analysis",
    "uflacs.elementtables",
    "uflacs.generation",
    "uflacs.representation",
    "uflacs.backends",
    "uflacs.backends.ffc",
    ]


scripts = [
    "uflacs",
    ]
scripts = [os.path.join("scripts", script) for script in scripts]


if platform.system() == "Windows" or "bdist_wininst" in sys.argv:
    # In the Windows command prompt we can't execute Python scripts
    # without a .py extension. A solution is to create batch files
    # that runs the different scripts.
    batch_files = []
    for script in scripts:
        batch_file = script + ".bat"
        with open(batch_file, "w") as f:
            f.write('python "%%~dp0\%s" %%*' % os.path.split(script)[1])
        batch_files.append(batch_file)
    scripts.extend(batch_files)


CLASSIFIERS = """
Development Status :: 3 - Alpha
Environment :: Console
Intended Audience :: Developers
Intended Audience :: Science/Research
Programming Language :: Python :: 2.5
License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)
Topic :: Scientific/Engineering :: Mathematics
Topic :: Software Development :: Compilers
Topic :: Software Development :: Libraries :: Python Modules
"""
classifiers = CLASSIFIERS.split('\n')[1:-1]


setup(name="uflacs",
      version=version,
      description="UFL Analyser and Compiler System",
      author="Martin Sandve Alnæs",
      author_email="martinal@simula.no",
      url="http://bitbucket.com/fenics-project/uflacs",
      classifiers=classifiers,
      scripts=scripts,
      packages=packages,
      package_dir={"uflacs": "uflacs"},
#     data_files=[(os.path.join("share", "man", "man1"),
#                  [os.path.join("doc", "man", "man1", "uflacs.1.gz")])]
    )
