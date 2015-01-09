#!/usr/bin/env python

from distutils.core import setup
from distutils import sysconfig
import sys
import os
import platform

if sys.version_info < (2, 7):
    print("Python 2.7 or higher required, please upgrade.")
    sys.exit(1)

# Version number
major = 1
minor = 4
maintenance = 0
#development = "" # Select this for release
development = "+" # Select this otherwise
version = "%d.%d.%d%s" % (major, minor, maintenance, development)


packages = [
    "uflacs",
    "uflacs.codeutils",
    "uflacs.datastructures",
    "uflacs.analysis",
    "uflacs.elementtables",
    "uflacs.generation",
    "uflacs.representation",
    "uflacs.backends",
    "uflacs.backends.ffc",
    ]


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
      author="Martin Sandve Alnaes",
      author_email="martinal@simula.no",
      url="http://bitbucket.com/fenics-project/uflacs",
      classifiers=classifiers,
      packages=packages,
      package_dir={"uflacs": "uflacs"},
#     data_files=[(os.path.join("share", "man", "man1"),
#                  [os.path.join("doc", "man", "man1", "uflacs.1.gz")])]
    )
