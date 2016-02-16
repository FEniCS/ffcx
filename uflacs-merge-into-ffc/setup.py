#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import sys
import re

if sys.version_info < (2, 7):
    print("Python 2.7 or higher required, please upgrade.")
    sys.exit(1)

version = re.findall('__version__ = "(.*)"',
                     open('uflacs/__init__.py', 'r').read())[0]

packages = [
    "uflacs",
    "uflacs.language",
    "uflacs.datastructures",
    "uflacs.analysis",
    "uflacs.elementtables",
    "uflacs.generation",
    "uflacs.representation",
    "uflacs.backends",
    "uflacs.backends.ffc",
    "uflacs.backends.ufc",
    ]


CLASSIFIERS = """
Development Status :: 3 - Alpha
Environment :: Console
Intended Audience :: Developers
Intended Audience :: Science/Research
Programming Language :: Python :: 2.7
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
      install_requires = ["numpy", "six", "ufl==1.7.0dev"],
      #data_files=[(os.path.join("share", "man", "man1"),
      #             [os.path.join("doc", "man", "man1", "uflacs.1.gz")])]
    )
