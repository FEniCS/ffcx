#!/usr/bin/env python
from distutils.core import setup
from distutils.core import Extension
import os

# the buggy swig-support in distutils doesn't manage to invoke g++, uses gcc...
os.system("make ufc_benchmark_wrap.cxx")
extension = Extension('_ufc_benchmark', ['ufc_benchmark.cpp', 'ufc_benchmark_wrap.cxx'], language="c++", include_dirs=["../../../ufc"])

setup(### metadata:
      name              = 'ufc_benchmark',
      version           = '1.1.2',
      author            = 'Martin Sandve Alnes',
      author_email      = 'martinal@simula.no',
      maintainer        = 'Martin Sandve Alnes',
      maintainer_email  = 'martinal@simula.no',
      url               = 'http://www.fenicsproject.org',
      description       = 'Benchmark utility for UFC implementations.',
      download_url      = 'https://launchpad.net/ufc',
      ### contents:
      py_modules   = ['ufc_benchmark'],
      ext_modules  = [extension],
      )

