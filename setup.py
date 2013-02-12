#!/usr/bin/env python

from distutils.core import setup
from distutils import sysconfig
from os.path import join as pjoin, split as psplit
import sys
import platform

# Version number
major = 0
minor = 2
maintenance = 0

scripts = [pjoin("scripts", "uflacs")]

if platform.system() == "Windows" or "bdist_wininst" in sys.argv:
    # In the Windows command prompt we can't execute Python scripts
    # without a .py extension. A solution is to create batch files
    # that runs the different scripts.
    batch_files = []
    for script in scripts:
        batch_file = script + ".bat"
        f = open(batch_file, "w")
        f.write('python "%%~dp0\%s" %%*' % psplit(script)[1])
        f.close()
        batch_files.append(batch_file)
    scripts.extend(batch_files)

setup(name = "uflacs",
      version = "%d.%d.%d" % (major, minor, maintenance),
      description = "UFL Analyser and Compiler System",
      author = "Martin Sandve Alnes",
      author_email = "ufl@lists.launchpad.net",
      url = 'https://launchpad.net/uflacs',
      classifiers=[
          'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python :: 2.5',
          'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Software Development :: Compilers',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Utilities',
          ],
      scripts = scripts,
      packages = ["uflacs",
                  "uflacs.commands",
                  "uflacs.utils",
                  "uflacs.codeutils",
                  "uflacs.geometry",
                  "uflacs.algorithms",
                  "uflacs.backends",
                  "uflacs.backends.latex",
                  "uflacs.backends.cpp2",
                  "uflacs.backends.ffc",
                  "uflacs.backends.dolfin",
                  ],
      package_dir = {"uflacs": "site-packages/uflacs"},
#     data_files = [(pjoin("share", "man", "man1"),
#                    [pjoin("doc", "man", "man1", "uflacs.1.gz")])]
    )

