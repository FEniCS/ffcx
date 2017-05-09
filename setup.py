# -*- coding: utf-8 -*-
from __future__ import print_function

import os
import sys
import subprocess
import string

from setuptools import setup, find_packages

if sys.version_info < (2, 7):
    print("Python 2.7 or higher required, please upgrade.")
    sys.exit(1)

VERSION = "2017.1.0"

URL = "https://bitbucket.org/fenics-project/ffc/"

if sys.version_info[0] == 2:
    ENTRY_POINTS = {'console_scripts': ['ffc = ffc.__main__:main',
                                        'ffc-2 = ffc.__main__:main']}
else:
    ENTRY_POINTS = {'console_scripts': ['ffc = ffc.__main__:main',
                                        'ffc-3 = ffc.__main__:main']}

AUTHORS = """\
Anders Logg, Kristian Oelgaard, Marie Rognes, Garth N. Wells,
Martin Sandve Alnæs, Hans Petter Langtangen, Kent-Andre Mardal,
Ola Skavhaug, et al.
"""

CLASSIFIERS = """\
Development Status :: 5 - Production/Stable
Intended Audience :: Developers
Intended Audience :: Science/Research
License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)
Operating System :: POSIX
Operating System :: POSIX :: Linux
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows
Programming Language :: Python
Programming Language :: Python :: 2
Programming Language :: Python :: 2.7
Programming Language :: Python :: 3
Programming Language :: Python :: 3.4
Programming Language :: Python :: 3.5
Programming Language :: Python :: 3.6
Topic :: Scientific/Engineering :: Mathematics
Topic :: Software Development :: Libraries :: Python Modules
Topic :: Software Development :: Code Generators
"""


def tarball():
    if "dev" in VERSION:
        return None
    return URL + "downloads/ffc-%s.tar.gz" % VERSION


def get_installation_prefix():
    "Get installation prefix"
    prefix = sys.prefix
    for arg in sys.argv[1:]:
        if "--user" in arg:
            import site
            prefix = site.getuserbase()
            break
        elif arg in ("--prefix", "--home", "--install-base"):
            prefix = sys.argv[sys.argv.index(arg) + 1]
            break
        elif "--prefix=" in arg or "--home=" in arg or "--install-base=" in arg:
            prefix = arg.split("=")[1]
            break

    return os.path.abspath(os.path.expanduser(prefix))


def get_git_commit_hash():
    """Return git commit hash of currently checked out revision
    or "unknown"
    """
    try:
        hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'])
    except (OSError, subprocess.CalledProcessError) as e:
        print('Retrieving git commit hash did not succeed with exception:')
        print('"%s"' % str(e))
        print()
        print('Stored git commit hash will be set to "unknown"!')
        return "unknown"
    else:
        return hash.strip()


def write_config_file(infile, outfile, variables={}):
    "Write config file based on template"
    class AtTemplate(string.Template):
        delimiter = "@"
    s = AtTemplate(open(infile, "r").read())
    s = s.substitute(**variables)
    with open(outfile, "w") as a:
        a.write(s)


def generate_git_hash_file(GIT_COMMIT_HASH):
    "Generate module with git hash"
    write_config_file(os.path.join("ffc", "git_commit_hash.py.in"),
                      os.path.join("ffc", "git_commit_hash.py"),
                      variables=dict(GIT_COMMIT_HASH=GIT_COMMIT_HASH))


def run_install():
    "Run installation"

    # Get common variables
    #INSTALL_PREFIX = get_installation_prefix()
    GIT_COMMIT_HASH = get_git_commit_hash()

    # Scripts list
    entry_points = ENTRY_POINTS

    # Generate module with git hash from template
    generate_git_hash_file(GIT_COMMIT_HASH)

    # FFC data files
    data_files = [(os.path.join("share", "man", "man1"),
                  [os.path.join("doc", "man", "man1", "ffc.1.gz")])]

    # Call distutils to perform installation
    setup(name="FFC",
          description="The FEniCS Form Compiler",
          version=VERSION,
          author=AUTHORS,
          classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
          license="LGPL version 3 or later",
          author_email="fenics-dev@googlegroups.com",
          maintainer_email="fenics-dev@googlegroups.com",
          url=URL,
          download_url=tarball(),
          platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
          packages=find_packages("."),
          package_dir={"ffc": "ffc"},
          package_data={"ffc" : [os.path.join('backends', 'ufc', '*.h')]},
          #scripts=scripts,  # Using entry_points instead
          entry_points=entry_points,
          data_files=data_files,
          install_requires=["numpy",
                            "six",
                            "fiat==%s" % VERSION,
                            "ufl==%s" % VERSION,
                            "dijitso==%s" % VERSION],
          zip_safe=False)


if __name__ == "__main__":
    run_install()
