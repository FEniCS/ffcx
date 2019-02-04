# -*- coding: utf-8 -*-
import os
import sys
import subprocess
import string
import setuptools

if sys.version_info < (3, 5):
    print("Python 3.5 or higher required, please upgrade.")
    sys.exit(1)

on_rtd = os.environ.get('READTHEDOCS') == 'True'

VERSION = "2018.2.0.dev0"
RESTRICT_REQUIREMENTS = ">=2018.2.0.dev0,<2018.3"

if on_rtd:
    REQUIREMENTS = []
else:
    REQUIREMENTS = [
        "numpy",
        "cffi",
        "fenics-fiat{}".format(RESTRICT_REQUIREMENTS),
        "fenics-ufl{}".format(RESTRICT_REQUIREMENTS),
    ]

URL = "https://bitbucket.org/fenics-project/ffc/"

ENTRY_POINTS = {'console_scripts': ['ffc = ffc.__main__:main', 'ffc-3 = ffc.__main__:main']}

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
Programming Language :: Python :: 3
Programming Language :: Python :: 3.5
Programming Language :: Python :: 3.6
Topic :: Scientific/Engineering :: Mathematics
Topic :: Software Development :: Libraries :: Python Modules
Topic :: Software Development :: Code Generators
"""


def tarball():
    if "dev" in VERSION:
        return None
    return URL + "downloads/fenics-ffc-{}.tar.gz".format(VERSION)


def get_git_commit_hash():
    """Return git commit hash of currently checked out revision
    or "unknown"
    """
    try:
        hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'])
    except (OSError, subprocess.CalledProcessError) as e:
        print('Retrieving git commit hash did not succeed with exception:')
        print('"{}"'.format(e))
        print()
        print('Stored git commit hash will be set to "unknown"!')
        return "unknown"
    else:
        return hash.decode("ascii").strip()


def write_config_file(infile, outfile, variables={}):
    """Write config file based on template"""

    class AtTemplate(string.Template):
        delimiter = "@"

    s = AtTemplate(open(infile, "r").read())
    s = s.substitute(**variables)
    with open(outfile, "w") as a:
        a.write(s)


def generate_git_hash_file(GIT_COMMIT_HASH):
    """Generate module with git hash"""
    write_config_file(
        os.path.join("ffc", "git_commit_hash.py.in"),
        os.path.join("ffc", "git_commit_hash.py"),
        variables=dict(GIT_COMMIT_HASH=GIT_COMMIT_HASH))


def run_install():
    """Run installation"""

    # Get common variables
    GIT_COMMIT_HASH = get_git_commit_hash()

    # Scripts list
    entry_points = ENTRY_POINTS

    # Generate module with git hash from template
    generate_git_hash_file(GIT_COMMIT_HASH)

    # Call distutils to perform installation
    setuptools.setup(
        name="fenics-ffc",
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
        packages=[
            "ffc",
            "ffc.codegeneration",
            "ffc.codegeneration.C",
            "ffc.ir",
            "ffc.ir.uflacs",
            "ffc.ir.uflacs.analysis",
        ],
        package_dir={"ffc": "ffc"},
        package_data={"ffc": [os.path.join('codegeneration', '*.h')]},
        #scripts=scripts,  # Using entry_points instead
        entry_points=entry_points,
        install_requires=REQUIREMENTS,
        zip_safe=False)


if __name__ == "__main__":
    run_install()
