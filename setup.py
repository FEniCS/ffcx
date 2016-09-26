# -*- coding: utf-8 -*-
from __future__ import print_function

import os
import sys
import platform
import re
import subprocess
import string
import tempfile
import shutil
import hashlib

from setuptools import setup, find_packages
from setuptools.command.install import install

from distutils.ccompiler import new_compiler

if sys.version_info < (2, 7):
    print("Python 2.7 or higher required, please upgrade.")
    sys.exit(1)

VERSION = re.findall('__version__ = "(.*)"',
                     open('ffc/__init__.py', 'r').read())[0]

URL = "https://bitbucket.org/fenics-project/ffc/"

SCRIPTS = [os.path.join("scripts", "ffc")]

AUTHORS = """\
Anders Logg, Kristian Oelgaard, Marie Rognes, Garth N. Wells,
Martin Sandve Alnaes, Hans Petter Langtangen, Kent-Andre Mardal,
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


def get_ufc_signature():
    """Compute SHA-1 hash of ufc.h"""
    with open(os.path.join('ffc', 'backends', 'ufc', 'ufc.h'), 'rb') as f:
        return hashlib.sha1(f.read()).hexdigest()


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


def get_cxx_flags():
    """Return flags needed for compilation of UFC C++11 program"""
    cc = new_compiler()
    CXX = os.environ.get("CXX")
    if CXX:
        cc.set_executables(compiler_so=CXX, compiler=CXX, compiler_cxx=CXX)
    CXX_FLAGS = os.environ.get("CXXFLAGS", "")
    if has_cxx_flag(cc, "-std=c++11"):
        CXX_FLAGS += " -std=c++11"
    elif has_cxx_flag(cc, "-std=c++0x"):
        CXX_FLAGS += " -std=c++0x"
    return CXX_FLAGS


def create_windows_batch_files(scripts):
    """Create Windows batch files, to get around problem that we
    cannot run Python scripts in the prompt without the .py
    extension."""
    batch_files = []
    for script in scripts:
        batch_file = script + ".bat"
        with open(batch_file, "w") as f:
            f.write(sys.executable + " \"%%~dp0\%s\" %%*\n" % os.path.split(script)[1])
        batch_files.append(batch_file)
    scripts.extend(batch_files)
    return scripts


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


def generate_ufc_config_py_file(CXX_FLAGS, UFC_SIGNATURE):
    "Generate module with UFC signature"
    write_config_file(os.path.join("ffc", "ufc_config.py.in"),
                      os.path.join("ffc", "ufc_config.py"),
                      variables=dict(CXX_FLAGS=CXX_FLAGS,
                                     UFC_SIGNATURE=UFC_SIGNATURE))


def has_cxx_flag(cc, flag):
    "Return True if compiler supports given flag"
    tmpdir = tempfile.mkdtemp(prefix="ffc-build-")
    devnull = oldstderr = None
    try:
        try:
            fname = os.path.join(tmpdir, "flagname.cpp")
            with open(fname, "w") as f:
                f.write("int main() { return 0; }")
            # Redirect stderr to /dev/null to hide any error messages
            # from the compiler.
            devnull = open(os.devnull, 'w')
            oldstderr = os.dup(sys.stderr.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            cc.compile([fname], output_dir=tmpdir, extra_preargs=[flag])
        except:
            return False
        return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()
        shutil.rmtree(tmpdir)


def run_install():
    "Run installation"

    # Get common variables
    INSTALL_PREFIX = get_installation_prefix()
    CXX_FLAGS = get_cxx_flags()
    UFC_SIGNATURE = get_ufc_signature()
    GIT_COMMIT_HASH = get_git_commit_hash()

    # Create batch files for Windows if necessary
    scripts = SCRIPTS
    if platform.system() == "Windows" or "bdist_wininst" in sys.argv:
        scripts = create_windows_batch_files(scripts)

    # Generate module with git hash from template
    generate_git_hash_file(GIT_COMMIT_HASH)

    # Generate ufc_config.py
    generate_ufc_config_py_file(CXX_FLAGS, UFC_SIGNATURE)

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
          package_dir={"ffc": "ffc"}
          package_data={"ffc" : [os.path.join('backends', 'ufc', '*.h')]},
          scripts=scripts,
          data_files=data_files,
          install_requires=["numpy",
                            "six",
                            "fiat==%s" % VERSION,
                            "ufl==%s" % VERSION,
                            "dijitso==%s" % VERSION],
          zip_safe=False)

if __name__ == "__main__":
    run_install()
