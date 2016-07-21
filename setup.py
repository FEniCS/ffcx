#!/usr/bin/env python

import os, sys, platform, re, subprocess, string, tempfile, shutil, hashlib

try:
    from setuptools import setup
    from setuptools.command.install import install
except ImportError:
    from distutils.core import setup
    from distutils.command.install import install

from distutils import sysconfig
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
License :: OSI Approved :: GNU General Public License v2 (GPLv2)
License :: Public Domain
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: POSIX :: Linux
Programming Language :: C++
Programming Language :: Python
Topic :: Scientific/Engineering :: Mathematics
Topic :: Software Development :: Libraries
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
    with open(os.path.join('ufc', 'ufc.h'), 'rb') as f:
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
        f = open(batch_file, "w")
        f.write("python \"%%~dp0\%s\" %%*\n" % os.path.split(script)[1])
        f.close()
        batch_files.append(batch_file)
    scripts.extend(batch_files)
    return scripts


def write_config_file(infile, outfile, variables={}):
    "Write config file based on template"
    class AtTemplate(string.Template):
        delimiter = "@"
    s = AtTemplate(open(infile, "r").read())
    s = s.substitute(**variables)
    a = open(outfile, "w")
    try:
        a.write(s)
    finally:
        a.close()


def find_library(package_name, lib_names):
    "Return the full path to the library (empty string if not found)"
    search_dirs = [
        "%s%slib" % (os.environ.get("%s_DIR" % package_name.upper(), ""),
                     os.path.sep),
        "%s" % sysconfig.get_config_vars().get("LIBDIR", ""),
        "/usr/lib/%s" % sysconfig.get_config_vars().get("MULTIARCH", ""),
        "/usr/local/lib",
        "/opt/local/lib",
        "/usr/lib",
        "/usr/lib64",
        ]
    lib = None
    cc = new_compiler()
    for name in lib_names:
        lib = cc.find_library_file(search_dirs, name)
        if lib is not None:
            break
    return lib or ""


def find_python_library():
    "Return the full path to the Python library (empty string if not found)"
    pyver = sysconfig.get_python_version()
    libpython_names = [
        "python%s" % pyver.replace(".", ""),
        "python%smu" % pyver,
        "python%sm" % pyver,
        "python%su" % pyver,
        "python%s" % pyver,
        ]
    return find_library("python", libpython_names)


def find_boost_math_library():
    "Return the full path to the Boost math library (empty string if not found)"
    return find_library("boost", ["boost_math_tr1", "boost_math_tr1-mt"])


def find_include_dir(package_name, header_file):
    "Return the path to the given header file (empty string if not found)"
    search_dirs = [
        "%s%sinclude" % (os.environ.get("%s_DIR" % package_name.upper(), ""),
                         os.path.sep),
        "/usr/local/include",
        "/opt/local/include",
        "/usr/include",
        ]
    for inc_dir in search_dirs:
        if os.path.isfile(os.path.join(inc_dir, header_file)):
            return inc_dir
    return ""


def find_boost_include_dir():
    "Return the path to the Boost include dir (empty string if not found)"
    return find_include_dir("boost", os.path.join("boost", "version.hpp"))


def generate_git_hash_file(GIT_COMMIT_HASH):
    "Generate module with git hash"
    write_config_file(os.path.join("ffc", "git_commit_hash.py.in"),
                      os.path.join("ffc", "git_commit_hash.py"),
                      variables=dict(GIT_COMMIT_HASH=GIT_COMMIT_HASH))


def generate_ufc_config_py_file(INSTALL_PREFIX, CXX_FLAGS, UFC_SIGNATURE):
    "Generate module with UFC signature"
    write_config_file(os.path.join("ffc", "ufc_config.py.in"),
                      os.path.join("ffc", "ufc_config.py"),
                      variables=dict(INSTALL_PREFIX=INSTALL_PREFIX,
                                     CXX_FLAGS=CXX_FLAGS,
                                     UFC_SIGNATURE=UFC_SIGNATURE))



def generate_ufc_config_files(INSTALL_PREFIX, CXX_FLAGS, UFC_SIGNATURE):
    "Generate and install UFC configuration files"

    # Get variables
    PYTHON_LIBRARY = os.environ.get("PYTHON_LIBRARY", find_python_library())

    # Generate UFCConfig.cmake
    write_config_file(os.path.join("cmake", "templates", "UFCConfig.cmake.in"),
                      os.path.join("cmake", "templates", "UFCConfig.cmake"),
                      variables=dict(INSTALL_PREFIX=INSTALL_PREFIX,
                                     CXX_FLAGS=CXX_FLAGS.strip(),
                                     PYTHON_INCLUDE_DIR=sysconfig.get_python_inc(),
                                     PYTHON_LIBRARY=PYTHON_LIBRARY,
                                     PYTHON_EXECUTABLE=sys.executable,
                                     FULLVERSION=VERSION,
                                     UFC_SIGNATURE=UFC_SIGNATURE,
                                     BOOST_INCLUDE_DIR=find_boost_include_dir(),
                                     BOOST_MATH_LIBRARY=find_boost_math_library()))

    # Generate UFCConfigVersion.cmake
    write_config_file(os.path.join("cmake", "templates",
                                   "UFCConfigVersion.cmake.in"),
                      os.path.join("cmake", "templates",
                                   "UFCConfigVersion.cmake"),
                      variables=dict(FULLVERSION=VERSION,
                                     MAJOR=VERSION.split(".")[0],
                                     MINOR=VERSION.split(".")[1],
                                     MICRO=VERSION.split(".")[2]))

    # Generate UseUFC.cmake
    write_config_file(os.path.join("cmake", "templates", "UseUFC.cmake.in"),
                      os.path.join("cmake", "templates", "UseUFC.cmake"))

    # FIXME: Generation of pkgconfig file may no longer be needed, so
    # FIXME: we may consider removing this.

    # Generate ufc-1.pc
    write_config_file(os.path.join("cmake", "templates", "ufc-1.pc.in"),
                      os.path.join("cmake", "templates", "ufc-1.pc"),
                      variables=dict(FULLVERSION=VERSION,
                                     INSTALL_PREFIX=INSTALL_PREFIX,
                                     CXX_FLAGS=CXX_FLAGS))


def has_cxx_flag(cc, flag):
    "Return True if compiler supports given flag"
    tmpdir = tempfile.mkdtemp(prefix="ffc-build-")
    devnull = oldstderr = None
    try:
        try:
            fname = os.path.join(tmpdir, "flagname.cpp")
            f = open(fname, "w")
            f.write("int main() { return 0;}")
            f.close()
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

    # Check if we're building inside a 'Read the Docs' container
    on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

    # Create batch files for Windows if necessary
    scripts = SCRIPTS
    if platform.system() == "Windows" or "bdist_wininst" in sys.argv:
        scripts = create_windows_batch_files(scripts)

    # Generate module with git hash from template
    generate_git_hash_file(GIT_COMMIT_HASH)

    # Generate config files
    generate_ufc_config_files(INSTALL_PREFIX, CXX_FLAGS, UFC_SIGNATURE)

    class my_install(install):
        def run(self):
            if not self.dry_run:
                # Generate ufc_config.py
                generate_ufc_config_py_file(INSTALL_PREFIX, CXX_FLAGS,
                                            UFC_SIGNATURE)

            # distutils uses old-style classes, so no super()
            install.run(self)

    # FFC data files
    data_files = [(os.path.join("share", "man", "man1"),
                  [os.path.join("doc", "man", "man1", "ffc.1.gz")])]

    # Add UFC data files (need to use complete path because setuptools
    # installs into the Python package directory, not --prefix). This
    # can be fixed when Swig, etc are removed from FFC).
    data_files_ufc = [(os.path.join(INSTALL_PREFIX, "include"),
                       [os.path.join("ufc", "ufc.h"),
                        os.path.join("ufc", "ufc_geometry.h")]),
                      (os.path.join(INSTALL_PREFIX, "share", "ufc"),
                       [os.path.join("cmake", "templates",
                                     "UFCConfig.cmake"),
                        os.path.join("cmake", "templates",
                                     "UFCConfigVersion.cmake"),
                        os.path.join("cmake", "templates",
                                     "UseUFC.cmake")]),
                      (os.path.join(INSTALL_PREFIX, "lib", "pkgconfig"),
                       [os.path.join("cmake", "templates", "ufc-1.pc")])]

    data_files = data_files + data_files_ufc

    # Call distutils to perform installation
    setup(name             = "FFC",
          description      = "The FEniCS Form Compiler",
          version          = VERSION,
          author           = AUTHORS,
          classifiers      = [_f for _f in CLASSIFIERS.split('\n') if _f],
          license          = "LGPL version 3 or later",
          author_email     = "fenics-dev@googlegroups.com",
          maintainer_email = "fenics-dev@googlegroups.com",
          url              = URL,
          download_url     = tarball(),
          platforms        = ["Windows", "Linux", "Solaris", "Mac OS-X",
                              "Unix"],
          packages         = ["ffc",
                              "ffc.quadrature",
                              "ffc.tensor",
                              "ffc.uflacsrepr",
                              "ffc.errorcontrol",
                              "ffc.backends",
                              "ffc.backends.dolfin",
                              "ffc.backends.ufc",
                              "uflacs",
                              "uflacs.analysis",
                              "uflacs.backends",
                              "uflacs.backends.ffc",
                              "uflacs.backends.ufc",
                              "uflacs.datastructures",
                              "uflacs.elementtables",
                              "uflacs.generation",
                              "uflacs.language",
                              "uflacs.representation",
                              "ufc"],
          package_dir      = {"ffc": "ffc",
                              "uflacs": "uflacs",
                              "ufc": "ufc"},
          scripts          = scripts,
          cmdclass         = {'install': my_install},
          data_files       = data_files,
          install_requires = ["numpy",
                              "six",
                              "fiat==%s" % VERSION,
                              "ufl==%s" % VERSION,
                              "dijitso==%s" % VERSION],
          zip_safe = False)

if __name__ == "__main__":
    run_install()
