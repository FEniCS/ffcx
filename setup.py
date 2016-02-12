#!/usr/bin/env python

import os, sys, platform, re, subprocess, string, tempfile, shutil, hashlib

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from distutils import sysconfig, spawn
from distutils.core import Extension
from distutils.command import build_ext
from distutils.command.build import build
from distutils.ccompiler import new_compiler
from distutils.version import LooseVersion

if sys.version_info < (2, 7):
    print("Python 2.7 or higher required, please upgrade.")
    sys.exit(1)

VERSION = re.findall('__version__ = "(.*)"',
                     open('ffc/__init__.py', 'r').read())[0]

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


def get_installation_prefix():
    "Get installation prefix"
    prefix = sys.prefix
    for arg in sys.argv[1:]:
        if "--user" in arg:
            import site
            prefix = site.USER_BASE
        elif arg in ("--prefix", "--home", "--root", "--install-base"):
            prefix = sys.argv[sys.argv.index(arg)+1]
        elif "--prefix=" in arg or "--home=" in arg or \
          "--root=" in arg or "--install-base=" in arg:
            prefix = arg.split("=")[1]
    return os.path.abspath(os.path.expanduser(prefix))


def get_swig_executable():
    "Get SWIG executable"
    # Find SWIG executable
    swig_executable = None
    swig_minimum_version = "3.0.3"
    for executable in ["swig", "swig3.0"]:
        swig_executable = spawn.find_executable(executable)
        if swig_executable is not None:
            # Check that SWIG version is ok
            output = subprocess.check_output([swig_executable, "-version"]).decode('utf-8')
            swig_version = re.findall(r"SWIG Version ([0-9.]+)", output)[0]
            if LooseVersion(swig_version) >= LooseVersion(swig_minimum_version):
                break
            swig_executable = None
    if swig_executable is None:
        raise OSError("Unable to find SWIG version %s or higher." % swig_minimum_version)
    print("Found SWIG: %s (version %s)" % (swig_executable, swig_version))

    return swig_executable


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
    dirs = [
        "%s/lib" % os.environ.get("PYTHON_DIR", ""),
        "%s" % sysconfig.get_config_vars().get("LIBDIR", ""),
        "/usr/lib/%s" % sysconfig.get_config_vars().get("MULTIARCH", ""),
        "/usr/local/lib",
        "/opt/local/lib",
        "/usr/lib",
        "/usr/lib64",
        ]
    libpython = None
    cc = new_compiler()
    for name in libpython_names:
        libpython = cc.find_library_file(dirs, name)
        if libpython is not None:
            break
    return libpython or ""


def generate_config_files(SWIG_EXECUTABLE, CXX_FLAGS):
    "Generate and install configuration files"

    # Get variables
    INSTALL_PREFIX = get_installation_prefix()
    PYTHON_LIBRARY = os.environ.get("PYTHON_LIBRARY", find_python_library())
    MAJOR, MINOR, MICRO = VERSION.split(".")
    UFC_SIGNATURE = get_ufc_signature()
    GIT_COMMIT_HASH = get_git_commit_hash()

    # Generate ufc_signature.py
    write_config_file(os.path.join("ffc", "ufc_signature.py.in"),
                      os.path.join("ffc", "ufc_signature.py"),
                      variables=dict(UFC_SIGNATURE=UFC_SIGNATURE))

    # Generate git_commit_hash.py
    write_config_file(os.path.join("ffc", "git_commit_hash.py.in"),
                      os.path.join("ffc", "git_commit_hash.py"),
                      variables=dict(GIT_COMMIT_HASH=GIT_COMMIT_HASH))

    # Generate UFCConfig.cmake
    write_config_file(os.path.join("cmake", "templates", "UFCConfig.cmake.in"),
                      os.path.join("cmake", "templates", "UFCConfig.cmake"),
                      variables=dict(INSTALL_PREFIX=INSTALL_PREFIX,
                                     CXX_FLAGS=CXX_FLAGS.strip(),
                                     PYTHON_INCLUDE_DIR=sysconfig.get_python_inc(),
                                     PYTHON_LIBRARY=PYTHON_LIBRARY,
                                     PYTHON_EXECUTABLE=sys.executable,
                                     SWIG_EXECUTABLE=SWIG_EXECUTABLE,
                                     FULLVERSION=VERSION,
                                     UFC_SIGNATURE=UFC_SIGNATURE))

    # Generate UFCConfigVersion.cmake
    write_config_file(os.path.join("cmake", "templates", \
                                   "UFCConfigVersion.cmake.in"),
                      os.path.join("cmake", "templates", \
                                   "UFCConfigVersion.cmake"),
                      variables=dict(FULLVERSION=VERSION,
                                     MAJOR=MAJOR, MINOR=MINOR, MICRO=MICRO))

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


def run_ufc_install():
    # Subclass extension building command to ensure that distutils to
    # finds the correct SWIG executable
    SWIG_EXECUTABLE = get_swig_executable()
    class my_build_ext(build_ext.build_ext):
        def find_swig(self):
            return SWIG_EXECUTABLE

    # Subclass the build command to ensure that build_ext produces
    # ufc.py before build_py tries to copy it.
    class my_build(build):
        def run(self):
            self.run_command('build_ext')
            build.run(self)

    cmdclass = {"build": my_build, "build_ext": my_build_ext}

    # Check that compiler supports C++11 features
    cc = new_compiler()
    CXX = os.environ.get("CXX")
    if CXX:
        cc.set_executables(compiler_so=CXX, compiler=CXX, compiler_cxx=CXX)
    CXX_FLAGS = os.environ.get("CXXFLAGS", "")
    if has_cxx_flag(cc, "-std=c++11"):
        CXX_FLAGS += " -std=c++11"
    elif has_cxx_flag(cc, "-std=c++0x"):
        CXX_FLAGS += " -std=c++0x"

    # Generate config files
    generate_config_files(SWIG_EXECUTABLE, CXX_FLAGS)

    # Setup extension module for UFC
    swig_options = ["-c++", "-shadow", "-modern",
                    "-modernargs", "-fastdispatch",
                    "-fvirtual", "-nosafecstrings",
                    "-noproxydel", "-fastproxy",
                    "-fastinit", "-fastunpack",
                    "-fastquery", "-nobuildnone"]
    if sys.version_info[0] > 2:
        swig_options.insert(0, "-py3")
    ext_module_ufc = Extension("ufc._ufc",
                               sources=[os.path.join("ufc", "ufc.i")],
                               depends=[os.path.join("ufc", "ufc.h"),
                                        os.path.join("ufc", "ufc_geometry.h")],
                               swig_opts=swig_options,
                               extra_compile_args=CXX_FLAGS.split(),
                               include_dirs=[os.path.join("ufc")])
    ext_modules = [ext_module_ufc]
    return cmdclass, ext_modules


def run_install():
    "Run installation"

    # Check if we're building inside a 'Read the Docs' container
    on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

    # Hack to skip ufc and avoid swig dependency on install
    # so readthedocs can install without the ufc wrapper
    numpy_include_dir = None
    if "--skip-ufc" in sys.argv:
        sys.argv.remove("--skip-ufc")
        skip_ufc_module = True
    elif on_rtd:
        skip_ufc_module = True
    else:
        skip_ufc_module = False
        import numpy
        numpy_include_dir = numpy.get_include()

    # Create batch files for Windows if necessary
    scripts = SCRIPTS
    if platform.system() == "Windows" or "bdist_wininst" in sys.argv:
        scripts = create_windows_batch_files(scripts)

    if skip_ufc_module:
        cmdclass = {}
        ext_modules = []
    else:
        cmdclass, ext_modules = run_ufc_install()

    # FFC data files
    data_files = [(os.path.join("share", "man", "man1"),
                  [os.path.join("doc", "man", "man1", "ffc.1.gz")])]

    # Add UFC data files
    if not skip_ufc_module:
        data_files_ufc = [(os.path.join("include"),
                           [os.path.join("ufc", "ufc.h"),
                            os.path.join("ufc", "ufc_geometry.h")]),
                          (os.path.join("share", "ufc"),
                           [os.path.join("cmake", "templates", \
                                         "UFCConfig.cmake"),
                            os.path.join("cmake", "templates", \
                                         "UFCConfigVersion.cmake"),
                            os.path.join("cmake", "templates", \
                                         "UseUFC.cmake")]),
                          (os.path.join("lib", "pkgconfig"),
                           [os.path.join("cmake", "templates", "ufc-1.pc")]),
                          (os.path.join("include", "swig"),
                           [os.path.join("ufc", "ufc.i"),
                            os.path.join("ufc", "ufc_shared_ptr_classes.i")])]

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
          url              = "http://fenicsproject.org/",
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
                              "ufc"],
          package_dir      = {"ffc": "ffc",
                              "ufc": "ufc"},
          scripts          = scripts,
          include_dirs     = [numpy_include_dir],
          ext_modules      = ext_modules,
          cmdclass         = cmdclass,
          data_files       = data_files,
          install_requires = ["numpy", "six", "fiat", "ufl", "instant", "uflacs"])


if __name__ == "__main__":
    run_install()
