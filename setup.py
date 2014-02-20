#!/usr/bin/env python

import os, sys, platform, re, subprocess, string, numpy
from distutils import sysconfig, spawn
from distutils.core import setup, Extension
from distutils.version import LooseVersion

VERSION   = "1.3.0+"
CXX_FLAGS = '-std=c++11 ' + os.environ.get('CXXFLAGS', '')

SCRIPTS = [os.path.join("scripts", "ffc")]

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
    try:
        prefix = [item for item in sys.argv[1:] \
                  if "--prefix=" in item][0].split("=")[1]
    except:
        try:
            prefix = sys.argv[sys.argv.index('--prefix')+1]
        except:
            prefix = sys.prefix
    return prefix

def get_swig_executable():
    "Get SWIG executable"

    # Find SWIG executable
    swig_executable = None
    for executable in ['swig', 'swig2.0']:
        swig_executable = spawn.find_executable(executable)
        if swig_executable is not None:
            break
    if swig_executable is None:
        raise OSError('Unable to find SWIG installation. Please install SWIG version 2.0 or higher.')

    # Check that SWIG version is ok
    output = subprocess.check_output([swig_executable, "-version"])
    swig_version = re.findall(r'SWIG Version ([0-9.]+)', output)[0]
    swig_version_ok = True
    swig_minimum_version = [2, 0, 0]
    for i, v in enumerate([int(v) for v in swig_version.split('.')]):
        if swig_minimum_version[i] < v:
            break
        elif swig_minimum_version[i] == v:
            continue
        else:
            swig_version_ok = False
    if not swig_version_ok:
        raise OSError('Unable to find SWIG version 2.0 or higher.')

    return swig_executable

def create_windows_batch_files(scripts):
    """Create Windows batch files, to get around problem that we
    cannot run Python scripts in the prompt without the .py
    extension."""
    batch_files = []
    for script in scripts:
        batch_file = script + ".bat"
        f = open(batch_file, "w")
        f.write('python "%%~dp0\%s" %%*\n' % os.path.split(script)[1])
        f.close()
        batch_files.append(batch_file)
    scripts.extend(batch_files)
    return scripts

def write_config_file(infile, outfile, variables={}):
    "Write config file based on template"
    class AtTemplate(string.Template):
        delimiter = '@'
    s = AtTemplate(open(infile, 'r').read())
    s = s.substitute(**variables)
    a = open(outfile, 'w')
    try:
        a.write(s)
    finally:
        a.close()

def generate_config_files():
    "Generate and install configuration files"

    # Get variables
    INSTALL_PREFIX  = get_installation_prefix()
    SWIG_EXECUTABLE = get_swig_executable()
    PYTHON_LIBRARY  = 'libpython2.7.so'
    MAJOR, MINOR, MICRO = VERSION.split(".")

    # Generate UFCConfig.cmake
    write_config_file(os.path.join('cmake/templates', 'UFCConfig.cmake.in'),
                      os.path.join('cmake/templates', 'UFCConfig.cmake'),
                      variables=dict(INSTALL_PREFIX=INSTALL_PREFIX,
                                     CXX_FLAGS=CXX_FLAGS,
                                     PYTHON_INCLUDE_DIR=sysconfig.get_python_inc(),
                                     PYTHON_LIBRARY=PYTHON_LIBRARY,
                                     PYTHON_EXECUTABLE=sys.executable,
                                     SWIG_EXECUTABLE=SWIG_EXECUTABLE,
                                     FULLVERSION=VERSION))

    # Generate UFCConfigVersion.cmake
    write_config_file(os.path.join('cmake/templates', 'UFCConfigVersion.cmake.in'),
                      os.path.join('cmake/templates', 'UFCConfigVersion.cmake'),
                      variables=dict(FULLVERSION=VERSION,
                                     MAJOR=MAJOR, MINOR=MINOR, MICRO=MICRO))

    # Generate UseUFC.cmake
    write_config_file(os.path.join('cmake/templates', 'UseUFC.cmake.in'),
                      os.path.join('cmake/templates', 'UseUFC.cmake'))

    # FIXME: Generation of pkgconfig file may no longer be needed, so
    # FIXME: we may consider removing this.

    # Generate ufc-1.pc
    write_config_file(os.path.join('cmake/templates', 'ufc-1.pc.in'),
                      os.path.join('cmake/templates', 'ufc-1.pc'),
                      variables=dict(FULLVERSION=VERSION,
                                     INSTALL_PREFIX=INSTALL_PREFIX,
                                     CXX_FLAGS=CXX_FLAGS))

def run_install():
    "Run installation"

    # Create batch files for Windows if necessary
    scripts = SCRIPTS
    if platform.system() == "Windows" or "bdist_wininst" in sys.argv:
        scripts = create_windows_batch_files(scripts)

    # Generate config files
    generate_config_files()

    # FIXME: Someone needs to explain what this does. Benjamin?
    ext_kwargs = dict(include_dirs=[numpy.get_include()])
    if LooseVersion(numpy.__version__) > LooseVersion("1.6.2"):
        ext_kwargs["define_macros"] = [("NPY_NO_DEPRECATED_API", "NPY_%s_%s_API_VERSION" \
                                            % tuple(numpy.__version__.split(".")[:2]))]

    # Setup extension module for FFC time elements
    ext_module_time = Extension("ffc.time_elements_ext",
                                ["ffc/ext/time_elements_interface.cpp",
                                 "ffc/ext/time_elements.cpp",
                                 "ffc/ext/LobattoQuadrature.cpp",
                                 "ffc/ext/RadauQuadrature.cpp",
                                 "ffc/ext/Legendre.cpp"],
                                **ext_kwargs)

    # Setup extension module for UFC
    ext_module_ufc = [Extension("_ufc",
                                sources=[os.path.join('src', 'ufc', 'ufc.i')],
                                swig_opts=['-c++', '-shadow', '-modern',
                                           '-modernargs', '-fastdispatch',
                                           '-fvirtual', '-nosafecstrings',
                                           '-noproxydel', '-fastproxy',
                                           '-fastinit', '-fastunpack',
                                           '-fastquery', '-nobuildnone'],
                                extra_compile_args=CXX_FLAGS.split(),
                                include_dirs=[os.path.join('src', 'ufc')])]

    # Call distutils to perform installation
    setup(name         = "FFC",
          version      = VERSION,
          description  = "The FEniCS Form Compiler",
          classifiers  = CLASSIFIERS,
          author       = "Anders Logg, Kristian Oelgaard, Marie Rognes, Garth N. Wells,  et al.",
          author_email = "fenics@fenicsproject.org",
          url          = "http://www.fenicsproject.org",
          packages     = ["ffc",
                          "ffc.quadrature", "ffc.tensor", "ffc.uflacsrepr",
                          "ffc.errorcontrol",
                          "ffc.dolfin"],
          package_dir  = {"ffc": "ffc"},
          scripts      = scripts,
          ext_modules  = [ext_module_time],
          data_files   = [(os.path.join("share", "man", "man1"),
                           [os.path.join("doc", "man", "man1", "ffc.1.gz")])])

def setup_package_ufc():

    # FIXME: Contents of this function should be moved:
    # - Merge call to setup with FFC call to setup below

    class my_build_ext(build_ext.build_ext):
        def find_swig(self):
            return SWIG_EXECUTABLE

    # Run setup
    setup(name='UFC',
          version=VERSION,
          maintainer='FEniCS Developers',
          maintainer_email='fenics@fenicsproject.org',
          description=DOCLINES[0],
          url='http://fenicsproject.org',
          author='Martin Sandve Alnaes, Hans Petter Langtangen, Anders Logg, Kent-Andre Mardal, Ola Skavhaug, et al.',
          download_url='https://bitbucket.org/fenics-project/ufc/downloads',
          license='Public Domain',
          platforms=['Windows', 'Linux', 'Solaris', 'Mac OS-X', 'Unix'],
          packages=['ufc', 'ufc_utils'],
          package_dir={'ufc': os.path.join('src', 'ufc'),
                       'ufc_utils': os.path.join('src', 'utils', 'python', 'ufc_utils')},
          data_files=[('include', [os.path.join('src', 'ufc', 'ufc.h'),
                                   os.path.join('src', 'ufc', 'ufc_geometry.h')]),
                      (os.path.join('share', 'ufc'),
                       [os.path.join('templates', 'UFCConfig.cmake'),
                        os.path.join('templates', 'UFCConfigVersion.cmake'),
                        os.path.join('templates', 'UseUFC.cmake')]),
                      (os.path.join('lib', 'pkgconfig'),
                       [os.path.join('templates', 'ufc-1.pc')]),
                      (os.path.join('include', 'swig'),
                       [os.path.join('src', 'ufc', 'ufc.i')])],
          ext_package="ufc",
          ext_modules=modules,
          cmdclass={'build_ext': my_build_ext},
          )

if __name__ == '__main__':
    run_install()
