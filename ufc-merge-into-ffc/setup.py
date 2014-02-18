#!/usr/bin/env python

"""UFC: unified code generation interface for form-compilers

UFC (Unified Form-assembly Code) is a unified framework for finite element
assembly. More precisely, it defines a fixed interface for communicating low
level routines (functions) for evaluating and assembling finite element
variational forms. The UFC interface consists of a single header file ufc.h
that specifies a C++ interface that must be implemented by code that complies
with the UFC specification. Examples of form compilers that support the UFC
interface are FFC and SyFi.
"""

from __future__ import division, print_function

import os
import re
import sys
import string
import subprocess

from distutils.core import setup, Extension
from distutils import sysconfig, spawn

MAJOR               = 2
MINOR               = 3
MICRO               = 0
ISRELEASED          = False
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)

DOCLINES = __doc__.split("\n")

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

def get_status_output(cmd, input=None, cwd=None, env=None):
    "Run command and return status and output"
    pipe = subprocess.Popen(cmd, shell=True, cwd=cwd, env=env,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    output, errout = pipe.communicate(input=input)
    assert not errout
    status = pipe.returncode
    return status, output

def git_version():
    "Return the git revision as a string"
    GIT_REVISION = "Unknown"
    failure, output = get_status_output('git rev-parse HEAD')
    if not failure:
        GIT_REVISION = output.strip().decode('ascii')

    return GIT_REVISION

def get_version_info():
    "Get version number and Git revision"
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        #FULLVERSION += '.dev-' + GIT_REVISION[:7]
        FULLVERSION += '+'

    return FULLVERSION, GIT_REVISION

def get_prefix():
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


    # Find SWIG executable
    swig_executable = None
    for executable in ['swig', 'swig2.0']:
        swig_executable = spawn.find_executable(executable)
        if swig_executable is not None:
            break
    if swig_executable is None:
        raise OSError('Unable to find SWIG installation. Please install SWIG version 2.0 or higher.')

    # Check that SWIG version is ok
    _, output = get_status_output('%s -version' % swig_executable)

    swig_version = re.findall(r'SWIG Version ([0-9.]+)', output)[0]
    swig_version_ok = True
    swig_minimum_version = [2,0,0]
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

def configure_file(infile, outfile, variables={}):
    class AtTemplate(string.Template):
        delimiter = '@'
    s = AtTemplate(open(infile, 'r').read())
    s = s.substitute(**variables)
    a = open(outfile, 'w')
    try:
        a.write(s)
    finally:
        a.close()

def setup_package():
    PREFIX = get_prefix()
    # FIXME: How can we get the path to the Python library? Do we need it?
    PYTHON_LIBRARY = 'libpython2.7.so'
    FULLVERSION, GIT_REVISION = get_version_info()
    CXX_FLAGS = '-std=c++11 ' + os.environ.get('CXXFLAGS', '')
    SWIG_EXECUTABLE = get_swig_executable()

    # Rewrite the CMake files every time
    configure_file(os.path.join('templates', 'UFCConfig.cmake.in'),
                   os.path.join('templates', 'UFCConfig.cmake'),
                   variables=dict(INSTALL_PREFIX=PREFIX,
                                  CXX_FLAGS=CXX_FLAGS,
                                  PYTHON_INCLUDE_DIR=sysconfig.get_python_inc(),
                                  PYTHON_LIBRARY=PYTHON_LIBRARY,
                                  PYTHON_EXECUTABLE=sys.executable,
                                  SWIG_EXECUTABLE=SWIG_EXECUTABLE,
                                  FULLVERSION=FULLVERSION))
    configure_file(os.path.join('templates', 'UFCConfigVersion.cmake.in'),
                   os.path.join('templates', 'UFCConfigVersion.cmake'),
                   variables=dict(FULLVERSION=FULLVERSION,
                                  MAJOR=MAJOR, MINOR=MINOR, MICRO=MICRO))
    configure_file(os.path.join('templates', 'UseUFC.cmake.in'),
                   os.path.join('templates', 'UseUFC.cmake'))
    # FIXME: Do we still need to generate the pkg-config file?
    configure_file(os.path.join('templates', 'ufc-1.pc.in'),
                   os.path.join('templates', 'ufc-1.pc'),
                   variables=dict(FULLVERSION=FULLVERSION,
                                  INSTALL_PREFIX=PREFIX,
                                  CXX_FLAGS=CXX_FLAGS))

    # Subclass build_ext to help distutils find the correct SWIG executable
    from distutils.command import build_ext
    class my_build_ext(build_ext.build_ext):
        def find_swig(self):
            return SWIG_EXECUTABLE

    modules = [Extension("_ufc",
                         sources=[os.path.join('src', 'ufc', 'ufc.i')],
                         swig_opts=['-c++', '-shadow', '-modern',
                                    '-modernargs', '-fastdispatch',
                                    '-fvirtual', '-nosafecstrings',
                                    '-noproxydel', '-fastproxy',
                                    '-fastinit', '-fastunpack',
                                    '-fastquery', '-nobuildnone'],
                         extra_compile_args=CXX_FLAGS.split(),
                         include_dirs=[os.path.join('src', 'ufc')])]

    setup(name='UFC',
          version=FULLVERSION,
          maintainer='FEniCS Developers',
          maintainer_email='fenics@fenicsproject.org',
          description=DOCLINES[0],
          long_description='\n'.join(DOCLINES[2:]),
          url='http://fenicsproject.org',
          author='Martin Sandve Alnaes, Hans Petter Langtangen, Anders Logg, Kent-Andre Mardal, Ola Skavhaug, et al.',
          download_url='https://bitbucket.org/fenics-project/ufc/downloads',
          license='Public Domain',
          classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
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
    setup_package()
