from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
import setuptools.command.build_py
import sys
import setuptools

import ffc
#mport generate_factory

__version__ = '0.0.1'


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


ext_modules = [
    Extension('ffc_test_factory.factory',
              ['build_pybind/pybind_wrapper.cpp',
               'build_pybind/elements_lagrange.cpp',
               'build_pybind/elements_vector.cpp',
               'build_pybind/elements_mixed0.cpp',
               'build_pybind/elements_mixed1.cpp',
               'build_pybind/elements_misc.cpp',
              ],
              include_dirs=[
                  # Path to pybind11 headers
                  ffc.get_include_path(),
                  get_pybind_include(),
                  get_pybind_include(user=True)
              ],
        language='c++'
    ),
]


# As of Python 3.6, CCompiler has a `has_flag` method.  cf
# http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')

class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        print("--------------------------------")
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=default'):
                opts.append('-fvisibility=default')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)


class BuildPyCommand(setuptools.command.build_py.build_py):
  """Custom build command to generate pybind11 wrappers."""

  def run(self):
      generate_factory.main()
      setuptools.command.build_py.build_py.run(self)


setup(
    name='ffc_test_factory',
    version=__version__,
    author='FEniCS Project',
    author_email='fenics-dev@googlegroups.com',
    url='https://fenicsproject.org',
    description='A factory module for creating obejcts from FFC\
    generate code for use in testing.',
    long_description='',
    packages=["ffc_test_factory"],
    package_dir={'ffc_test_factory' : 'ffc-test-factory'},
    package_data={'ffc_test_factory': ['element_factory_data.p']},
    ext_modules=ext_modules,
    install_requires=['pybind11>=1.7',
                      'ufc_wrappers'],
#    cmdclass={'build_py': BuildPyCommand,
#              'build_ext': BuildExt},
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
)
