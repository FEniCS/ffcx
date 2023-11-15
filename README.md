# FFCx: The FEniCSx Form Compiler

[![FFCx CI](https://github.com/FEniCS/ffcx/actions/workflows/pythonapp.yml/badge.svg)](https://github.com/FEniCS/ffcx/actions/workflows/pythonapp.yml)
[![Spack install](https://github.com/FEniCS/ffcx/actions/workflows/spack.yml/badge.svg)](https://github.com/FEniCS/ffcx/actions/workflows/spack.yml)
[![Coverage Status](https://coveralls.io/repos/github/FEniCS/ffcx/badge.svg?branch=main)](https://coveralls.io/github/FEniCS/ffcx?branch=main)

FFCx is a new version of the FEniCS Form Compiler. It is being actively
developed and is compatible with DOLFINx.

FFCx is a compiler for finite element variational forms. From a
high-level description of the form in the Unified Form Language (UFL),
it generates efficient low-level C code that can be used to assemble the
corresponding discrete operator (tensor). In particular, a bilinear form
may be assembled into a matrix and a linear form may be assembled into a
vector.  FFCx may be used either from the command line (by invoking the
`ffcx` command) or as a Python module (`import ffcx`).

FFCx is part of the FEniCS Project. For more information, visit
https://www.fenicsproject.org


## Installation

To install FFCx from PyPI:
```
$ pip install fenics-ffcx
```

To install FFCx from the source directory:
```
$ pip install .
```

## Documentation

Documentation can be viewed at https://docs.fenicsproject.org/ffcx/main


## Interface file installation only

FFCx provides the `ufcx.h` interface header for finite element kernels,
used by DOLFINx. `ufcx.h` is installed by FFCx within the Python site
packages, but it is sometimes helpful to install only the header file.
This can be done using `cmake`:
```
$ cmake -B build-dir -S cmake/
$ cmake --build build-dir
$ cmake --install build-dir
```

## License

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this program. If not, see <https://www.gnu.org/licenses/>.
