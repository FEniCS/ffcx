# Minimal C example

Calls an FFCx-generated kernel from C, with no other FEniCS component
involved: `poisson.py` is compiled by FFCx into `poisson.c`/`poisson.h`,
and `main.c` assembles the P1 Poisson stiffness matrix on a single
triangle and checks it against the exact result.

Building needs the UFCx interface header installed as a CMake package:

```console
cmake -S cmake -B build-ufcx -DCMAKE_INSTALL_PREFIX=<prefix>
cmake --install build-ufcx
```

from the repository root, after which:

```console
cmake -S demo/c -B build -DCMAKE_PREFIX_PATH=<prefix>
cmake --build build
./build/poisson
```

prints the element matrix and exits non-zero if it is wrong.
