# -*- coding: utf-8 -*-
# Copyright (C) 2018 Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
#from cffi import FFI

import ffc
import ffc.backends.ufc.jit
import ufl

# Build list of UFL elements
cell = ufl.triangle
elements = [ufl.FiniteElement("Lagrange", cell, p) for p in range(1, 3)]

# Compile elements
compiled_elements, module = ffc.backends.ufc.jit.compile_elements(elements)

# Test
for e, compiled_e in zip(elements, compiled_elements):
    print(compiled_e.geometric_dimension)
    print(compiled_e.degree)
    print(compiled_e)
    print(compiled_e.family)
    test = module.ffi.string(compiled_e.family)

    tdim = compiled_e.topological_dimension
    space_dim = compiled_e.space_dimension
    X = np.zeros([space_dim, tdim])

    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    compiled_e.tabulate_reference_dof_coordinates(X_ptr)
    print(X)
