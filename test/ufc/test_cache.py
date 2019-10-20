# Copyright (C) 2019 Chris Richardson
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import sys
import ffc.codegeneration.jit
import ufl


def test_cache_modes():
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]

    # Load form from /tmp
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(forms)
    tmpname = module.__name__
    tmpfile = module.__file__
    print(tmpname, tmpfile)
    del sys.modules[tmpname]

    # Load form from cache
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(forms, cache_dir="./compile-cache")
    newname = module.__name__
    newfile = module.__file__
    print(newname, newfile)

    assert(newname == tmpname)
    assert(newfile != tmpfile)
