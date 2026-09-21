# Copyright (C) 2019 Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import importlib._bootstrap_external
import os
import sys
from pathlib import Path
from unittest import mock

import basix.ufl
import pytest
import ufl

import ffcx.codegeneration.jit


def test_cache_modes(compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]

    # Load form from /tmp
    _compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )
    tmpname = module.__name__
    tmpfile = module.__file__
    print(tmpname, tmpfile)
    del sys.modules[tmpname]

    # Load form from cache
    _compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, cache_dir="./compile-cache", cffi_extra_compile_args=compile_args
    )
    newname = module.__name__
    newfile = module.__file__
    print(newname, newfile)

    assert newname == tmpname
    assert newfile != tmpfile


def test_cache_hit_does_not_scan_cache_dir(compile_args, tmp_path):
    """A cache hit must not list the cache directory.

    The finder previously used listed the whole directory on every hit,
    and the cache grows by four files per form ever compiled.

    ``importlib._bootstrap_external`` binds ``listdir`` at import time,
    so intercept it there rather than on ``os``.
    """
    bootstrap = importlib._bootstrap_external
    if not hasattr(bootstrap, "_os") or not hasattr(bootstrap._os, "listdir"):
        pytest.skip("cannot intercept the import system's directory listing")

    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    forms = [ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx]

    cache_dir = tmp_path / "cache"
    _, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cache_dir=cache_dir, cffi_extra_compile_args=compile_args
    )
    del sys.modules[module.__name__]

    listed = []
    real_listdir = bootstrap._os.listdir

    def counting_listdir(path="."):
        listed.append(os.fspath(path))
        return real_listdir(path)

    with mock.patch.object(bootstrap._os, "listdir", counting_listdir):
        _, cached, _ = ffcx.codegeneration.jit.compile_forms(
            forms, cache_dir=cache_dir, cffi_extra_compile_args=compile_args
        )

    assert cached.__name__ == module.__name__
    scanned = [p for p in listed if Path(p).resolve() == cache_dir.resolve()]
    assert not scanned, f"cache hit listed the cache directory {len(scanned)} time(s)"
