# Copyright (C) 2020 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Test configuration."""

import os
import sys

import pytest


@pytest.fixture(autouse=True, scope="session")
def add_cwd_to_syspath():
    """Ensure the current working directory is in sys.path.

    CFFI and ffcx.main compile/generate files into CWD by default. Without
    this, importlib.import_module cannot find those files when ``pytest``
    is invoked from a directory that is not already on sys.path
    (e.g. the project root rather than the test/ subdirectory).
    """
    cwd = os.getcwd()
    if cwd not in sys.path:
        sys.path.append(cwd)


@pytest.fixture(scope="module")
def compile_args():
    """Compiler arguments."""
    if sys.platform.startswith("win32"):
        return ["-Od"]
    else:
        return ["-O1", "-Wall", "-Werror"]
