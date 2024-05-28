# Copyright (C) 2020 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Test configuration."""

import sys

import pytest


@pytest.fixture(scope="module")
def compile_args():
    """Compiler arguments."""
    if sys.platform.startswith("win32"):
        return ["-Od"]
    else:
        return ["-O1", "-Wall", "-Werror"]
