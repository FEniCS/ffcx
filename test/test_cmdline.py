# Copyright (C) 2018 Chris N. Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import os
import os.path
import subprocess

import pytest


def test_cmdline_simple():
    os.chdir(os.path.dirname(__file__))
    subprocess.run(["ffcx", "Poisson.py"])


def test_visualise():
    try:
        import pygraphviz  # noqa: F401
    except ImportError:
        pytest.skip("pygraphviz not installed")

    os.chdir(os.path.dirname(__file__))
    subprocess.run(["ffcx", "--visualise", "Poisson.py"])
    assert os.path.isfile("S.pdf")
    assert os.path.isfile("F.pdf")
