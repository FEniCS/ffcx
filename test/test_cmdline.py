# Copyright (C) 2018-2025 Chris N. Richardson and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import importlib
import subprocess
from pathlib import Path

import pytest


def test_cmdline_simple():
    dir = Path(__file__).parent.resolve()
    subprocess.run(["ffcx", dir / "poisson.py", "-o", dir], check=True)


@pytest.mark.skipif(
    importlib.util.find_spec("pygraphviz") is None, reason="pygraphviz not installed."
)
def test_visualise():
    dir = Path(__file__).parent.resolve()
    subprocess.run(["ffcx", "--visualise", dir / "poisson.py", "-o", dir])
    assert (dir / "S.pdf").is_file()
    assert (dir / "F.pdf").is_file()
