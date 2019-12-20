# Copyright (C) 2018 Chris N. Richardson
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import os
import os.path
import subprocess


def test_cmdline_simple():
    os.chdir(os.path.join("ufl_forms", os.path.dirname(__file__)))
    subprocess.run("ffc")
    subprocess.run(["ffc", "-v", "Poisson.ufl"])
    subprocess.run(["ffc", "Poisson.ufl"])


def test_visualise():
    os.chdir(os.path.join("ufl_forms", os.path.dirname(__file__)))
    subprocess.run(["ffc", "--visualise", "Poisson.ufl"])
    assert os.path.isfile("S.pdf")
    assert os.path.isfile("F.pdf")
