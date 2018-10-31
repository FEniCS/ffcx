# -*- coding: utf-8 -*-
# Copyright (C) 2018 Chris N. Richardson
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import subprocess
import os
import os.path
import ffc


def test_cmdline_simple():
    os.chdir(os.path.dirname(__file__))
    subprocess.run("ffc")
    subprocess.run(["ffc", "-v", "Poisson.ufl"])
    subprocess.run(["ffc", "-l", "dolfin", "Poisson.ufl"])


def test_visualise():
    os.chdir(os.path.dirname(__file__))
    subprocess.run(["ffc", "-f", "visualise", "1", "Poisson.ufl"])
    assert os.path.isfile("S.pdf")
    assert os.path.isfile("F.pdf")


def test_cmdline_python():
    ffc.main(["-v", "Poisson.ufl"])
    ffc.main(["-f", "visualise", "1", "Poisson.ufl"])
    ffc.main(["-l", "dolfin", "Poisson.ufl"])
    ffc.main(["-v", "PoissonDG.ufl"])
