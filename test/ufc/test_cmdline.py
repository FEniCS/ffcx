# -*- coding: utf-8 -*-
# Copyright (C) 2018 Chris N. Richardson
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import subprocess
import os.path


def test_cmdline_simple():
    subprocess.run("ffc")
    subprocess.run(["ffc", "-v", "Poisson.ufl"])
    subprocess.run(["ffc", "-l", "dolfin", "Poisson.ufl"])


def test_visualise():
    subprocess.run(["ffc", "-f", "visualise", "1", "Poisson.ufl"])
    assert os.path.isfile("S.pdf")
    assert os.path.isfile("F.pdf")
