# Copyright (C) 2018 Chris N. Richardson
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import os
import os.path

import ffc


def test_forms():
    os.chdir(os.path.dirname(__file__))
    ffc.main(["-v", "Poisson.ufl"])
    ffc.main(["--visualise", "1", "Poisson.ufl"])
    ffc.main(["-v", "PoissonDG.ufl"])
    ffc.main(["-v", "Conditional.ufl"])
    ffc.main(["-v", "HyperElasticity.ufl"])
    ffc.main(["-v", "VectorLaplaceGradCurl.ufl"])
    ffc.main(["-v", "ProjectionManifold.ufl"])
    ffc.main(["-v", "Symmetry.ufl"])
