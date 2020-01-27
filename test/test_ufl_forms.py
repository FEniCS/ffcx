# Copyright (C) 2018 Chris N. Richardson
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import os
import os.path

import ffcx


def test_forms():
    os.chdir(os.path.join(os.path.dirname(__file__), "ufl_forms"))
    ffcx.main(["-v", "Poisson.ufl"])
    ffcx.main(["--visualise", "1", "Poisson.ufl"])
    ffcx.main(["-v", "PoissonDG.ufl"])
    ffcx.main(["-v", "Conditional.ufl"])
    ffcx.main(["-v", "HyperElasticity.ufl"])
    ffcx.main(["-v", "VectorLaplaceGradCurl.ufl"])
    ffcx.main(["-v", "ProjectionManifold.ufl"])
    ffcx.main(["-v", "Symmetry.ufl"])
    ffcx.main(["-v", "MixedGradient.ufl"])
