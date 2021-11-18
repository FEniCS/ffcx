# Copyright (C) 2021 Igor A. Baratta
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import ffcx.parameters
import ufl
from ffcx.analysis import analyze_ufl_objects
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.integrals import IntegralGenerator
from ffcx.ir.representation import compute_ir
from typing import Optional


def count_flops(form: ufl.Form, parameters: Optional[dict] = {}):
    """Return a list with the estimated number of floating point operations
    for each kernel in the ufl.Form.
    """
    parameters = ffcx.parameters.get_parameters(parameters)
    assert(isinstance(form, ufl.Form))
    analysis = analyze_ufl_objects([form], parameters)
    ir = compute_ir(analysis, {}, "flops", parameters, False)

    flops = []

    for integral_ir in ir.integrals:
        # Create FFCx C backend
        backend = FFCXBackend(integral_ir, parameters)
        # Configure kernel generator
        ig = IntegralGenerator(integral_ir, backend)
        # Generate code ast for the tabulate_tensor body
        ast = ig.generate()
        flops.append(ast.flops())

    return flops
