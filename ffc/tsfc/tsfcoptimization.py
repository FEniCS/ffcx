# Copyright (C) 2016 Jan Blechta
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later


def optimize_integral_ir(ir, parameters):
    ir = ir.copy()
    integral_data, form_data, prefix, parameters = ir["compile_integral"]
    parameters = parameters.copy()

    # Spectral mode not currently working with complex scalars
    # Related to the issue https://github.com/firedrakeproject/tsfc/issues/166
    if parameters.get("scalar_type", "double") == "double":
        parameters.setdefault("mode", "spectral")  # default optimization mode
    else:
        # FIXME: What sensible optimized mode should used? Keep vanilla?
        pass
    ir["compile_integral"] = (integral_data, form_data, prefix, parameters)
    return ir
