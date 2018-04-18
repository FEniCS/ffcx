# Copyright (C) 2016 Jan Blechta
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later


def optimize_integral_ir(ir, parameters):
    ir = ir.copy()
    integral_data, form_data, prefix, parameters = ir["compile_integral"]
    parameters = parameters.copy()
    parameters.setdefault("mode", "spectral")  # default optimization mode
    ir["compile_integral"] = (integral_data, form_data, prefix, parameters)
    return ir
