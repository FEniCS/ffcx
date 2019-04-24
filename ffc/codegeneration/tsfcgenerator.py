# Copyright (C) 2016 Jan Blechta
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging

import coffee.base as coffee
import tsfc.kernel_interface.ufc as ufc_interface
from coffee.plan import ASTKernel
from coffee.visitors import Find
from ffc.ir.representationutils import initialize_integral_code
from tsfc.driver import compile_integral

logger = logging.getLogger(__name__)


def generate_integral_code(ir, parameters):
    """Generate code for integral from intermediate representation."""

    logger.info("Generating code from tsfc representation")

    # Generate generic ffc code snippets
    code = initialize_integral_code(ir, parameters)

    # Go unoptimized if TSFC mode has not been set yet
    integral_data, form_data, prefix, parameters = ir["compile_integral"]
    parameters = parameters.copy()
    parameters.setdefault("mode", "vanilla")

    # Generate tabulate_tensor body
    ast = compile_integral(integral_data, form_data, None, parameters, interface=ufc_interface)

    # COFFEE vectorize
    knl = ASTKernel(ast)
    knl.plan_cpu(dict(optlevel='Ov'))

    tsfc_code = "".join(b.gencode() for b in ast.body)
    tsfc_code = tsfc_code.replace("#pragma coffee", "//#pragma coffee")  # FIXME
    code["tabulate_tensor"] = tsfc_code

    includes = set()
    includes.update(ir.get("additional_includes_set", ()))
    includes.update(ast.headers)
    includes.add("#include <cstring>")  # memset
    if any(
            node.funcall.symbol.startswith("boost::math::")
            for node in Find(coffee.FunCall).visit(ast)[coffee.FunCall]):
        includes.add("#include <boost/math/special_functions.hpp>")
    code["additional_includes_set"] = includes

    return code
