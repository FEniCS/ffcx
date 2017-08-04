# Copyright (C) 2016 Jan Blechta
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.

import coffee.base as coffee
from coffee.plan import ASTKernel
from coffee.visitors import Find

from ffc.log import info
from ffc.representationutils import initialize_integral_code

from tsfc.driver import compile_integral
import tsfc.kernel_interface.ufc as ufc_interface


def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."

    info("Generating code from tsfc representation")

    # Generate generic ffc code snippets
    code = initialize_integral_code(ir, prefix, parameters)

    # Go unoptimized if TSFC mode has not been set yet
    integral_data, form_data, prefix, parameters = ir["compile_integral"]
    parameters = parameters.copy()
    parameters.setdefault("mode", "vanilla")

    # Generate tabulate_tensor body
    ast = compile_integral(integral_data, form_data, None, parameters,
                           interface=ufc_interface)

    # COFFEE vectorize
    knl = ASTKernel(ast)
    knl.plan_cpu(dict(optlevel='Ov'))

    tsfc_code = "".join(b.gencode() for b in ast.body)
    tsfc_code = tsfc_code.replace("#pragma coffee", "//#pragma coffee") # FIXME
    code["tabulate_tensor"] = tsfc_code

    includes = set()
    includes.update(ir.get("additional_includes_set", ()))
    includes.update(ast.headers)
    includes.add("#include <cstring>")  # memset
    if any(node.funcall.symbol.startswith("boost::math::")
           for node in Find(coffee.FunCall).visit(ast)[coffee.FunCall]):
        includes.add("#include <boost/math/special_functions.hpp>")
    code["additional_includes_set"] = includes

    return code
