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

from ffc.log import info
from ffc.representationutils import initialize_integral_code


def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."

    info("Generating code from tsfc representation")

    # Generate generic ffc code snippets
    code = initialize_integral_code(ir, prefix, parameters)

    # Generate tabulate_tensor body from ast
    ast = ir["tsfc"].ast
    tsfc_code = "".join(b.gencode() for b in ast.body)
    tsfc_code = tsfc_code.replace("#pragma coffee", "//#pragma coffee") # FIXME
    code["tabulate_tensor"] = tsfc_code

    code["additional_includes_set"] = set()
    code["additional_includes_set"].update(ir.get("additional_includes_set",()))
    code["additional_includes_set"].update(ast.headers)
    code["additional_includes_set"].add("#include <cstring>") # memset

    return code
