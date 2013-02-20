# Copyright (C) 2013 Martin Alnaes
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
#
# First added:  2013-02-12
# Last changed: 2013-02-14

from ffc.representationutils import initialize_integral_code
from ffc.log import info, error, begin, end, debug_ir, ffc_assert, warning
from ffc.cpp import format

from ffc.quadrature.quadraturegenerator import _tabulate_psis

def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."

    info("Generating code from uflacs representation")

    # Generate generic ffc code snippets
    code = initialize_integral_code(ir, prefix, parameters)

    # Generate code for basis function tables
    used_psi_tables = set(ir["unique_tables"].keys()) # TODO: Build this from required set in uflacs?
    used_nzcs     = set() # TODO: Can this be empty until we decide to optimise?
    psi_tables_code = _tabulate_psis(ir["unique_tables"], used_psi_tables, ir["name_map"], used_nzcs, ir["optimise_parameters"])

    # Delegate to uflacs to generate tabulate_tensor body
    import uflacs.backends.ffc
    ttcode = uflacs.backends.ffc.generate_tabulate_tensor_code(ir, parameters)

    # TODO: Indent psi_tables_code? Or pass psi_tables_code to uflacs for insertion there?
    code["tabulate_tensor"] = '\n\n'.join(psi_tables_code + [ttcode])

    code["tabulate_tensor_quadrature"] = format["do nothing"] # TODO: Remove
    return code
