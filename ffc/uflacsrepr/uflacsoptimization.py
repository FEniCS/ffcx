# Copyright (C) 2013-2014 Martin Alnaes
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

def optimize_integral_ir(ir, parameters):
    "Compute optimized intermediate representation of integral."

    info("Optimizing uflacs representation")

    # Call upon uflacs to optimize ssa representation prior to code generation.
    # Should be possible to skip this step.
    import uflacs.backends.ffc  # Don't laugh on the next line here :-p
    oir = uflacs.backends.ffc.ffc_uflacs_representation.optimize_tabulate_tensor_ir(ir, parameters)

    return oir
