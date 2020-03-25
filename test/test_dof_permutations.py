# Copyright (C) 2020 Matthew W. Scroggs
#
# This file is part of FFCX.
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
"Unit tests for FFCX"


import pytest

from ffcx.ir import dof_permutations


@pytest.mark.parametrize("s", range(1, 5))
@pytest.mark.parametrize("blocksize", range(1, 5))
@pytest.mark.parametrize("perm_f, dof_f", [
    (dof_permutations.edge_flip, lambda s: s),
    (dof_permutations.triangle_rotation, lambda s: s * (s + 1) // 2),
    (dof_permutations.triangle_reflection, lambda s: s * (s + 1) // 2),
    (dof_permutations.quadrilateral_rotation, lambda s: s ** 2),
    (dof_permutations.quadrilateral_reflection, lambda s: s ** 2)])
def test_permutations(s, perm_f, dof_f, blocksize):
    """Test that permutations are valid."""
    dofs = dof_f(s) * blocksize
    perms = perm_f(list(range(dofs)), blocksize)
    if not isinstance(perms[0], list):
        perms = [perms]
    for p in perms:
        for i in p:
            # Each number must appear exactly once in the permutation
            assert p.count(i) == 1
            # Each number must be less than the number of dofs
            assert i < dofs
