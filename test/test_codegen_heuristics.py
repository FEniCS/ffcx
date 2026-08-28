# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Tests for code-generation size and storage heuristics."""

import pytest

from ffcx.codegeneration.definitions import _should_unroll_coordinate_dofs
from ffcx.codegeneration.integral_generator import _fw_cache_tile_size


def test_coordinate_dof_unrolling_is_limited_to_small_non_tensor_geometry():
    """Curved simplex geometry must not be expanded into an unbounded sum."""
    assert _should_unroll_coordinate_dofs(10, False)
    assert not _should_unroll_coordinate_dofs(17, False)
    assert not _should_unroll_coordinate_dofs(4, True)


@pytest.mark.parametrize(
    "n_fw, num_points, expected",
    [
        (1, 10000, 2048),
        (256, 10000, 8),
        (257, 10000, 7),
        (2048, 10000, 1),
        (2049, 10000, None),
        (1000, 4, 2),
    ],
)
def test_fw_cache_tiles_do_not_exceed_the_stack_budget(n_fw, num_points, expected):
    """The cache must remain bounded even with many fw intermediates."""
    tile = _fw_cache_tile_size(n_fw, num_points)
    assert tile == expected
    if tile is not None:
        assert n_fw * tile * 16 <= 32768
