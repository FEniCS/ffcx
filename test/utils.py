# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Shared helpers for the FFCx test suite."""

import basix
import numpy as np

from ffcx.ir.elementtables import (
    permute_quadrature_interval,
    permute_quadrature_quadrilateral,
    permute_quadrature_triangle,
)

ENTITY_ORIENTATION_VARIANTS = {
    # Per entity cell type: the `permute_quadrature_*` function and its
    # kwargs, in the order `build_optimized_tables` permutes over --
    # which is what `quadrature_permutation` indexes, not basix's own
    # `permute_subentity_closure` numbering.
    basix.CellType.interval: (
        permute_quadrature_interval,
        [{"reflections": ref} for ref in range(2)],
    ),
    basix.CellType.triangle: (
        permute_quadrature_triangle,
        [{"reflections": ref, "rotations": rot} for rot in range(3) for ref in range(2)],
    ),
    basix.CellType.quadrilateral: (
        permute_quadrature_quadrilateral,
        [{"reflections": ref, "rotations": rot} for rot in range(4) for ref in range(2)],
    ),
}


def reference_vertex_permutations(entity_celltype):
    """Return the reference-vertex permutation induced by each orientation variant.

    One entry per variant of `ENTITY_ORIENTATION_VARIANTS`, each a list
    `perm` with `perm[new_local_vertex] == old_local_vertex`.

    Independent ground truth for the closure tables' vertex ordering,
    derived only from `permute_quadrature_*` and deliberately not from
    basix's `permute_subentity_closure`, which the production tables use.
    The two diverge from variant 2 of a triangle facet onwards, which
    once produced a silently wrong submesh gradient; comparing against
    this catches a recurrence.
    """
    permute_fn, variants = ENTITY_ORIENTATION_VARIANTS[entity_celltype]
    ref_vertices = np.asarray(basix.geometry(entity_celltype))
    dim = ref_vertices.shape[1]
    permutations = []
    for kwargs in variants:
        mapped = []
        for k in range(len(ref_vertices)):
            permuted_point = np.asarray(permute_fn(ref_vertices[k : k + 1].copy(), **kwargs))[0][
                :dim
            ]
            (old_vertex,) = (
                j for j, v in enumerate(ref_vertices) if np.allclose(permuted_point, v)
            )
            mapped.append(old_vertex)
        permutations.append(mapped)
    return permutations
