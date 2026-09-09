# Copyright (C) 2020 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Test configuration."""

import os
import sys

import basix
import numpy as np
import pytest

from ffcx.ir.elementtables import (
    permute_quadrature_interval,
    permute_quadrature_quadrilateral,
    permute_quadrature_triangle,
)


@pytest.fixture(autouse=True, scope="session")
def add_cwd_to_syspath():
    """Ensure the current working directory is in sys.path.

    CFFI and ffcx.main compile/generate files into CWD by default. Without
    this, importlib.import_module cannot find those files when ``pytest``
    is invoked from a directory that is not already on sys.path
    (e.g. the project root rather than the test/ subdirectory).
    """
    cwd = os.getcwd()
    if cwd not in sys.path:
        sys.path.append(cwd)


@pytest.fixture(scope="module")
def compile_args():
    """Compiler arguments."""
    if sys.platform.startswith("win32"):
        return ["-Od"]
    else:
        return ["-O1", "-Wall", "-Werror"]


_ENTITY_ORIENTATION_VARIANTS = {
    # For each entity (facet/ridge) reference cell type: the
    # `permute_quadrature_*` function and the (reflections, rotations)
    # keyword combinations for it, enumerated in the exact order
    # `ffcx/ir/elementtables.py`'s `build_optimized_tables` already
    # permutes mixed-dimensional-submesh element tables over for that
    # same entity type -- `quadrature_permutation` must index into
    # *this* ordering, which is the one DOLFINx and FFCx already agree
    # on, not basix's own (different) `permute_subentity_closure`
    # numbering.
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


def _reference_vertex_permutations(entity_celltype):
    """Return the reference-vertex permutation induced by each orientation variant.

    Returned as a list (one entry per variant, see
    `_ENTITY_ORIENTATION_VARIANTS`) of lists `perm` with
    `perm[new_local_vertex] == old_local_vertex`.

    This is the tests' independent ground truth for
    `facet_closure_dofs`/`ridge_closure_dofs`'s vertex ordering, derived
    only from the already-trusted `permute_quadrature_*` functions --
    with no dependence on basix's `permute_subentity_closure`/`_inv`,
    which `geometry._closure_dofs_table` actually uses to build the
    production tables. It exists to catch exactly the regression that
    motivated it: basix's own `permute_subentity_closure` (forward, not
    `_inv`) numbering diverges from this convention starting at variant
    index 2 for a triangle facet, which silently produced a wrong
    mixed-dimensional-submesh gradient (caught by a 3D
    tetrahedron/hexahedron DOLFINx test) before `_closure_dofs_table`
    was switched to `permute_subentity_closure_inv`. Comparing the
    production table's vertex ordering against this function -- rather
    than only checking the table is internally self-consistent -- is
    what would catch a similar mismatch again, e.g. if a future basix
    release changed that numbering.
    """
    permute_fn, variants = _ENTITY_ORIENTATION_VARIANTS[entity_celltype]
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
