# Copyright (C) 2020 Matthew W. Scroggs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import math
import numpy as np
from ffcx.fiatinterface import create_element

# TODO: currently these dof types are not correctly handled:
#       FrobeniusIntegralMoment
# TODO: currently these dof types are not handled at all:
#       PointEdgeTangent
#       PointFaceTangent
#       IntegralMomentOfNormalDerivative


def base_permutations_and_reflection_entities(ufl_element):
    # If this element has no sub elements
    if ufl_element.num_sub_elements() == 0:
        return base_permutations_from_subdofmap(ufl_element)

    # If this element has sub elements, combine the permutations from them all
    perms = None
    reflections = []
    rotations = []
    for e in ufl_element.sub_elements():
        bp, re, ro = base_permutations_and_reflection_entities(e)
        if perms is None:
            perms = [[] for i in bp]
        if len(bp) > 0:
            start = len(perms[0])
            for i, b in enumerate(bp):
                perms[i] += [a + start for a in b]
        reflections += re
        rotations += ro
    return perms, reflections, rotations


def base_permutations_from_subdofmap(ufl_element):
    # FIXME: This information should be added to FIAT instead of reverse engineered here
    cname = ufl_element.cell().cellname()
    if cname == 'point':
        return []

    if cname == 'interval':
        entity_counts = [2, 1, 0, 0]
        entity_functions = [None, permute_edge, None, None]
    elif cname == 'triangle':
        entity_counts = [3, 3, 1, 0]
        entity_functions = [None, permute_edge, permute_triangle, None]
    elif cname == 'tetrahedron':
        entity_counts = [4, 6, 4, 1]
        entity_functions = [None, permute_edge, permute_triangle, permute_tetrahedron]
    elif cname == 'quadrilateral':
        entity_counts = [4, 4, 1, 0]
        entity_functions = [None, permute_edge, permute_quadrilateral, None]
    elif cname == 'hexahedron':
        entity_counts = [8, 12, 6, 1]
        entity_functions = [None, permute_edge, permute_quadrilateral, permute_hexahedron]
    else:
        raise ValueError("Unrecognised cell type")

    fiat_element = create_element(ufl_element)

    num_dofs = len(fiat_element.dual_basis())
    dof_types = [e.functional_type for e in fiat_element.dual_basis()]
    entity_dofs = fiat_element.entity_dofs()
    num_perms = entity_counts[1] + 2 * entity_counts[2] + 4 * entity_counts[3]

    perms = empty_permutations(num_perms, num_dofs)
    reflections = [None for i in range(num_dofs)]
    rotations = []
    perm_n = 0
    for dim in range(1, 4):
        for n in range(entity_counts[dim]):
            dofs = entity_dofs[dim][n]
            types = [dof_types[i] for i in dofs]
            unique_types = []
            for t in types:
                if t not in unique_types:
                    unique_types.append(t)
            for t in unique_types:
                type_dofs = [i for i, j in zip(dofs, types) if j == t]
                if t in ["PointScaledNormalEval", "ComponentPointEval", "PointEdgeTangent",
                         "PointScaledNormalEval", "PointNormalEval",
                         "IntegralMoment"]:
                    for i in type_dofs:
                        reflections[i] = [(dim, n)]
                elif t == "PointFaceTangent" and cname in ["triangle","tetrahedron"]:
                    assert dim == 2
                    for dofs in zip(type_dofs[::2], type_dofs[1::2]):
                        # entity dimension, entity number, first dof, second dof,
                        #   matrix that represents a rotation, order
                        rotations.append((dim, n, dofs, np.array([[0, -1], [1, -1]]), 3))
                elif t == "FrobeniusIntegralMoment":
                    if dim == 2:
                        if len(type_dofs) != 3:
                            # FIXME
                            raise ValueError("Not yet implemented")
                        for a, b in zip(type_dofs, fiat_element.ref_el.connectivity[dim, dim - 1][n]):
                            reflections[a] = [(dim, n), (dim - 1, b)]

                if t in ["PointEval", "PointNormalDeriv", "PointEdgeTangent",
                         "PointDeriv", "PointNormalEval", "PointScaledNormalEval"]:
                    # Dof is a point evaluation, use blocksize 1
                    permuted = entity_functions[dim](type_dofs, 1)
                elif t in ["ComponentPointEval", "IntegralMoment"]:
                    # Dof blocksize is equal to entity dimension
                    permuted = entity_functions[dim](type_dofs, dim)
                elif t == "PointFaceTangent":
                    # Dof blocksize is 2
                    permuted = entity_functions[dim](type_dofs, 2, True)
                elif t in ["FrobeniusIntegralMoment"] and dim == 2:
                    permuted = permute_frobenius_face(type_dofs, 1)
                elif t in ["FrobeniusIntegralMoment", "PointwiseInnerProductEval"]:
                    # FIXME: temporarily does no permutation; needs replacing
                    permuted = [type_dofs for i in range(2 ** (dim - 1))]
                else:
                    # TODO: What to do with other dof types
                    raise ValueError("Permutations are not currently implemented for this dof type (" + t + ").")

                # Apply these permutations
                for p in range(2 ** (dim - 1)):
                    for i, j in zip(type_dofs, permuted[p]):
                        perms[perm_n + p][i] = j
            perm_n += 2 ** (dim - 1)

    return perms, reflections, rotations


def permute_frobenius_face(dofs, blocksize):
    n = len(dofs) // blocksize
    s = 0
    while 6 * s + s * (s - 1) < 2 * n:
        s += 1
    assert 6 * s + s * (s - 1) == 2 * n

    rot = []
    for dof in range(2 * s, 3 * s):
        rot += [dof * blocksize + k for k in range(blocksize)]
    for dof in range(s - 1, -1, -1):
        rot += [dof * blocksize + k for k in range(blocksize)]
    for dof in range(2 * s - 1, s - 1, -1):
        rot += [dof * blocksize + k for k in range(blocksize)]
    rot += triangle_rotation(list(range(3 * s * blocksize, n)), blocksize)
    assert len(rot) == len(dofs)

    ref = []
    for dof in range(s - 1, -1, -1):
        ref += [dof * blocksize + k for k in range(blocksize)]
    for dof in range(2 * s, 3 * s):
        ref += [dof * blocksize + k for k in range(blocksize)]
    for dof in range(s, 2 * s):
        ref += [dof * blocksize + k for k in range(blocksize)]
    ref += triangle_reflection(list(range(3 * s * blocksize, n)), blocksize)
    assert len(ref) == len(dofs)

    return [[dofs[i] for i in rot],
            [dofs[i] for i in ref]]


def permute_edge(dofs, blocksize, reverse_blocks=False):
    return [edge_flip(dofs, blocksize, reverse_blocks)]


def permute_triangle(dofs, blocksize, reverse_blocks=False):
    return [triangle_rotation(dofs, blocksize), triangle_reflection(dofs, blocksize, reverse_blocks)]


def permute_quadrilateral(dofs, blocksize, reverse_blocks=False):
    return [quadrilateral_rotation(dofs, blocksize), quadrilateral_reflection(dofs, blocksize, reverse_blocks)]


def permute_tetrahedron(dofs, blocksize, reverse_blocks=False):
    return tetrahedron_rotations(dofs, blocksize) + [tetrahedron_reflection(dofs, blocksize, reverse_blocks)]


def permute_hexahedron(dofs, blocksize, reverse_blocks=False):
    return hexahedron_rotations(dofs, blocksize) + [hexahedron_reflection(dofs, blocksize, reverse_blocks)]


def empty_permutations(num_perms, num_dofs, reverse_blocks=False):
    return [list(range(num_dofs)) for i in range(num_perms)]


def edge_flip(dofs, blocksize=1, reverse_blocks=False):
    n = len(dofs) // blocksize

    perm = []
    for dof in range(n - 1, -1, -1):
        if reverse_blocks:
            perm += [dof * blocksize + k for k in range(blocksize)][::-1]
        else:
            perm += [dof * blocksize + k for k in range(blocksize)]

    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def triangle_rotation(dofs, blocksize=1, reverse_blocks=False):
    n = len(dofs) // blocksize
    s = (math.floor(math.sqrt(1 + 8 * n)) - 1) // 2
    assert s * (s + 1) == 2 * n

    perm = []
    st = n - 1
    for i in range(1, s + 1):
        dof = st
        for sub in range(i, s + 1):
            if reverse_blocks:
                perm += [dof * blocksize + k for k in range(blocksize)][::-1]
            else:
                perm += [dof * blocksize + k for k in range(blocksize)]
            dof -= sub + 1
        st -= i
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def triangle_reflection(dofs, blocksize=1, reverse_blocks=False):
    n = len(dofs) // blocksize
    s = (math.floor(math.sqrt(1 + 8 * n)) - 1) // 2
    assert s * (s + 1) == 2 * n

    perm = []
    for st in range(s):
        dof = st
        for add in range(s, st, -1):
            if reverse_blocks:
                perm += [dof * blocksize + k for k in range(blocksize)][::-1]
            else:
                perm += [dof * blocksize + k for k in range(blocksize)]
            dof += add
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def quadrilateral_rotation(dofs, blocksize=1, reverse_blocks=False):
    n = len(dofs) // blocksize
    s = math.floor(math.sqrt(n))
    assert s ** 2 == n

    perm = []
    for st in range(n - s, n):
        for dof in range(st, -1, -s):
            if reverse_blocks:
                perm += [dof * blocksize + k for k in range(blocksize)][::-1]
            else:
                perm += [dof * blocksize + k for k in range(blocksize)]
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def quadrilateral_reflection(dofs, blocksize=1, reverse_blocks=False):
    n = len(dofs) // blocksize
    s = math.floor(math.sqrt(n))
    assert s ** 2 == n
    if s == 0:
        return dofs

    perm = []
    for st in range(s):
        for dof in range(st, n, s):
            if reverse_blocks:
                perm += [dof * blocksize + k for k in range(blocksize)][::-1]
            else:
                perm += [dof * blocksize + k for k in range(blocksize)]
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def tetrahedron_rotations(dofs, blocksize=1, reverse_blocks=False):
    n = len(dofs) // blocksize
    s = 0
    while s * (s + 1) * (s + 2) < 6 * n:
        s += 1
    assert s * (s + 1) * (s + 2) == 6 * n

    rot1 = []
    for side in range(s, 0, -1):
        face_dofs = list(range(len(rot1), len(rot1) + blocksize * side * (side + 1) // 2))
        rot1 += triangle_rotation(face_dofs, blocksize, reverse_blocks)
    assert len(rot1) == len(dofs)

    rot2 = []
    for side in range(s, -1, -1):
        face_dofs = []
        start = side * s - 1 - side * (side - 1) // 2
        for row in range(side):
            dof = start
            for k in range(blocksize):
                face_dofs.append(dof * blocksize + k)
            for sub in range(2 + (s - side), s - row + 1):
                dof -= sub
                for k in range(blocksize):
                    face_dofs.append(dof * blocksize + k)
            start += (s - row - 1) * (s - row) // 2

        rot2 += triangle_rotation(face_dofs, blocksize, reverse_blocks)
    assert len(rot2) == len(dofs)

    rot3 = [rot2[rot2[j]] for j in rot1]
    assert len(rot3) == len(dofs)

    return [[dofs[i] for i in rot1], [dofs[i] for i in rot2], [dofs[i] for i in rot3]]


def tetrahedron_reflection(dofs, blocksize=1, reverse_blocks=False):
    n = len(dofs) // blocksize
    s = 0
    while s * (s + 1) * (s + 2) < 6 * n:
        s += 1
    assert s * (s + 1) * (s + 2) == 6 * n

    perm = []
    layerst = 0
    for layer in range(s):
        st = layerst
        for i in range(s, layer - 1, -1):
            for dof in range(st, st + i - layer):
                if reverse_blocks:
                    perm += [dof * blocksize + k for k in range(blocksize)][::-1]
                else:
                    perm += [dof * blocksize + k for k in range(blocksize)]
            st += i * (i + 1) // 2 - layer
        layerst += s - layer
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def hexahedron_rotations(dofs, blocksize=1, reverse_blocks=False):
    n = len(dofs) // blocksize
    s = 0
    while s ** 3 < n:
        s += 1
    area = s ** 2
    assert s ** 3 == n

    rot1 = []
    for lst in range(area - s, n, area):
        for st in range(lst, lst + s):
            for dof in range(st, lst - area + s - 1, -s):
                if reverse_blocks:
                    rot1 += [dof * blocksize + k for k in range(blocksize)][::-1]
                else:
                    rot1 += [dof * blocksize + k for k in range(blocksize)]
    assert len(rot1) == len(dofs)

    rot2 = []
    for lst in range(s - 1, -1, -1):
        for st in range(lst, area, s):
            for dof in range(st, n, area):
                if reverse_blocks:
                    rot2 += [dof * blocksize + k for k in range(blocksize)][::-1]
                else:
                    rot2 += [dof * blocksize + k for k in range(blocksize)]
    assert len(rot2) == len(dofs)

    rot3 = []
    for st in range(0, area):
        for dof in range(st, n, area):
            if reverse_blocks:
                rot3 += [dof * blocksize + k for k in range(blocksize)][::-1]
            else:
                rot3 += [dof * blocksize + k for k in range(blocksize)]
    assert len(rot3) == len(dofs)

    return [[dofs[i] for i in rot1], [dofs[i] for i in rot2], [dofs[i] for i in rot3]]


def hexahedron_reflection(dofs, blocksize=1, reverse_blocks=False):
    n = len(dofs) // blocksize
    s = 0
    while s ** 3 < n:
        s += 1
    area = s ** 2
    assert s ** 3 == n

    perm = []
    for lst in range(0, area, s):
        for st in range(lst, n, area):
            for dof in range(st, st + s):
                if reverse_blocks:
                    perm += [dof * blocksize + k for k in range(blocksize)][::-1]
                else:
                    perm += [dof * blocksize + k for k in range(blocksize)]
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]
