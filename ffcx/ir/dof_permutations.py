# Copyright (C) 2020 Matthew W. Scroggs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import math
from ffcx.fiatinterface import create_element

# TODO: currently these dof types are not correctly handled:
#       FrobeniusIntegralMoment
# TODO: currently these dof types are not handled at all:
#       PointEdgeTangent
#       PointFaceTangent
#       IntegralMomentOfNormalDerivative
#       PointwiseInnerProductEval


def base_permutations(ufl_element):
    if ufl_element.num_sub_elements() == 0:
        return base_permutations_from_subdofmap(ufl_element)

    perms = None
    for e in ufl_element.sub_elements():
        bp = base_permutations(e)
        if perms is None:
            perms = [[] for i in bp]
        if len(bp) > 0:
            start = len(perms[0])
            for i, b in enumerate(bp):
                perms[i] += [a + start for a in b]
    return perms


def base_permutations_from_subdofmap(ufl_element):
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
                if t in ["PointEval", "PointNormalDeriv", "PointEdgeTangent",
                         "PointScaledNormalEval", "PointDeriv", "PointNormalEval"]:
                    # Dof is a point evaluation, use blocksize 1
                    permuted = entity_functions[dim](dofs, 1)
                elif t in ["ComponentPointEval", "IntegralMoment"]:
                    # Dof blocksize is equal to entity dimension
                    permuted = entity_functions[dim](dofs, dim)
                elif t == "PointFaceTangent":
                    # Dof blocksize is 2
                    permuted = entity_functions[dim](dofs, 2)
                elif t == "FrobeniusIntegralMoment":
                    # FIXME: temporarily does no permutation; needs replacing
                    permuted = [dofs for i in range(2 ** (dim - 1))]
                else:
                    # TODO: What to do with other dof types
                    raise ValueError("Permutations are not currently implemented for this dof type (" + t + ").")

                # Apply these permutations
                for p in range(2 ** (dim - 1)):
                    for i, j in zip(type_dofs, permuted[p]):
                        perms[perm_n + p][i] = j
            perm_n += 2 ** (dim - 1)

    return perms


def permute_edge(dofs, blocksize):
    return [edge_flip(dofs, blocksize)]


def permute_triangle(dofs, blocksize):
    return [triangle_rotation(dofs, blocksize), triangle_reflection(dofs, blocksize)]


def permute_quadrilateral(dofs, blocksize):
    return [quadrilateral_rotation(dofs, blocksize), quadrilateral_reflection(dofs, blocksize)]


def permute_tetrahedron(dofs, blocksize):
    return tetrahedron_rotations(dofs, blocksize) + [tetrahedron_reflection(dofs, blocksize)]


def permute_hexahedron(dofs, blocksize):
    return hexahedron_rotations(dofs, blocksize) + [hexahedron_reflection(dofs, blocksize)]


def empty_permutations(num_perms, num_dofs):
    return [list(range(num_dofs)) for i in range(num_perms)]


def edge_flip(dofs, blocksize=1):
    n = len(dofs) // blocksize

    perm = []
    for i in range(n - 1, -1, -1):
        for k in range(blocksize):
            perm.append(i * blocksize + k)
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def triangle_rotation(dofs, blocksize=1):
    n = len(dofs) // blocksize
    s = (math.floor(math.sqrt(1 + 8 * n)) - 1) // 2
    assert s * (s + 1) == 2 * n

    perm = []
    st = n - 1
    for i in range(1, s + 1):
        dof = st
        for sub in range(i, s + 1):
            for k in range(blocksize):
                perm.append(blocksize * dof + k)
            dof -= sub + 1
        st -= i
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def triangle_reflection(dofs, blocksize=1):
    n = len(dofs) // blocksize
    s = (math.floor(math.sqrt(1 + 8 * n)) - 1) // 2
    assert s * (s + 1) == 2 * n

    perm = []
    for st in range(s):
        dof = st
        for add in range(s, st, -1):
            for k in range(blocksize):
                perm.append(blocksize * dof + k)
            dof += add
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def quadrilateral_rotation(dofs, blocksize=1):
    n = len(dofs) // blocksize
    s = math.floor(math.sqrt(n))
    assert s ** 2 == n

    perm = []
    for st in range(n - s, n):
        for dof in range(st, -1, -s):
            for k in range(blocksize):
                perm.append(blocksize * dof + k)
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def quadrilateral_reflection(dofs, blocksize=1):
    n = len(dofs) // blocksize
    s = math.floor(math.sqrt(n))
    assert s ** 2 == n
    if s == 0:
        return dofs

    perm = []
    for st in range(s):
        for dof in range(st, n, s):
            for k in range(blocksize):
                perm.append(blocksize * dof + k)
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def tetrahedron_rotations(dofs, blocksize=1):
    n = len(dofs) // blocksize
    s = 0
    while s * (s + 1) * (s + 2) < 6 * n:
        s += 1
    assert s * (s + 1) * (s + 2) == 6 * n

    rot1 = []
    for side in range(s, 0, -1):
        face_dofs = list(range(len(rot1), len(rot1) + blocksize * side * (side + 1) // 2))
        rot1 += triangle_rotation(face_dofs, blocksize)
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

        rot2 += triangle_rotation(face_dofs, blocksize)
    assert len(rot2) == len(dofs)

    rot3 = [rot2[rot2[j]] for j in rot1]
    assert len(rot3) == len(dofs)

    return [[dofs[i] for i in rot1], [dofs[i] for i in rot2], [dofs[i] for i in rot3]]


def tetrahedron_reflection(dofs, blocksize=1):
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
                for k in range(blocksize):
                    perm.append(dof * blocksize + k)
            st += i * (i + 1) // 2 - layer
        layerst += s - layer
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def hexahedron_rotations(dofs, blocksize=1):
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
                for k in range(blocksize):
                    rot1.append(blocksize * dof + k)
    assert len(rot1) == len(dofs)

    rot2 = []
    for lst in range(s - 1, -1, -1):
        for st in range(lst, area, s):
            for dof in range(st, n, area):
                for k in range(blocksize):
                    rot2.append(blocksize * dof + k)
    assert len(rot2) == len(dofs)

    rot3 = []
    for st in range(0, area):
        for dof in range(st, n, area):
            for k in range(blocksize):
                rot3.append(blocksize * dof + k)
        print(rot3)
    print(rot3, dofs)
    assert len(rot3) == len(dofs)

    return [[dofs[i] for i in rot1], [dofs[i] for i in rot2], [dofs[i] for i in rot3]]


def hexahedron_reflection(dofs, blocksize=1):
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
                for k in range(blocksize):
                    perm.append(blocksize * dof + k)
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]
