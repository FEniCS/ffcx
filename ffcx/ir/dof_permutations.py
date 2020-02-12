# Copyright (C) 2020 Matthew Scroggs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import math


def base_permutations(ufl_element, fiat_element):
    cname = ufl_element.cell().cellname()
    if cname == 'point':
        return []
    if cname == 'interval':
        return base_permutations_interval(ufl_element, fiat_element)
    if cname == 'triangle':
        return base_permutations_triangle(ufl_element, fiat_element)
    if cname == 'tetrahedron':
        return base_permutations_tetrahedron(ufl_element, fiat_element)
    if cname == 'quadrilateral':
        return base_permutations_quadrilateral(ufl_element, fiat_element)
    if cname == 'hexahedron':
        return base_permutations_hexahedron(ufl_element, fiat_element)

    raise ValueError("Unrecognised cell type")


def base_permutations_interval(ufl_element, fiat_element):
    num_dofs = len(fiat_element.dual_basis())
    dof_types = [e.functional_type for e in fiat_element.dual_basis()]
    if all_equal(dof_types, "PointEval"):
        return [edge_flip(range(num_dofs))]

    # TODO: What to do if the dofs are not all PointEvals

    return empty_permutations(1, num_dofs)


def base_permutations_triangle(ufl_element, fiat_element):
    num_dofs = len(fiat_element.dual_basis())
    dof_types = [e.functional_type for e in fiat_element.dual_basis()]

    perms = empty_permutations(5, num_dofs)
    if all_equal(dof_types, "PointEval"):
        entity_dofs = fiat_element.entity_dofs()

        # Edge flips
        for edge, dofs in entity_dofs[1].items():
            for i, j in zip(dofs, edge_flip(dofs)):
                perms[edge][i] = j

        # Face rotation and reflection
        dofs = entity_dofs[2][0]
        for i, j in zip(dofs, triangle_rotation(dofs)):
            perms[3][i] = j
        for i, j in zip(dofs, triangle_reflection(dofs)):
            perms[4][i] = j

    # TODO: What to do if the dofs are not all PointEvals
    return perms


def base_permutations_quadrilateral(ufl_element, fiat_element):
    num_dofs = len(fiat_element.dual_basis())
    dof_types = [e.functional_type for e in fiat_element.dual_basis()]

    perms = empty_permutations(6, num_dofs)
    if all_equal(dof_types, "PointEval"):
        entity_dofs = fiat_element.entity_dofs()

        # Edge flips
        for edge, dofs in entity_dofs[1].items():
            for i, j in zip(dofs, edge_flip(dofs)):
                perms[edge][i] = j

        # Face rotation and reflection
        dofs = entity_dofs[2][0]
        for i, j in zip(dofs, quadrilateral_rotation(dofs)):
            perms[4][i] = j
        for i, j in zip(dofs, quadrilateral_reflection(dofs)):
            perms[5][i] = j

    # TODO: What to do if the dofs are not all PointEvals
    return perms


def base_permutations_tetrahedron(ufl_element, fiat_element):
    num_dofs = len(fiat_element.dual_basis())
    dof_types = [e.functional_type for e in fiat_element.dual_basis()]

    perms = empty_permutations(18, num_dofs)
    if all_equal(dof_types, "PointEval"):
        entity_dofs = fiat_element.entity_dofs()

        # Edge flips
        for edge, dofs in entity_dofs[1].items():
            for i, j in zip(dofs, edge_flip(dofs)):
                perms[edge][i] = j

        # Face rotations and reflections
        for face, dofs in entity_dofs[2].items():
            for i, j in zip(dofs, triangle_rotation(dofs)):
                perms[6 + 2 * face][i] = j
            for i, j in zip(dofs, triangle_reflection(dofs)):
                perms[7 + 2 * face][i] = j

        # Volume rotations and reflection
        dofs = entity_dofs[3][0]
        for n, r in enumerate(tetrahedron_rotations(dofs)):
            for i, j in zip(dofs, r):
                perms[14 + n][i] = j
        for i, j in zip(dofs, tetrahedron_reflection(dofs)):
            perms[17][i] = j

    # TODO: What to do if the dofs are not all PointEvals
    return perms


def base_permutations_hexahedron(ufl_element, fiat_element):
    num_dofs = len(fiat_element.dual_basis())
    dof_types = [e.functional_type for e in fiat_element.dual_basis()]

    perms = empty_permutations(28, num_dofs)
    if all_equal(dof_types, "PointEval"):
        entity_dofs = fiat_element.entity_dofs()

        # Edge flips
        for edge, dofs in entity_dofs[1].items():
            for i, j in zip(dofs, edge_flip(dofs)):
                perms[edge][i] = j

        # Face rotations and reflections
        for face, dofs in entity_dofs[2].items():
            for i, j in zip(dofs, quadrilateral_rotation(dofs)):
                perms[12 + 2 * face][i] = j
            for i, j in zip(dofs, quadrilateral_reflection(dofs)):
                perms[13 + 2 * face][i] = j

        # Volume rotations and reflection
        dofs = entity_dofs[3][0]
        for n, r in enumerate(hexahedron_rotations(dofs)):
            for i, j in zip(dofs, r):
                perms[24 + n][i] = j
        for i, j in zip(dofs, hexahedron_reflection(dofs)):
            perms[27][i] = j

    # TODO: What to do if the dofs are not all PointEvals
    return perms


def empty_permutations(num_perms, num_dofs):
    return [list(range(num_dofs)) for i in range(num_perms)]


def all_equal(ls, val):
    for i in ls:
        if i != val:
            return False
    return True


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
    s = (math.floor(math.sqrt(1 + 4 * n ** 2)) - 1) // 2
    assert s * (s + 1) == 2 * n

    perm = []
    for i in range(n):
        dof = n - 1 - i
        for sub in range(i, s):
            for k in range(blocksize):
                perm.append(blocksize * dof + k)
            dof -= sub + 1
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def triangle_reflection(dofs, blocksize=1):
    n = len(dofs) // blocksize
    s = (math.floor(math.sqrt(1 + 4 * n ** 2)) - 1) // 2
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
        for dof in range(st, -1, -1):
            for k in range(blocksize):
                perm.append(blocksize * dof + k)
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def quadrilateral_reflection(dofs, blocksize=1):
    n = len(dofs) // blocksize
    s = math.floor(math.sqrt(n))
    assert s ** 2 == n

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
    for side in range(s, -1, -1):
        face_dofs = list(range(len(rot1), len(rot1) + n * s * (s + 1) // 2))
        rot1 += triangle_rotation(face_dofs, blocksize)
    assert len(rot1) == len(dofs)

    rot2 = []
    for side in range(s, -1, -1):
        face_dofs = []
        start = side * s - 1 - s * (s - 1) // 2
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
            st += i * (i + 1) / 2 - layer
        layerst += s - layer
    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def hexahedron_rotations(dofs, blocksize=1):
    n = len(dofs) // blocksize
    s = math.floor(math.cbrt(n))
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
    for st in range(area):
        for dof in range(st, n, area):
            for k in range(blocksize):
                rot2.append(blocksize * dof + k)
    assert len(rot3) == len(dofs)

    return [[dofs[i] for i in rot1], [dofs[i] for i in rot2], [dofs[i] for i in rot3]]


def hexahedron_reflection(dofs, blocksize=1):
    n = len(dofs) // blocksize
    s = math.floor(math.cbrt(n))
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
