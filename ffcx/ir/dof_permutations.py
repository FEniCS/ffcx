# Copyright (C) 2020 Matthew W. Scroggs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import warnings
import math
import ufl
from ffcx.fiatinterface import create_element

# TODO: This information should be moved to FIAT instead of being reverse engineered here
# TODO: Currently none of the vector-valued stuff has been tested on quads and hexes
# TODO: currently these dof types are not handled at all:
#       PointwiseInnerProductEval
#       IntegralMomentOfNormalDerivative


def base_permutations(ufl_element):
    """Returns the base permutations."""
    if ufl_element.num_sub_elements() == 0:
        # If the element has no sub elements, return its permutations
        return base_permutations_from_subdofmap(ufl_element)

    # If the element has sub elements, combine their permutations
    perms = None

    if isinstance(ufl_element, ufl.VectorElement) or isinstance(ufl_element, ufl.TensorElement):
        block_size = ufl_element.num_sub_elements()
        return [
            [block_size * j + i for j in perm for i in range(block_size)]
            for perm in base_permutations(ufl_element.sub_elements()[0])
        ]

    for e in ufl_element.sub_elements():
        bp = base_permutations(e)
        if perms is None:
            perms = [[] for i in bp]
        for i, b in enumerate(bp):
            perms[i] += [a + len(perms[i]) for a in b]
    return perms


def reflection_entities(ufl_element):
    """Returns the entities that the direction of vector-valued functions depend on."""
    if ufl_element.num_sub_elements() == 0:
        # If the element has no sub elements, return its reflection entities
        return reflection_entities_from_subdofmap(ufl_element)

    if isinstance(ufl_element, ufl.VectorElement) or isinstance(ufl_element, ufl.TensorElement):
        block_size = ufl_element.num_sub_elements()
        return [
            ref
            for ref in reflection_entities(ufl_element.sub_elements()[0])
            for i in range(block_size)
        ]
    # If the element has sub elements, combine their reflections
    reflections = []
    for e in ufl_element.sub_elements():
        reflections += reflection_entities(e)
    return reflections


def face_tangents(ufl_element):
    """Returns the rotations that rotate the direction of vector-valued face tangent dofs."""
    if ufl_element.num_sub_elements() == 0:
        # If the element has no sub elements, return its rotations
        return face_tangents_from_subdofmap(ufl_element)

    if isinstance(ufl_element, ufl.VectorElement) or isinstance(ufl_element, ufl.TensorElement):
        block_size = ufl_element.num_sub_elements()
        return [
            tan
            for tan in face_tangents(ufl_element.sub_elements()[0])
            for i in range(block_size)
        ]

    # If the element has sub elements, combine their rotations
    rotations = []
    for e in ufl_element.sub_elements():
        rotations += face_tangents(e)
    return rotations


def base_permutations_from_subdofmap(ufl_element):
    """Calculate permutations and reflection entites for a root element.
    Calculates the base permutations and the entities that the direction of vector-valued
    functions depend on for an element with no sub elements."""
    fiat_element = create_element(ufl_element)
    dual = fiat_element.dual_basis()
    num_dofs = len(dual)

    tdim = ufl_element.cell().topological_dimension()

    # Get the entity counts and shape of each entity for the cell type
    entity_counts = get_entity_counts(fiat_element)
    entity_functions = get_entity_functions(ufl_element)

    # There is 1 permutation for a 1D entity, 2 for a 2D entity and 4 for a 3D entity
    num_perms = sum([0, 1, 2, 4][i] * j for i, j in enumerate(entity_counts[:tdim]))

    dof_types = [e.functional_type for e in dual]
    entity_dofs = fiat_element.entity_dofs()

    perms = identity_permutations(num_perms, num_dofs)
    perm_n = 0
    # Iterate through the entities of the reference element
    for dim in range(1, tdim):
        for entity_n in range(entity_counts[dim]):
            dofs = entity_dofs[dim][entity_n]
            types = [dof_types[i] for i in dofs]
            # Find the unique dof types
            unique_types = []
            for t in types:
                if t not in unique_types:
                    unique_types.append(t)
            # Permute the dofs of each entity type separately
            for t in unique_types:
                permuted = None
                type_dofs = [i for i, j in zip(dofs, types) if j == t]
                if t in ["PointEval", "PointNormalDeriv", "PointEdgeTangent",
                         "PointDeriv", "PointNormalEval", "PointScaledNormalEval"]:
                    # Dof is a point evaluation, use blocksize 1
                    permuted = entity_functions[dim](type_dofs, 1)
                elif t in ["ComponentPointEval", "IntegralMoment"]:
                    # Dof blocksize is equal to entity dimension
                    permuted = entity_functions[dim](type_dofs, dim)
                elif t == "PointFaceTangent":
                    # Dof blocksize is 2
                    permuted = entity_functions[dim](type_dofs, 2)
                elif t in ["FrobeniusIntegralMoment"] and dim == 2:
                    permuted = permute_frobenius_face(type_dofs, 1)
                else:
                    if dim < tdim:
                        # TODO: What to do with other dof types
                        warnings.warn("Permutations of " + t + " dofs not yet "
                                      "implemented. Results on unordered meshes may be incorrect")
                    continue

                # Apply these permutations
                for n, p in enumerate(permuted):
                    for i, j in zip(type_dofs, p):
                        perms[perm_n + n][i] = j
            perm_n += 2 ** (dim - 1)
    return perms


def reflection_entities_from_subdofmap(ufl_element):
    """Calculate permutations and reflection entites for a root element.
    Calculates the base permutations and the entities that the direction of vector-valued
    functions depend on for an element with no sub elements."""
    fiat_element = create_element(ufl_element)
    dual = fiat_element.dual_basis()
    num_dofs = len(dual)

    cname = ufl_element.cell().cellname()
    tdim = ufl_element.cell().topological_dimension()

    # Get the entity counts for the cell type
    entity_counts = get_entity_counts(fiat_element)

    dof_types = [e.functional_type for e in dual]
    entity_dofs = fiat_element.entity_dofs()

    reflections = [None for i in range(num_dofs)]
    # Iterate through the entities of the reference element
    for dim in range(1, tdim):
        for entity_n in range(entity_counts[dim]):
            dofs = entity_dofs[dim][entity_n]
            types = [dof_types[i] for i in dofs]
            # Find the unique dof types
            unique_types = []
            for t in types:
                if t not in unique_types:
                    unique_types.append(t)
            # Permute the dofs of each entity type separately
            for t in unique_types:
                type_dofs = [i for i, j in zip(dofs, types) if j == t]
                if t in ["PointScaledNormalEval", "ComponentPointEval", "PointEdgeTangent",
                         "PointScaledNormalEval", "PointNormalEval", "IntegralMoment"]:
                    for i in type_dofs:
                        reflections[i] = [(dim, entity_n)]
                elif t == "FrobeniusIntegralMoment" and cname in ["triangle", "tetrahedron"]:
                    if dim == 2:
                        s = get_frobenius_side_length(len(type_dofs))
                        for i, b in enumerate(fiat_element.ref_el.connectivity[dim, dim - 1][entity_n]):
                            for a in type_dofs[i * s:(i + 1) * s]:
                                reflections[a] = [(dim, entity_n), (dim - 1, b)]
    return reflections


def face_tangents_from_subdofmap(ufl_element):
    """Calculate permutations and reflection entites for a root element.
    Calculates the base permutations and the entities that the direction of vector-valued
    functions depend on for an element with no sub elements."""
    fiat_element = create_element(ufl_element)
    dual = fiat_element.dual_basis()
    cname = ufl_element.cell().cellname()

    dof_types = [e.functional_type for e in dual]
    entity_dofs = fiat_element.entity_dofs()

    rotations = []
    if cname == "tetrahedron":
        # Iterate through faces
        for entity_n in range(len(entity_dofs[2])):
            dofs = entity_dofs[2][entity_n]
            types = [dof_types[i] for i in dofs]

            # PointFaceTangent dofs
            if "PointFaceTangent" in types:
                type_dofs = [i for i, t in zip(dofs, types) if t == "PointFaceTangent"]
                for dof_pair in zip(type_dofs[::2], type_dofs[1::2]):
                    # (entity_dim, entity_number), dofs
                    rotations.append(((2, entity_n), dof_pair))
            # FrobeniusIntegralMoment dofs
            if "FrobeniusIntegralMoment" in types:
                type_dofs = [i for i, t in zip(dofs, types) if t == "FrobeniusIntegralMoment"]
                s = get_frobenius_side_length(len(type_dofs))
                for dof_pair in zip(type_dofs[3 * s::2], type_dofs[3 * s + 1::2]):
                    # (entity_dim, entity_number), dofs
                    rotations.append(((2, entity_n), dof_pair))

    return rotations


def get_entity_counts(fiat_element):
    topology = fiat_element.ref_el.topology
    return [len(topology[i]) for i in range(len(topology))]


def get_entity_functions(ufl_element):
    cname = ufl_element.cell().cellname()
    if cname == 'point':
        return [None]
    elif cname == 'interval':
        return [None, None]
    elif cname == 'triangle':
        return [None, permute_edge, None]
    elif cname == 'tetrahedron':
        return [None, permute_edge, permute_triangle, None]
    elif cname == 'quadrilateral':
        return [None, permute_edge, None]
    elif cname == 'hexahedron':
        return [None, permute_edge, permute_quadrilateral, None]
    else:
        raise ValueError("Unrecognised cell type")


def get_frobenius_side_length(n):
    """Get the side length the arrangement of FrobeniusIntegralMoment dofs of a face of a N2curl space."""
    s = 0
    while 3 * s + s * (s - 1) < n:
        s += 1
    assert 3 * s + s * (s - 1) == n
    return s


def permute_frobenius_face(dofs, blocksize):
    """Permute the FrobeniusIntegralMoment dofs of a face of a N2curl space."""
    n = len(dofs) // blocksize
    s = get_frobenius_side_length(n)

    # Make the rotation
    rot = []
    for dof in range(2 * s, 3 * s):
        rot += [dof * blocksize + k for k in range(blocksize)]
    for dof in range(s - 1, -1, -1):
        rot += [dof * blocksize + k for k in range(blocksize)]
    for dof in range(2 * s - 1, s - 1, -1):
        rot += [dof * blocksize + k for k in range(blocksize)]
    rot += triangle_rotation(list(range(3 * s * blocksize, n)), 2 * blocksize)
    assert len(rot) == len(dofs)

    # Make the reflection
    ref = []
    for dof in range(s - 1, -1, -1):
        ref += [dof * blocksize + k for k in range(blocksize)]
    for dof in range(2 * s, 3 * s):
        ref += [dof * blocksize + k for k in range(blocksize)]
    for dof in range(s, 2 * s):
        ref += [dof * blocksize + k for k in range(blocksize)]
    ref += triangle_reflection(list(range(3 * s * blocksize, n)), 2 * blocksize)
    assert len(ref) == len(dofs)

    return [[dofs[i] for i in rot],
            [dofs[i] for i in ref]]


def permute_edge(dofs, blocksize, reverse_blocks=False):
    """Permute the dofs on an edge."""
    return [edge_flip(dofs, blocksize, reverse_blocks)]


def permute_triangle(dofs, blocksize, reverse_blocks=False):
    """Permute the dofs on a triangle."""
    return [triangle_rotation(dofs, blocksize), triangle_reflection(dofs, blocksize, reverse_blocks)]


def permute_quadrilateral(dofs, blocksize, reverse_blocks=False):
    """Permute the dofs on a quadrilateral."""
    return [quadrilateral_rotation(dofs, blocksize), quadrilateral_reflection(dofs, blocksize, reverse_blocks)]


def identity_permutations(num_perms, num_dofs):
    """Return identity permutations of the given shape."""
    return [list(range(num_dofs)) for i in range(num_perms)]


def edge_flip(dofs, blocksize=1, reverse_blocks=False):
    """Flip the dofs on an edge."""
    n = len(dofs) // blocksize

    perm = []
    for dof in range(n - 1, -1, -1):
        if reverse_blocks:
            # Reverse the dofs within a block
            perm += [dof * blocksize + k for k in range(blocksize)][::-1]
        else:
            perm += [dof * blocksize + k for k in range(blocksize)]

    assert len(perm) == len(dofs)

    return [dofs[i] for i in perm]


def triangle_rotation(dofs, blocksize=1, reverse_blocks=False):
    """Rotate the dofs in a triangle."""
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
    """Reflect the dofs in a triangle."""
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
    """Rotate the dofs in a quadrilateral."""
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
    """Reflect the dofs in a quadrilateral."""
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
