# Copyright (C) 2015-2017 Martin Sandve Alnæs
# Modified by Matthew Scroggs, 2020-2021
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# TODO: Move these to ffcx.language utils?
import numpy
index_type = "int"


def generate_return_new(L, classname):
    return L.Return(L.Call("create_" + classname))


def generate_return_new_switch(L, i, classnames, args=None):

    if isinstance(i, str):
        i = L.Symbol(i)

    def create(classname):
        return L.Call("create_" + classname)

    default = L.Return(L.Null())
    if classnames:
        cases = []
        if args is None:
            args = list(range(len(classnames)))
        for j, classname in zip(args, classnames):
            if classname:
                cases.append((j, L.Return(create(classname))))
        return L.Switch(i, cases, default=default)
    else:
        return default


def generate_return_literal_switch(L,
                                   i,
                                   values,
                                   default,
                                   literal_type,
                                   typename=None):
    # TODO: UFC functions of this type could be replaced with return vector<T>{values}.

    if isinstance(i, str):
        i = L.Symbol(i)
    return_default = L.Return(literal_type(default))

    if values and typename is not None:
        # Store values in static table and return from there
        V = L.Symbol("return_values")
        decl = L.ArrayDecl("static const %s" % typename, V, len(values),
                           [literal_type(k) for k in values])
        return L.StatementList(
            [decl,
             L.If(L.GE(i, len(values)), return_default),
             L.Return(V[i])])
    elif values:
        # Need typename to create static array, fallback to switch
        cases = [(j, L.Return(literal_type(k))) for j, k in enumerate(values)]
        return L.Switch(i, cases, default=return_default)
    else:
        # No values, just return default
        return return_default


def generate_return_int_switch(L, i, values, default):
    return generate_return_literal_switch(L, i, values, default, L.LiteralInt,
                                          "int")


def make_perm_data(L, base_perms, cell_shape, rotations_first=False, reversed_rotations=False):
    if cell_shape == "interval":
        entities = {}
    elif cell_shape == "triangle":
        entities = {1: 3}
    elif cell_shape == "quadrilateral":
        entities = {1: 4}
    elif cell_shape == "tetrahedron":
        entities = {1: 6, 2: 4}
        face_rotation_order = 3
    elif cell_shape == "hexahedron":
        entities = {1: 12, 2: 6}
        face_rotation_order = 4
    else:
        raise NotImplementedError

    perm_n = 0
    perm_data = []
    if 1 in entities:
        for edge in range(entities[1]):
            perm_data.append((
                entity_reflection(L, (1, edge), cell_shape),
                None,
                base_perms[perm_n]
            ))
            perm_n += 1
    if 2 in entities:
        for face in range(entities[2]):
            reflection = (
                entity_reflection(L, (2, face), cell_shape),
                None,
                base_perms[perm_n + 1]
            )
            if not rotations_first:
                perm_data.append(reflection)
            for rot in range(1, face_rotation_order):
                if reversed_rotations:
                    rot = face_rotation_order - rot
                perm_data.append((
                    entity_rotations(L, (2, face), cell_shape),
                    rot,
                    numpy.linalg.matrix_power(base_perms[perm_n], rot)
                ))
            if rotations_first:
                perm_data.append(reflection)
            perm_n += 2

    assert perm_n == len(base_perms)

    return perm_data


def apply_permutations_to_data(L, base_permutations, cell_shape, data,
                               rotations_first=False, reversed_rotations=False,
                               indices=lambda dof: dof, ranges=None, dtype="double"):
    perm_data = make_perm_data(
        L, base_permutations, cell_shape,
        rotations_first=rotations_first, reversed_rotations=reversed_rotations)

    # Apply entity permutations
    apply_permutations = []
    temporary_variables = 0
    for entity_perm, value, perm in perm_data:
        body = []

        # Use temporary variables t0, t1, ... to store current data
        temps = {}
        for index, row in enumerate(perm):
            if not numpy.allclose(row, [1 if i == index else 0 for i, j in enumerate(row)]):
                for dof, w in enumerate(row):
                    if not numpy.isclose(w, 0) and dof not in temps:
                        temps[dof] = L.Symbol("t" + str(len(temps)))
                body.append(L.Assign(data[indices(index)],
                                     sum(temps[dof] if numpy.isclose(w, 1) else w * temps[dof]
                                         for dof, w in enumerate(row) if not numpy.isclose(w, 0))))
        temporary_variables = max(temporary_variables, len(temps))

        # If no changes would be made, continue to next entity
        if len(body) == 0:
            continue

        if value is None:
            condition = entity_perm
        else:
            condition = L.EQ(entity_perm, value)

        body = [L.Assign(t, data[indices(dof)]) for dof, t in temps.items()] + body
        if ranges is None:
            apply_permutations.append(L.If(condition, body))
        else:
            apply_permutations.append(L.If(condition,
                                           L.ForRanges(*ranges, index_type=index_type, body=body)))

    if len(apply_permutations) > 0:
        apply_permutations = [L.VariableDecl(dtype, L.Symbol("t" + str(i)), 0)
                              for i in range(temporary_variables)] + apply_permutations

    return apply_permutations


def entity_reflection(L, i, cell_shape):
    """Returns the bool that says whether or not an entity has been reflected."""
    cell_info = L.Symbol("cell_permutation")
    if cell_shape in ["triangle", "quadrilateral"]:
        num_faces = 0
        face_bitsize = 1
        assert i[0] == 1
    if cell_shape == "tetrahedron":
        num_faces = 4
        face_bitsize = 3
    if cell_shape == "hexahedron":
        num_faces = 6
        face_bitsize = 3
    if i[0] == 1:
        return L.NE(L.BitwiseAnd(cell_info, L.BitShiftL(1, face_bitsize * num_faces + i[1])), 0)
    elif i[0] == 2:
        return L.NE(L.BitwiseAnd(cell_info, L.BitShiftL(1, face_bitsize * i[1])), 0)
    return L.LiteralBool(False)


def entity_rotations(L, i, cell_shape):
    """Returns number of times an entity has been rotated."""
    cell_info = L.Symbol("cell_permutation")
    assert cell_shape in ["tetrahedron", "hexahedron"]
    assert i[0] == 2
    return L.BitwiseAnd(L.BitShiftR(cell_info, 3 * i[1] + 1), 3)
