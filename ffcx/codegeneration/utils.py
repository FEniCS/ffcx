# Copyright (C) 2015-2017 Martin Sandve Alnæs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# TODO: Move these to ffcx.language utils?


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


_vnames_to_reflect = {}
_table_dofmaps = {}


def get_vector_reflection_array(L, dof_types, space_dimension, entity_dofs, vname="reflected_dofs"):
    # List of dof types that require multiplying by -1 if their entity has been reflected
    # TODO: check that these are all vector and that no other types are vector
    vector_types = ["PointScaledNormalEval", "ComponentPointEval", "PointEdgeTangent",
                    # "PointFaceTangent",
                    "PointScaledNormalEval", "PointNormalEval",
                    "IntegralMoment", "FrobeniusIntegralMoment"
                    ]
    _vnames_to_reflect[vname] = False
    for v in vector_types:
        if v in dof_types:
            _vnames_to_reflect[vname] = True
            break

    if not _vnames_to_reflect[vname]:
        return []

    # For each dof that needs reflection this will contain the edge of face reflection that the dof is associated
    # with, or False if the dof is on a Point or inside a volume.
    # Entities are numbered: Point0, Point1, ..., Edge0, Edge1, ..., Face0, [Face1, ..., Volume]
    # If no reflection needed, this will be -1
    reflect_dofs = [False for i in range(space_dimension)]
    dof_n = 0
    face_reflections = L.Symbol("face_reflections")
    edge_reflections = L.Symbol("edge_reflections")

    # Run through the entities and mark vector dofs on each entity
    for dim, e_dofs in entity_dofs.items():
        for n, dofs in e_dofs.items():
            if dim == 1 or dim == 2:
                for d in dofs:
                    t = dof_types[d]
                    if t in vector_types:
                        if dim == 1:
                            reflect_dofs[dof_n] = edge_reflections[n]
                        if dim == 2:
                            reflect_dofs[dof_n] = face_reflections[n]
                    dof_n += 1
            else:
                dof_n += len(dofs)

    # If at least one vector dof needs reflecting
    return [L.ArrayDecl(
        "const bool", L.Symbol(vname), (space_dimension, ), values=reflect_dofs)]


def get_vector_reflection(L, idof, vname="reflected_dofs", tablename=None):
    assert vname in _vnames_to_reflect
    # If at least one vector dof needs reflecting
    if _vnames_to_reflect[vname]:
        if tablename is None:
            return L.Conditional(L.Symbol(vname)[idof], 1, -1)
        return L.Conditional(L.Symbol(vname)[get_table_dofmap(L, tablename, idof)], 1, -1)
    # If no dofs need reflecting
    else:
        return 1


def get_table_dofmap_array(L, dofmap, pname):
    for i, j in enumerate(dofmap):
        if i != j:
            _table_dofmaps[pname] = L.Symbol(pname + "_dofmap")
            return [L.ArrayDecl(
                "const int", _table_dofmaps[pname], (len(dofmap), ), values=dofmap)]
    _table_dofmaps[pname] = None
    return []


def get_table_dofmap(L, pname, idof):
    assert pname in _table_dofmaps
    if _table_dofmaps[pname] is None:
        return idof
    return _table_dofmaps[pname][idof]
