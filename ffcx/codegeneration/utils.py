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


def get_reflection(L, i):
    """Returns the bool that says whether or not an entity has been reflected."""
    if i[0] == 1:
        edge_reflections = L.Symbol("edge_reflections")
        return edge_reflections[i[1]]
    elif i[0] == 2:
        face_reflections = L.Symbol("face_reflections")
        return face_reflections[i[1]]
    else:
        return L.LiteralBool(False)


def get_vector_reflection_array(L, dof_reflection_entities, vname="reflected_dofs"):
    """Returns array containing true for dofs that require multiplying by -1."""
    _vnames_to_reflect[vname] = False
    reflect_dofs = []
    c_false = L.LiteralBool(False)
    for i in dof_reflection_entities:
        if i is None:
            # Dof does not need reflecting, so put false in array
            reflect_dofs.append(c_false)
        else:
            # Loop through entities that the direction of the dof depends on to
            # make a conditional
            ref = c_false
            for j in i:
                if ref == c_false:
                    # No condition has been added yet, so overwrite false
                    ref = get_reflection(L, j)
                else:
                    # This is not the first condition, so XOR
                    ref = L.Conditional(get_reflection(L, j), L.Not(ref), ref)
            reflect_dofs.append(ref)
            if ref != c_false:
                # Mark this space as needing reflections
                _vnames_to_reflect[vname] = True

    # If no dofs need reflecting, don't write any array
    if not _vnames_to_reflect[vname]:
        return []

    # If at least one vector dof needs reflecting
    return [L.ArrayDecl(
        "const bool", L.Symbol(vname), (len(reflect_dofs), ), values=reflect_dofs)]


def get_vector_reflection(L, idof, vname="reflected_dofs", tablename=None):
    """Gives 1 or -1 to scale a vector."""
    assert vname in _vnames_to_reflect
    if _vnames_to_reflect[vname]:
        # If at least one vector dof needs reflecting, return a conditional that gives -1
        # if the dof needs negating
        if tablename is None:
            return L.Conditional(L.Symbol(vname)[idof], 1, -1)
        return L.Conditional(L.Symbol(vname)[get_table_dofmap(L, tablename, idof)], 1, -1)
    else:
        # If no dofs need reflecting, return 1
        return 1


def get_table_dofmap_array(L, dofmap, pname):
    """If some dofs have been stripped from a table, return an array to map table entries
    to dofs."""
    for i, j in enumerate(dofmap):
        if i != j:
            # If a dof has been removed, write the data
            _table_dofmaps[pname] = L.Symbol(pname + "_dofmap")
            return [L.ArrayDecl(
                "const int", _table_dofmaps[pname], (len(dofmap), ), values=dofmap)]
    # If no dofs  have been removed, return nothing
    _table_dofmaps[pname] = None
    return []


def get_table_dofmap(L, pname, idof):
    """Returns the dof number corresponding to a row in the table."""
    assert pname in _table_dofmaps
    if _table_dofmaps[pname] is None:
        # If no dofs have been removed, then idof is the dof number
        return idof
    return _table_dofmaps[pname][idof]
