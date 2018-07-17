# -*- coding: utf-8 -*-
# Copyright (C) 2009-2018 Anders Logg, Martin Sandve Alnæs and Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
#
# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

import ffc.backends.ufc.dofmap_template as ufc_dofmap
from ffc.backends.ufc.utils import generate_return_new_switch


def tabulate_dofs(L, ir):
    # Input arguments
    entity_indices = L.Symbol("entity_indices")
    num_mesh_entities = L.Symbol("num_global_entities")

    # Output arguments
    dofs_variable = L.Symbol("dofs")

    ir = ir["tabulate_dofs"]
    if ir is None:
        # Special case for SpaceOfReals, ir returns None

        # FIXME: This is how the old code did in this case, is it correct?
        # Is there only 1 dof? I guess VectorElement(Real) handled
        # elsewhere?
        code = [L.Assign(dofs_variable[0], 0)]
        return L.StatementList(code)

    # Extract representation
    # (dofs_per_element, num_dofs_per_element, need_offset, fakes) = ir # names from original ffc code
    (subelement_dofs, num_dofs_per_subelement, need_offset, is_subelement_real) = ir

    # len(is_subelement_real) == len(subelement_dofs) == len(num_dofs_per_subelement)
    #    == number of combined (flattened) subelements?
    # is_subelement_real[subelement_number] is bool, is subelement Real?
    # num_dofs_per_subelement[subelement_number] is int, number of dofs for each (flattened) subelement
    # subelement_dofs[subelement_number] is list of entity_dofs for each (flattened) subelement
    # len(entity_dofs) == tdim+1
    # len(entity_dofs[cell_entity_dim]) == "num_cell_entities[dim]" (i.e. 3 vertices, 3 edges, 1 face for triangle)
    # len(entity_dofs[cell_entity_dim][entity_index[dim][k]]) == number of dofs for this entity
    # num = entity_dofs[dim]

    # entity_dofs =? entity_dofs[d][i][:] # dofs on entity (d,i)

    # Collect code pieces in list
    code = []

    # Declare offset if needed
    if need_offset:
        offset = L.Symbol("offset")
        code.append(L.VariableDecl("int64_t", offset, value=0))
    else:
        offset = 0

    # Generate code for each element
    subelement_offset = 0
    for (subelement_index, entity_dofs) in enumerate(subelement_dofs):

        # Handle is_subelement_real (Space of reals)
        if is_subelement_real[subelement_index]:
            assert num_dofs_per_subelement[subelement_index] == 1
            code.append(L.Assign(dofs_variable[subelement_offset], offset))
            if need_offset:
                code.append(L.AssignAdd(offset, 1))
            subelement_offset += 1
            continue

        # Generate code for each degree of freedom for each dimension
        for (cell_entity_dim, dofs_on_cell_entity) in enumerate(entity_dofs):
            num_dofs_per_mesh_entity = len(dofs_on_cell_entity[0])
            assert all(num_dofs_per_mesh_entity == len(dofs) for dofs in dofs_on_cell_entity)

            # Ignore if no dofs for this dimension
            if num_dofs_per_mesh_entity == 0:
                continue

            # For each cell entity of this dimension
            for (cell_entity_index, dofs) in enumerate(dofs_on_cell_entity):
                # dofs is a list of the local dofs that live on this cell entity

                # find offset for this particular mesh entity
                entity_offset = len(dofs) * entity_indices[cell_entity_dim, cell_entity_index]

                for (j, dof) in enumerate(dofs):
                    # dof is the local dof index on the subelement
                    # j is the local index of dof among the dofs
                    # on this particular cell/mesh entity
                    local_dof_index = subelement_offset + dof
                    global_dof_index = offset + entity_offset + j
                    code.append(L.Assign(dofs_variable[local_dof_index], global_dof_index))

            # Update offset corresponding to mesh entity:
            if need_offset:
                value = num_dofs_per_mesh_entity * num_mesh_entities[cell_entity_dim]
                code.append(L.AssignAdd(offset, value))

        subelement_offset += num_dofs_per_subelement[subelement_index]

    return L.StatementList(code)


def tabulate_facet_dofs(L, ir):
    all_facet_dofs = ir["tabulate_facet_dofs"]

    # Input arguments
    facet = L.Symbol("facet")
    dofs = L.Symbol("dofs")

    # For each facet, copy all_facet_dofs[facet][:] into output argument array dofs[:]
    cases = []
    for f, single_facet_dofs in enumerate(all_facet_dofs):
        assignments = [L.Assign(dofs[i], dof) for (i, dof) in enumerate(single_facet_dofs)]
        if assignments:
            cases.append((f, L.StatementList(assignments)))
    if cases:
        return L.Switch(facet, cases, autoscope=False)
    else:
        return L.NoOp()


def tabulate_dof_permutations(L, ir):
    edge_perms, facet_perms, cell = ir["dof_permutations"]
    tdim = cell.topological_dimension()

    perm = L.Symbol("perm")
    ndofs = L.Symbol("ndofs")  # FIXME: get ndofs from somewhere
    i = L.Symbol("i")
    code = [L.ForRange(i, 0, ndofs, body=[L.Assign(perm[i], i)])]

    if tdim == 1 or (len(edge_perms) == 0 and len(facet_perms) == 0):
        code += [L.Return()]
        return L.StatementList(code)

    global_indices = L.Symbol("global_indices")

    # Copied from ufc_geometry.h
    edge_vertices = {'triangle': ((1, 2), (0, 2), (0, 1)),
                     'tetrahedron': ((2, 3), (1, 3), (1, 2), (0, 3), (0, 2), (0, 1)),
                     'quadrilateral': ((0, 1), (2, 3), (0, 2), (1, 3)),
                     'hexahedron': ((0, 1), (2, 3), (4, 5), (6, 7), (0, 2), (1, 3),
                                    (4, 6), (5, 7), (0, 4), (1, 5), (2, 6), (3, 7))}
    facet_edge_vertices = {'tetrahedron': (((2, 3), (1, 3), (1, 2)),
                                           ((2, 3), (0, 3), (0, 2)),
                                           ((1, 3), (0, 3), (0, 1)),
                                           ((1, 2), (0, 2), (0, 1))),
                           'hexahedron': (((0, 1), (2, 3), (0, 2), (1, 3)),
                                          ((4, 5), (6, 7), (4, 6), (5, 7)),
                                          ((0, 1), (4, 5), (0, 4), (1, 5)),
                                          ((2, 3), (6, 7), (2, 6), (3, 7)),
                                          ((0, 2), (4, 6), (0, 4), (2, 6)),
                                          ((1, 3), (5, 7), (1, 5), (3, 7)))}

    # These are also fixed tables, but calculated here from data from ufc_geometry.h above
    facet_edges = {}
    facet_edges['tetrahedron'] = tuple(tuple(edge_vertices['tetrahedron'].index(edge) for edge in facet)
                                       for facet in facet_edge_vertices['tetrahedron'])
    facet_edges['hexahedron'] = tuple(tuple(edge_vertices['hexahedron'].index(edge) for edge in facet)
                                      for facet in facet_edge_vertices['hexahedron'])

    celltype = cell.cellname()
    num_edges = cell.num_edges()
    num_facets = cell.num_facets()
    edge_ordering = L.Symbol("edge_ordering")
    facet_ordering = L.Symbol("facet_ordering")

    code += [L.ArrayDecl("int", edge_ordering, [num_edges])]
    # Figure out the ordering of the global vertices on each edge
    # 0 = reference cell order, 1 = reversed
    for i in range(num_edges):
        code += [L.Assign(edge_ordering[i],
                          L.GT(global_indices[edge_vertices[celltype][i][0]],
                               global_indices[edge_vertices[celltype][i][1]]))]

    if celltype == 'tetrahedron':
        code += [L.ArrayDecl("int", facet_ordering, [num_facets])]
        for i in range(num_facets):
            code += [L.Assign(facet_ordering[i],
                              edge_ordering[facet_edges[celltype][i][0]]
                              + 2 * (edge_ordering[facet_edges[celltype][i][1]]
                                     + edge_ordering[facet_edges[celltype][i][2]]))]
    elif celltype == 'hexahedron':
        # Work out some permutations on quadrilateral facets of hexahedron
        # There are 8 possible permutations with Z-ordering
        # FIXME: more work to be done on this
        #
        f_edge_verts = facet_edge_vertices['hexahedron']
        cross_facet_order = L.Symbol("xf")
        t0 = L.Symbol('t0')
        t1 = L.Symbol('t1')
        code += [L.ArrayDecl("int", cross_facet_order, [2]),
                 L.VariableDecl("int", t0),
                 L.VariableDecl("int", t1)]
        for i in range(num_facets):
            code += [L.Assign(cross_facet_order[0], L.GT(global_indices[f_edge_verts[i][0][0]],
                                                         global_indices[f_edge_verts[i][1][1]])),
                     L.Assign(cross_facet_order[1], L.GT(global_indices[f_edge_verts[i][0][1]],
                                                         global_indices[f_edge_verts[i][1][0]]))]
            # Figure out which vertex has the lowest global index and then which neighbour is next
            tcases = [(0, [L.Assign(t1, cross_facet_order[1]),
                           L.If(edge_ordering[facet_edges[celltype][i][2]],
                                L.Assign(t1, 4 + cross_facet_order[0]))]),
                      (1, [L.Assign(t1, cross_facet_order[1]),
                           L.If(cross_facet_order[0],
                                L.Assign(t1, 6 + cross_facet_order[1]))]),
                      (2, [L.Assign(t1, 2 + cross_facet_order[0]),
                           L.If(cross_facet_order[1],
                                L.Assign(t1, 4 + cross_facet_order[0]))]),
                      (3, [L.Assign(t1, 2 + cross_facet_order[0]),
                           L.If(edge_ordering[facet_edges[celltype][i][3]],
                                L.Assign(t1, 6 + cross_facet_order[1]))])]

            code += [L.Assign(t0, 2 * edge_ordering[facet_edges[celltype][i][0]]
                              + edge_ordering[facet_edges[celltype][i][1]]),
                     L.Switch(t0, tcases),
                     L.VerbatimStatement('printf("t0=%d t1=%d\\n", t0, t1);')]

    # Sanity check that edge and face dofs are distinct - do this test in representation.py
    assert len(set(edge_perms.keys()).intersection(facet_perms.keys())) == 0

    # Change map indexing so it is listed by edge/facet. FIXME: go back and do this in representation.py instead

    edge_perms2 = [{} for i in range(num_edges)]
    for idx, p in edge_perms.items():
        edge_perms2[p[0]][idx] = p[1]

    facet_perms2 = [{} for i in range(num_facets)]
    for idx, p in facet_perms.items():
        facet_perms2[p[0]][idx] = p[1]

    # Mkae changes to the identity mapping where required for specific edges/facets

    for i, q in enumerate(edge_perms2):
        assignments = [L.Assign(perm[j], k) for j, k in q.items()]
        code += [L.If(edge_ordering[i], assignments)]

    for i, q in enumerate(facet_perms2):
        assignments = [(w, [L.Assign(perm[j], k[w]) for j, k in q.items()]) for w in range(6)]
        # FIXME - get 6 from somewhere (tet specific)
        code += [L.Switch(facet_ordering[i], assignments)]

    return L.StatementList(code)


def tabulate_entity_dofs(L, ir):
    entity_dofs, num_dofs_per_entity = ir["tabulate_entity_dofs"]

    # Output argument array
    dofs = L.Symbol("dofs")

    # Input arguments
    d = L.Symbol("d")
    i = L.Symbol("i")

    # TODO: Removed check for (d <= tdim + 1)
    tdim = len(num_dofs_per_entity) - 1

    # Generate cases for each dimension:
    all_cases = []
    for dim in range(tdim + 1):

        # Ignore if no entities for this dimension
        if num_dofs_per_entity[dim] == 0:
            continue

        # Generate cases for each mesh entity
        cases = []
        for entity in range(len(entity_dofs[dim])):
            casebody = []
            for (j, dof) in enumerate(entity_dofs[dim][entity]):
                casebody += [L.Assign(dofs[j], dof)]
            cases.append((entity, L.StatementList(casebody)))

        # Generate inner switch
        # TODO: Removed check for (i <= num_entities-1)
        inner_switch = L.Switch(i, cases, autoscope=False)
        all_cases.append((dim, inner_switch))

    if all_cases:
        return L.Switch(d, all_cases, autoscope=False)
    else:
        return L.NoOp()


def tabulate_entity_closure_dofs(L, ir):
    # Extract variables from ir
    entity_closure_dofs, entity_dofs, num_dofs_per_entity = \
        ir["tabulate_entity_closure_dofs"]

    # Output argument array
    dofs = L.Symbol("dofs")

    # Input arguments
    d = L.Symbol("d")
    i = L.Symbol("i")

    # TODO: Removed check for (d <= tdim + 1)
    tdim = len(num_dofs_per_entity) - 1

    # Generate cases for each dimension:
    all_cases = []
    for dim in range(tdim + 1):
        num_entities = len(entity_dofs[dim])

        # Generate cases for each mesh entity
        cases = []
        for entity in range(num_entities):
            casebody = []
            for (j, dof) in enumerate(entity_closure_dofs[(dim, entity)]):
                casebody += [L.Assign(dofs[j], dof)]
            cases.append((entity, L.StatementList(casebody)))

        # Generate inner switch unless empty
        if cases:
            # TODO: Removed check for (i <= num_entities-1)
            inner_switch = L.Switch(i, cases, autoscope=False)
            all_cases.append((dim, inner_switch))

    if all_cases:
        return L.Switch(d, all_cases, autoscope=False)
    else:
        return L.NoOp()


def create_sub_dofmap(L, ir):
    classnames = ir["create_sub_dofmap"]
    return generate_return_new_switch(L, "i", classnames, factory=True)


def sub_dofmap_declaration(L, ir):
    classnames = set(ir["create_sub_dofmap"])
    code = ""
    for name in classnames:
        code += "ufc_dofmap* create_{name}(void);\n".format(name=name)
    return code


def ufc_dofmap_generator(ir, parameters):
    """Generate UFC code for a dofmap"""

    d = {}

    # Attributes
    d["factory_name"] = ir["classname"]
    d["signature"] = "\"{}\"".format(ir["signature"])
    d["num_global_support_dofs"] = ir["num_global_support_dofs"]
    d["num_element_support_dofs"] = ir["num_element_support_dofs"]
    d["num_element_dofs"] = ir["num_element_dofs"]
    d["num_facet_dofs"] = ir["num_facet_dofs"]
    d["num_sub_dofmaps"] = ir["num_sub_dofmaps"]
    d["num_entity_dofs"] = ir["num_entity_dofs"] + [0, 0, 0, 0]
    d["num_entity_closure_dofs"] = ir["num_entity_closure_dofs"] + [0, 0, 0, 0]

    import ffc.uflacs.language.cnodes as L

    # Functions
    d["tabulate_dofs"] = tabulate_dofs(L, ir)
    d["tabulate_dof_permutations"] = tabulate_dof_permutations(L, ir)
    d["tabulate_facet_dofs"] = tabulate_facet_dofs(L, ir)
    d["tabulate_entity_dofs"] = tabulate_entity_dofs(L, ir)
    d["tabulate_entity_closure_dofs"] = tabulate_entity_closure_dofs(L, ir)
    d["sub_dofmap_declaration"] = sub_dofmap_declaration(L, ir)
    d["create_sub_dofmap"] = create_sub_dofmap(L, ir)

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fields = [fname for _, fname, _, _ in Formatter().parse(ufc_dofmap.factory) if fname]
    # Remove square brackets from any field names
    fields = [f.split("[")[0] for f in fields]
    assert set(fields) == set(d.keys()), "Mismatch between keys in template and in formattting dict."

    # Format implementation code
    implementation = ufc_dofmap.factory.format_map(d)

    # Format declaration
    declaration = ufc_dofmap.declaration.format(factory_name=ir["classname"])

    return declaration, implementation
