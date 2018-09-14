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


def tabulate_dof_permutations(L, ir):
    """Generate code for a permutation vector for the dofmap, depending on the cell orientation, as
    determined by its global vertex indices. For triangles and quads, this requires reversing the
    ordering of some edge dofs; for 3D cells, the facet dofs may also
    be permuted by rotation or reflection in the plane. """

    edge_perms, facet_perms, cell, cell_topology = ir["dof_permutations"]
    ndofs = ir["num_element_support_dofs"] + ir["num_global_support_dofs"]

    perm = L.Symbol("perm")
    i = L.Symbol("i")
    # Fill up an identity permutation, for any dofs which are not permuted.
    code = [L.ForRange(i, 0, ndofs, body=[L.Assign(perm[i], i)])]

    tdim = cell.topological_dimension()
    if tdim == 1 or (len(edge_perms) == 0 and len(facet_perms) == 0):
        code += [L.Return()]
        return L.StatementList(code)

    celltype = cell.cellname()
    num_edges = cell.num_edges()
    num_facets = cell.num_facets()

    global_indices = L.Symbol("global_indices")
    edge_ordering = L.Symbol("edge_ordering")

    code += [L.ArrayDecl("int", edge_ordering, [num_edges])]
    # Calculate the edge order for each edge, may also be needed later for facet orientation
    edge_vertices = cell_topology['edge_vertices']
    for i in range(num_edges):
        # Figure out the ordering of the global vertices on each edge
        # 0 = reference cell order, 1 = reversed
        code += [
            L.Assign(edge_ordering[i],
                     L.GT(global_indices[edge_vertices[i][0]], global_indices[edge_vertices[i][1]]))
        ]

    # Make changes to the identity mapping where required for specific edges
    for i, q in enumerate(edge_perms):
        assignments = [L.Assign(perm[j], k) for j, k in q.items()]
        code += [L.If(edge_ordering[i], assignments)]

    if celltype == 'tetrahedron' and len(facet_perms) > 0:
        facet_edges = cell_topology['facet_edges']
        facet_ordering = L.Symbol("facet_ordering")
        code += [L.VariableDecl("int", facet_ordering, 0)]
        # Six possible orientations for each triangular facet
        assert len(facet_perms) == num_facets
        for i, q in enumerate(facet_perms):
            code += [
                L.Comment("DOF reordering for facet %d" % i),
                L.Assign(facet_ordering, edge_ordering[facet_edges[i][0]] + 2 *
                         (edge_ordering[facet_edges[i][1]] + edge_ordering[facet_edges[i][2]]))
            ]
            # Make changes to the identity mapping where required for specific facets
            cases = []
            for w in range(6):
                assignments = []
                for j, k in q.items():
                    if j != k[w]:  # Ignore cases where there is no permutation
                        assignments += [L.Assign(perm[j], k[w])]
                cases += [(w, assignments)]
            code += [L.Switch(facet_ordering, cases)]

    elif celltype == 'hexahedronxxx':  # FIXME - disabled for now
        # Work out some permutations on quadrilateral facets of hexahedron
        # There are 8 possible permutations with Z-ordering
        # FIXME: more work to be done on this
        #
        f_edge_verts = cell_topology['facet_edge_vertices']
        facet_edges = cell_topology['facet_edges']
        cross_facet_order = L.Symbol("xf")
        t0 = L.Symbol('t0')
        t1 = L.Symbol('t1')
        code += [
            L.ArrayDecl("int", cross_facet_order, [2]),
            L.VariableDecl("int", t0),
            L.VariableDecl("int", t1)
        ]
        for i in range(num_facets):
            code += [
                L.Assign(cross_facet_order[0],
                         L.GT(global_indices[f_edge_verts[i][0][0]],
                              global_indices[f_edge_verts[i][1][1]])),
                L.Assign(cross_facet_order[1],
                         L.GT(global_indices[f_edge_verts[i][0][1]],
                              global_indices[f_edge_verts[i][1][0]]))
            ]
            # Figure out which vertex has the lowest global index and then which neighbour is next
            tcases = [(0, [
                L.Assign(t1, cross_facet_order[1]),
                L.If(edge_ordering[facet_edges[i][2]], L.Assign(t1, 4 + cross_facet_order[0]))
            ]), (1, [
                L.Assign(t1, cross_facet_order[1]),
                L.If(cross_facet_order[0], L.Assign(t1, 6 + cross_facet_order[1]))
            ]), (2, [
                L.Assign(t1, 2 + cross_facet_order[0]),
                L.If(cross_facet_order[1], L.Assign(t1, 4 + cross_facet_order[0]))
            ]), (3, [
                L.Assign(t1, 2 + cross_facet_order[0]),
                L.If(edge_ordering[facet_edges[i][3]], L.Assign(t1, 6 + cross_facet_order[1]))
            ])]

            code += [
                L.Assign(t0,
                         2 * edge_ordering[facet_edges[i][0]] + edge_ordering[facet_edges[i][1]]),
                L.Switch(t0, tcases),
                L.VerbatimStatement('printf("t0=%d t1=%d\\n", t0, t1);')
            ]

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
    d["num_sub_dofmaps"] = ir["num_sub_dofmaps"]
    d["num_entity_dofs"] = ir["num_entity_dofs"] + [0, 0, 0, 0]
    d["num_entity_closure_dofs"] = ir["num_entity_closure_dofs"] + [0, 0, 0, 0]

    import ffc.language.cnodes as L

    # Functions
    d["tabulate_dof_permutations"] = tabulate_dof_permutations(L, ir)
    d["tabulate_entity_dofs"] = tabulate_entity_dofs(L, ir)
    d["tabulate_entity_closure_dofs"] = tabulate_entity_closure_dofs(L, ir)
    d["sub_dofmap_declaration"] = sub_dofmap_declaration(L, ir)
    d["create_sub_dofmap"] = create_sub_dofmap(L, ir)

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fields = [fname for _, fname, _, _ in Formatter().parse(ufc_dofmap.factory) if fname]
    # Remove square brackets from any field names
    fields = [f.split("[")[0] for f in fields]
    assert set(fields) == set(
        d.keys()), "Mismatch between keys in template and in formattting dict."

    # Format implementation code
    implementation = ufc_dofmap.factory.format_map(d)

    # Format declaration
    declaration = ufc_dofmap.declaration.format(factory_name=ir["classname"])

    return declaration, implementation
