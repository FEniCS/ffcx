# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FFCx/UFC specific variable access."""

import logging
import warnings

import basix.ufl
import ffcx.codegeneration.lnodes as L
import ufl
from ffcx.ir.analysis.modified_terminals import ModifiedTerminal
from ffcx.ir.elementtables import UniqueTableReferenceT
from ffcx.ir.representationutils import QuadratureRule
from ffcx.codegeneration.backend import ReprManager

logger = logging.getLogger("ffcx")


def get(repr: ReprManager,
        mt: ModifiedTerminal,
        tabledata: UniqueTableReferenceT,
        quadrature_rule: QuadratureRule):

    # Call appropriate handler, depending on the type of terminal
    terminal = mt.terminal
    ttype = type(terminal)

    # Look for parent class of ttype or direct handler
    while ttype not in call_lookup and ttype.__bases__:  # noqa: E721
        ttype = ttype.__bases__[0]

    # Get the handler from the lookup, or None if not found
    handler = call_lookup.get(ttype)

    if handler is None:
        raise NotImplementedError(f"No handler for terminal type: {ttype}")

    # Call the handler
    return handler(repr, mt, tabledata, quadrature_rule)


def coefficient(repr: ReprManager,
                mt: ModifiedTerminal,
                tabledata: UniqueTableReferenceT,
                quadrature_rule: QuadratureRule):

    ttype = tabledata.ttype
    assert ttype != "zeros"

    num_dofs = tabledata.values.shape[3]
    begin = tabledata.offset
    end = begin + tabledata.block_size * (num_dofs - 1) + 1

    if ttype == "ones" and (end - begin) == 1:
        # f = 1.0 * f_{begin}, just return direct reference to dof
        # array at dof begin (if mt is restricted, begin contains
        # cell offset)
        return repr.symbols.coefficient_dof_access(mt.terminal, begin)
    else:
        # Return symbol, see definitions for computation
        return repr.symbols.coefficient_value(mt)


def constant(repr: ReprManager,
             mt: ModifiedTerminal,
             tabledata: UniqueTableReferenceT,
             quadrature_rule: QuadratureRule):
    """Access to a constant is handled trivially, directly through constants symbol."""
    return repr.symbols.constant_index_access(mt.terminal, mt.flat_component)


def spatial_coordinate(repr: ReprManager,
                       mt: ModifiedTerminal,
                       tabledata: UniqueTableReferenceT,
                       quadrature_rule: QuadratureRule):
    if mt.global_derivatives:
        raise RuntimeError("Not expecting global derivatives of SpatialCoordinate.")
    if mt.averaged is not None:
        raise RuntimeError("Not expecting average of SpatialCoordinates.")

    integral_type = repr.ir.integral_type
    if integral_type in ufl.custom_integral_types:
        if mt.local_derivatives:
            raise RuntimeError("FIXME: Jacobian in custom integrals is not implemented.")

        # Access predefined quadrature points table
        x = repr.symbols.custom_points_table
        iq = repr.symbols.quadrature_loop_index
        gdim, = mt.terminal.ufl_shape
        if gdim == 1:
            index = iq
        else:
            index = iq * gdim + mt.flat_component
        return x[index]
    elif integral_type == "expression":
        # Physical coordinates are computed by code generated in
        # definitions
        return repr.symbols.x_component(mt)
    else:
        # Physical coordinates are computed by code generated in
        # definitions
        return repr.symbols.x_component(mt)


def cell_coordinate(repr: ReprManager,
                    mt: ModifiedTerminal,
                    tabledata: UniqueTableReferenceT,
                    quadrature_rule: QuadratureRule):
    if mt.global_derivatives:
        raise RuntimeError("Not expecting derivatives of CellCoordinate.")
    if mt.local_derivatives:
        raise RuntimeError("Not expecting derivatives of CellCoordinate.")
    if mt.averaged is not None:
        raise RuntimeError("Not expecting average of CellCoordinate.")

    integral_type = repr.ir.integral_type
    if integral_type == "cell" and not mt.restriction:
        # Access predefined quadrature points table
        X = repr.symbols.points_table(quadrature_rule)
        tdim, = mt.terminal.ufl_shape
        iq = repr.symbols.quadrature_loop_index()
        if quadrature_rule.points.size == 1:
            index = mt.flat_component
        elif tdim == 1:
            index = iq
        else:
            index = iq * tdim + mt.flat_component
        return X[index]
    else:
        # X should be computed from x or Xf symbolically instead of
        # getting here
        raise RuntimeError("Expecting reference cell coordinate to be symbolically rewritten.")


def facet_coordinate(repr: ReprManager,
                     mt: ModifiedTerminal,
                     tabledata: UniqueTableReferenceT,
                     quadrature_rule: QuadratureRule):
    if mt.global_derivatives:
        raise RuntimeError("Not expecting derivatives of FacetCoordinate.")
    if mt.local_derivatives:
        raise RuntimeError("Not expecting derivatives of FacetCoordinate.")
    if mt.averaged is not None:
        raise RuntimeError("Not expecting average of FacetCoordinate.")
    if mt.restriction:
        raise RuntimeError("Not expecting restriction of FacetCoordinate.")

    integral_type = repr.ir.integral_type
    if integral_type in ("interior_facet", "exterior_facet"):
        tdim, = mt.terminal.ufl_shape
        if tdim == 0:
            raise RuntimeError("Vertices have no facet coordinates.")
        elif tdim == 1:
            warnings.warn(
                "Vertex coordinate is always 0, should get rid of this in ufl geometry lowering."
            )
            return L.LiteralFloat(0.0)
        Xf = repr.symbols.points_table(quadrature_rule)
        iq = repr.symbols.quadrature_loop_index()
        assert 0 <= mt.flat_component < (tdim - 1)
        if quadrature_rule.points.size == 1:
            index = mt.flat_component
        elif tdim == 2:
            index = iq
        else:
            index = iq * (tdim - 1) + mt.flat_component
        return Xf[index]
    else:
        # Xf should be computed from X or x symbolically instead of
        # getting here
        raise RuntimeError("Expecting reference facet coordinate to be symbolically rewritten.")


def jacobian(repr: ReprManager,
             mt: ModifiedTerminal,
             tabledata: UniqueTableReferenceT,
             quadrature_rule: QuadratureRule):
    if mt.averaged is not None:
        raise RuntimeError("Not expecting average of Jacobian.")
    return repr.symbols.J_component(mt)


def reference_cell_volume(repr: ReprManager,
                          mt: ModifiedTerminal,
                          tabledata: UniqueTableReferenceT,
                          quadrature_rule: QuadratureRule):
    cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
    if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
        return L.Symbol(f"{cellname}_reference_cell_volume", dtype=L.DataType.REAL)
    else:
        raise RuntimeError(f"Unhandled cell types {cellname}.")


def reference_facet_volume(repr: ReprManager,
                           mt: ModifiedTerminal,
                           tabledata: UniqueTableReferenceT,
                           quadrature_rule: QuadratureRule):

    cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
    if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
        return L.Symbol(f"{cellname}_reference_facet_volume", dtype=L.DataType.REAL)
    else:
        raise RuntimeError(f"Unhandled cell types {cellname}.")


def reference_normal(repr: ReprManager,
                     mt: ModifiedTerminal,
                     tabledata: UniqueTableReferenceT,
                     quadrature_rule: QuadratureRule):
    cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
    if cellname in ("interval", "triangle", "tetrahedron", "quadrilateral", "hexahedron"):
        table = L.Symbol(f"{cellname}_reference_facet_normals", dtype=L.DataType.REAL)
        facet = repr.symbols.entity("facet", mt.restriction)
        return table[facet][mt.component[0]]
    else:
        raise RuntimeError(f"Unhandled cell types {cellname}.")


def cell_facet_jacobian(repr: ReprManager,
                        mt: ModifiedTerminal,
                        tabledata: UniqueTableReferenceT,
                        quadrature_rule: QuadratureRule):
    cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
    if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
        table = L.Symbol(f"{cellname}_reference_facet_jacobian", dtype=L.DataType.REAL)
        facet = repr.symbols.entity("facet", mt.restriction)
        return table[facet][mt.component[0]][mt.component[1]]
    elif cellname == "interval":
        raise RuntimeError("The reference facet jacobian doesn't make sense for interval cell.")
    else:
        raise RuntimeError(f"Unhandled cell types {cellname}.")


def reference_cell_edge_vectors(repr: ReprManager,
                                mt: ModifiedTerminal,
                                tabledata: UniqueTableReferenceT,
                                quadrature_rule: QuadratureRule):
    cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
    if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
        table = L.Symbol(f"{cellname}_reference_edge_vectors", dtype=L.DataType.REAL)
        return table[mt.component[0]][mt.component[1]]
    elif cellname == "interval":
        raise RuntimeError("The reference cell edge vectors doesn't make sense for interval cell.")
    else:
        raise RuntimeError(f"Unhandled cell types {cellname}.")


def reference_facet_edge_vectors(repr: ReprManager,
                                 mt: ModifiedTerminal,
                                 tabledata: UniqueTableReferenceT,
                                 quadrature_rule: QuadratureRule):
    cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
    if cellname in ("tetrahedron", "hexahedron"):
        table = L.Symbol(f"{cellname}_reference_edge_vectors", dtype=L.DataType.REAL)
        facet = repr.symbols.entity("facet", mt.restriction)
        return table[facet][mt.component[0]][mt.component[1]]
    elif cellname in ("interval", "triangle", "quadrilateral"):
        raise RuntimeError(
            "The reference cell facet edge vectors doesn't make sense for interval or triangle cell."
        )
    else:
        raise RuntimeError(f"Unhandled cell types {cellname}.")


def facet_orientation(repr: ReprManager,
                      mt: ModifiedTerminal,
                      tabledata: UniqueTableReferenceT,
                      quadrature_rule: QuadratureRule):
    cellname = ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
    if cellname not in ("interval", "triangle", "tetrahedron"):
        raise RuntimeError(f"Unhandled cell types {cellname}.")

    table = L.Symbol(f"{cellname}_facet_orientations", dtype=L.DataType.INT)
    facet = repr.symbols.entity("facet", mt.restriction)
    return table[facet]


def cell_vertices(repr: ReprManager,
                  mt: ModifiedTerminal,
                  tabledata: UniqueTableReferenceT,
                  quadrature_rule: QuadratureRule):
    # Get properties of domain
    domain = ufl.domain.extract_unique_domain(mt.terminal)
    gdim = domain.geometric_dimension()
    coordinate_element = domain.ufl_coordinate_element()

    # Get dimension and dofmap of scalar element
    assert isinstance(coordinate_element, basix.ufl._BlockedElement)
    assert coordinate_element.value_shape == (gdim, )
    ufl_scalar_element, = set(coordinate_element.sub_elements)
    scalar_element = ufl_scalar_element
    assert scalar_element.value_size == 1 and scalar_element.block_size == 1

    vertex_scalar_dofs = scalar_element.entity_dofs[0]
    num_scalar_dofs = scalar_element.dim

    # Get dof and component
    dof, = vertex_scalar_dofs[mt.component[0]]
    component = mt.component[1]

    expr = repr.symbols.domain_dof_access(dof, component, gdim, num_scalar_dofs, mt.restriction)
    return expr


def cell_edge_vectors(repr: ReprManager,
                      mt: ModifiedTerminal,
                      tabledata: UniqueTableReferenceT,
                      quadrature_rule: QuadratureRule):
    # Get properties of domain
    domain = ufl.domain.extract_unique_domain(mt.terminal)
    cellname = domain.ufl_cell().cellname()
    gdim = domain.geometric_dimension()
    coordinate_element = domain.ufl_coordinate_element()

    if cellname in ("triangle", "tetrahedron", "quadrilateral", "hexahedron"):
        pass
    elif cellname == "interval":
        raise RuntimeError("The physical cell edge vectors doesn't make sense for interval cell.")
    else:
        raise RuntimeError(f"Unhandled cell types {cellname}.")

    # Get dimension and dofmap of scalar element
    assert isinstance(coordinate_element, basix.ufl._BlockedElement)
    assert coordinate_element.value_shape == (gdim, )
    ufl_scalar_element, = set(coordinate_element.sub_elements)
    scalar_element = ufl_scalar_element
    assert scalar_element.value_size == 1 and scalar_element.block_size == 1

    vertex_scalar_dofs = scalar_element.entity_dofs[0]
    num_scalar_dofs = scalar_element.dim

    # Get edge vertices
    edge = mt.component[0]
    vertex0, vertex1 = scalar_element.reference_topology[1][edge]

    # Get dofs and component
    dof0, = vertex_scalar_dofs[vertex0]
    dof1, = vertex_scalar_dofs[vertex1]
    component = mt.component[1]

    return repr.symbols.domain_dof_access(
        dof0, component, gdim, num_scalar_dofs, mt.restriction
    ) - repr.symbols.domain_dof_access(
        dof1, component, gdim, num_scalar_dofs, mt.restriction
    )


def facet_edge_vectors(repr: ReprManager,
                       mt: ModifiedTerminal,
                       tabledata: UniqueTableReferenceT,
                       quadrature_rule: QuadratureRule):

    # Get properties of domain
    domain = ufl.domain.extract_unique_domain(mt.terminal)
    cellname = domain.ufl_cell().cellname()
    gdim = domain.geometric_dimension()
    coordinate_element = domain.ufl_coordinate_element()

    if cellname in ("tetrahedron", "hexahedron"):
        pass
    elif cellname in ("interval", "triangle", "quadrilateral"):
        raise RuntimeError(
            f"The physical facet edge vectors doesn't make sense for {cellname} cell.")
    else:
        raise RuntimeError(f"Unhandled cell types {cellname}.")

    # Get dimension and dofmap of scalar element
    assert isinstance(coordinate_element, basix.ufl._BlockedElement)
    assert coordinate_element.value_shape == (gdim, )
    ufl_scalar_element, = set(coordinate_element.sub_elements)
    scalar_element = ufl_scalar_element
    assert scalar_element.value_size == 1 and scalar_element.block_size == 1

    scalar_element = ufl_scalar_element
    num_scalar_dofs = scalar_element.dim

    # Get edge vertices
    facet = repr.symbols.entity("facet", mt.restriction)
    facet_edge = mt.component[0]
    facet_edge_vertices = L.Symbol(f"{cellname}_facet_edge_vertices", dtype=L.DataType.INT)
    vertex0 = facet_edge_vertices[facet][facet_edge][0]
    vertex1 = facet_edge_vertices[facet][facet_edge][1]

    # Get dofs and component
    component = mt.component[1]
    assert coordinate_element.embedded_superdegree == 1, "Assuming degree 1 element"
    dof0 = vertex0
    dof1 = vertex1
    expr = (
        repr.symbols.domain_dof_access(dof0, component, gdim, num_scalar_dofs, mt.restriction)
        - repr.symbols.domain_dof_access(dof1, component, gdim, num_scalar_dofs, mt.restriction))

    return expr


def _pass(repr: ReprManager,
          mt: ModifiedTerminal,
          tabledata: UniqueTableReferenceT,
          quadrature_rule: QuadratureRule):
    """Return one."""
    return 1


def table_access(repr, tabledata: UniqueTableReferenceT, entitytype: str, restriction: str,
                 quadrature_index: L.MultiIndex, dof_index: L.MultiIndex):
    """
    Access element table for given entity, quadrature point, and dof index.

    Args:
        tabledata: Table data object
        entitytype: Entity type ("cell", "facet", "vertex")
        restriction: Restriction ("+", "-")
        quadrature_index: Quadrature index
        dof_index: Dof index
    """
    entity = repr.symbols.entity(entitytype, restriction)
    iq_global_index = quadrature_index.global_index
    ic_global_index = dof_index.global_index
    qp = 0  # quadrature permutation

    symbols = []
    if tabledata.is_uniform:
        entity = L.LiteralInt(0)

    if tabledata.is_piecewise:
        iq_global_index = L.LiteralInt(0)

    # FIXME: Hopefully tabledata is not permuted when applying sum
    # factorization
    if tabledata.is_permuted:
        qp = repr.symbols.quadrature_permutation[0]
        if restriction == "-":
            qp = repr.symbols.quadrature_permutation[1]

    if dof_index.dim == 1 and quadrature_index.dim == 1:
        symbols += [L.Symbol(tabledata.name, dtype=L.DataType.REAL)]
        return repr.symbols.element_tables[tabledata.name][qp][entity][iq_global_index][ic_global_index], symbols
    else:
        FE = []
        for i in range(dof_index.dim):
            factor = tabledata.tensor_factors[i]
            iq_i = quadrature_index.local_index(i)
            ic_i = dof_index.local_index(i)
            table = repr.symbols.element_tables[factor.name][qp][entity][iq_i][ic_i]
            symbols += [repr.symbols.element_tables[factor.name]]
            FE.append(table)
        return L.Product(FE), symbols


# Lookup table for handler to call when the "get" method (below) is
# called, depending on the first argument type.
call_lookup = {ufl.coefficient.Coefficient: coefficient,
               ufl.constant.Constant: constant,
               ufl.geometry.Jacobian: jacobian,
               ufl.geometry.CellCoordinate: cell_coordinate,
               ufl.geometry.FacetCoordinate: facet_coordinate,
               ufl.geometry.CellVertices: cell_vertices,
               ufl.geometry.FacetEdgeVectors: facet_edge_vectors,
               ufl.geometry.CellEdgeVectors: cell_edge_vectors,
               ufl.geometry.CellFacetJacobian: cell_facet_jacobian,
               ufl.geometry.ReferenceCellVolume: reference_cell_volume,
               ufl.geometry.ReferenceFacetVolume: reference_facet_volume,
               ufl.geometry.ReferenceCellEdgeVectors: reference_cell_edge_vectors,
               ufl.geometry.ReferenceFacetEdgeVectors: reference_facet_edge_vectors,
               ufl.geometry.ReferenceNormal: reference_normal,
               ufl.geometry.CellOrientation: _pass,
               ufl.geometry.FacetOrientation: facet_orientation,
               ufl.geometry.SpatialCoordinate: spatial_coordinate}
