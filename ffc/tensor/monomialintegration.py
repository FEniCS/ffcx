"This module implements efficient integration of monomial forms."

# Copyright (C) 2004-2011 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Thanks to Robert C. Kirby for suggesting the initial algorithm that
# this implementation is based on.
#
# Modified by Garth N. Wells, 2006
# Modified by Marie E. Rognes, 2008
# Modified by Kristian B. Oelgaard, 2009
#
# First added:  2004-11-03
# Last changed: 2011-11-28

# Python modules
import numpy
import time

# UFL modules
from ufl.classes import Measure
from ufl.geometry import cellname2facetname

# FFC modules
from ffc.log import info, debug, error
from ffc.fiatinterface import create_element
from ffc.fiatinterface import map_facet_points
from ffc.quadrature_schemes import create_quadrature

# FFC tensor representation modules
from multiindex import build_indices
from monomialextraction import MonomialException
from monomialtransformation import MonomialIndex

def integrate(monomial,
              domain_type,
              facet0, facet1,
              quadrature_degree,
              quadrature_rule,
              cellname,
              facet_cellname):
    """Compute the reference tensor for a given monomial term of a
    multilinear form"""

    info("Precomputing integrals on reference element")

    # Start timing
    tic = time.time()

    # Initialize quadrature points and weights
    (points, weights) = _init_quadrature(monomial.arguments,
                                         domain_type,
                                         quadrature_degree,
                                         quadrature_rule,
                                         cellname,
                                         facet_cellname)

    # Initialize quadrature table for basis functions
    table = _init_table(monomial.arguments,
                        domain_type,
                        points,
                        facet0, facet1)

    # Compute table Psi for each factor
    psis = [_compute_psi(v, table, len(points), domain_type) \
                for v in monomial.arguments]

    # Compute product of all Psis
    A0 = _compute_product(psis, monomial.float_value * weights)

    # Report elapsed time and number of entries
    toc = time.time() - tic
    num_entries = numpy.prod(numpy.shape(A0))
    debug("%d entries computed in %.3g seconds" % (num_entries, toc))
    debug("Shape of reference tensor: " + str(numpy.shape(A0)))

    return A0

def _init_quadrature(arguments, domain_type, quadrature_degree, quadrature_rule, cellname, facet_cellname):
    "Initialize quadrature for given monomial."
    # Create quadrature rule and get points and weights
    if domain_type == Measure.CELL:
        (points, weights) = create_quadrature(cellname, quadrature_degree, quadrature_rule)
    else:
        (points, weights) = create_quadrature(facet_cellname, quadrature_degree, quadrature_rule)

    return (points, weights)

def _init_table(arguments, domain_type, points, facet0, facet1):
    """Initialize table of basis functions and their derivatives at
    the given quadrature points for each element."""

    # Compute maximum number of derivatives for each element
    num_derivatives = {}
    for v in arguments:
        ufl_element = v.element
        order = len(v.derivatives)
        if ufl_element in num_derivatives:
            num_derivatives[ufl_element] = max(order, num_derivatives[ufl_element])
        else:
            num_derivatives[ufl_element] = order

    # Call FIAT to tabulate the basis functions for each element
    table = {}
    for (ufl_element, order) in num_derivatives.items():
        fiat_element = create_element(ufl_element)
        if domain_type == Measure.CELL:
            table[(ufl_element, None)] = fiat_element.tabulate(order, points)
        elif domain_type == Measure.EXTERIOR_FACET:
            x = map_facet_points(points, facet0)
            table[(ufl_element, None)] = fiat_element.tabulate(order, x)
        elif domain_type == Measure.INTERIOR_FACET:
            x0 = map_facet_points(points, facet0)
            x1 = map_facet_points(points, facet1)
            table[(ufl_element, "+")] = fiat_element.tabulate(order, x0)
            table[(ufl_element, "-")] = fiat_element.tabulate(order, x1)

    return table

def _compute_psi(v, table, num_points, domain_type):
    "Compute the table Psi for the given basis function v."

    # We just need to pick the values for Psi from the table, which is
    # somewhat tricky since the table created by tabulate_jet() is a
    # mix of list, dictionary and numpy.array.
    #
    # The dimensions of the resulting table are ordered as follows:
    #
    #     one dimension  corresponding to quadrature points
    #     all dimensions corresponding to internal  Indices
    #     all dimensions corresponding to primary   Indices
    #     all dimensions corresponding to secondary Indices
    #
    # All fixed Indices are removed here. The first set of dimensions
    # corresponding to quadrature points and internal Indices are removed
    # later when we sum over these dimensions.

    # Get topological dimension of cell
    tdim = v.element.cell().topological_dimension()

    # Get indices and shapes for components
    if len(v.components) ==  0:
        cindex = []
        cshape = []
    elif len(v.components) == 1:
        cindex = [v.components[0]]
        cshape = [len(v.components[0].index_range)]
    else:
        raise MonomialException, "Can only handle rank 0 or rank 1 tensors."

    # Get indices and shapes for derivatives
    dindex = [d for d in v.derivatives]
    dshape = [len(d.index_range) for d in v.derivatives]

    # Get indices and shapes for basis functions
    vindex = [v.index]
    vshape = [len(v.index.index_range)]

    # Create list of indices that label the dimensions of the tensor Psi
    indices = cindex + dindex + vindex
    shapes = cshape + dshape + vshape + [num_points]

    # Initialize tensor Psi: component, derivatives, basis function, points
    Psi = numpy.zeros(shapes, dtype = numpy.float)

    # Iterate over derivative indices
    dlists = build_indices([index.index_range for index in dindex]) or [[]]
    if len(cindex) > 0:
        etable = table[(v.element, v.restriction)]
        for component in range(len(cindex[0].index_range)):
            for dlist in dlists:
                # Translate derivative multiindex to lookup tuple
                dtuple = _multiindex_to_tuple(dlist, tdim)
                # Get values from table
                Psi[component][tuple(dlist)] = \
                    etable[dtuple][:, cindex[0].index_range[component], :]
    else:
        etable = table[(v.element, v.restriction)]
        for dlist in dlists:
            # Translate derivative multiindex to lookup tuple
            dtuple = _multiindex_to_tuple(dlist, tdim)
            # Get values from table
            Psi[tuple(dlist)] = etable[dtuple]

    # Rearrange Indices as (fixed, internal, primary, secondary)
    (rearrangement, num_indices) = _compute_rearrangement(indices)
    indices = [indices[i] for i in rearrangement]
    Psi = numpy.transpose(Psi, rearrangement + (len(indices),))

    # Remove fixed indices
    for i in range(num_indices[0]):
        Psi = Psi[0, ...]
    indices = [index for index in indices \
                   if not index.index_type == MonomialIndex.FIXED]

    # Put quadrature points first
    rank = numpy.rank(Psi)
    Psi = numpy.transpose(Psi, (rank - 1,) + tuple(range(0, rank - 1)))

    # Compute internal index positions for current Psi
    bpart = [i.index_id for i in indices if i.index_type == MonomialIndex.INTERNAL]

    return (Psi, indices, bpart)

def _compute_product(psis, weights):
    "Compute special product of list of Psis."

    # The reference tensor is obtained by summing over quadrature
    # points and internal Indices the outer product of all the Psis
    # with the first dimension (corresponding to quadrature points)
    # and all internal dimensions removed.

    # Initialize zero reference tensor (will be rearranged later)
    (shape, indices) = _compute_shape(psis)
    A0 = numpy.zeros(shape, dtype= numpy.float)

    # Initialize list of internal multiindices
    bshape = _compute_internal_shape(psis)
    bindices = build_indices([range(b) for b in bshape]) or [[]]

    # Sum over quadrature points and internal indices
    num_points = len(weights)
    for q in range(num_points):
        for b in bindices:

            # Compute outer products of subtables for current (q, b)
            B = weights[q]
            for (Psi, index, bpart) in psis:
                B = numpy.multiply.outer(B, Psi[ tuple([q] + [b[i] for i in bpart])])

            # Add product to reference tensor
            numpy.add(A0, B, A0)

    # Rearrange Indices as (primary, secondary)
    (rearrangement, num_indices) = _compute_rearrangement(indices)
    A0 = numpy.transpose(A0, rearrangement)

    return A0

def _compute_rearrangement(indices):
    """
    Compute rearrangement tuple for given list of Indices, so that the
    tuple reorders the given list of Indices with fixed, primary,
    secondary and internal Indices in rising order.
    """
    fixed     = _find_indices(indices, MonomialIndex.FIXED)
    internal  = _find_indices(indices, MonomialIndex.INTERNAL)
    primary   = _find_indices(indices, MonomialIndex.PRIMARY)
    secondary = _find_indices(indices, MonomialIndex.SECONDARY)
    assert len(fixed + internal + primary + secondary) == len(indices)
    return (tuple(fixed + internal + primary + secondary), \
            (len(fixed), len(internal), len(primary), len(secondary)))

def _compute_shape(psis):
    "Compute shape of reference tensor from given list of tables."
    shape, indices = [], []
    for (Psi, index, bpart) in psis:
        num_internal = len([0 for i in index if i.index_type == MonomialIndex.INTERNAL])
        shape += numpy.shape(Psi)[1 + num_internal:]
        indices += index[num_internal:]
    return (shape, indices)

def _compute_internal_shape(psis):
    """
    Compute shape for internal indices from given list of tables.
    Also compute a list of mappings from each table to the internal
    dimensions associated with that table.
    """
    # First find the number of different internal indices (check maximum)
    bs = [b for (Psi, index, bpart) in psis for b in bpart]
    if len(bs) == 0: return []
    bmax = max(bs)
    # Find the dimension for each internal index
    bshape = [0 for i in range(bmax + 1)]
    for (Psi, index, bpart) in psis:
        for i in range(len(bpart)):
            bshape[bpart[i]] = numpy.shape(Psi)[i + 1]
    # Check that we found the shape for each internal index
    if 0 in bshape:
        error("Unable to compute the shape for each internal index.")
    return bshape

def _find_indices(indices, index_type):
    "Return sorted list of positions for given index type."
    pos = [i for i in range(len(indices)) if indices[i].index_type == index_type]
    val = [indices[i].index_id for i in range(len(indices)) if indices[i].index_type == index_type]
    return [pos[i] for i in numpy.argsort(val)]

def _multiindex_to_tuple(dindex, cell_dimension):
    """
    Compute lookup tuple from given derivative multiindex. Necessary
    since the table we get from FIAT is a dictionary with the tuples
    as keys. A derivative tuple specifies the number of derivatives in
    each space dimension, rather than listing the space dimensions for
    the derivatives.
    """
    dtuple = [0 for i in range(cell_dimension)]
    for d in dindex:
        dtuple[d] += 1
    return tuple(dtuple)
