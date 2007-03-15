"This module implements efficient integration of monomial forms"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2007-03-14"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Thanks to Robert C. Kirby for suggesting the initial algorithm that
# this implementation is based on.

# Modified by Garth N. Wells 2006

# Python modules
import numpy
import time

# FIAT modules
from FIAT.quadrature import *
from FIAT.shapes import *

# FFC common modules
from ffc.common.debug import *
from ffc.common.progress import *

# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.algebra import *
from ffc.compiler.language.integral import *

# FFC tensor representation modules
from multiindex import *
from pointreordering import *

def integrate(monomial, facet0, facet1):
    """Compute the reference tensor for a given monomial term of a
    multilinear form"""

    tic = time.time()

    # Check for integral type
    integral_type = monomial.integral.type

    # Initialize quadrature points and weights
    (points, weights, vscaling, dscaling) = __init_quadrature(monomial.basisfunctions, integral_type)

    # Initialize quadrature table for basis functions
    table = __init_table(monomial.basisfunctions, integral_type, points, facet0, facet1)

    # Compute table Psi for each factor
    psis = [__compute_psi(v, table, len(points), dscaling) for v in monomial.basisfunctions]

    # Compute product of all Psis
    A0 = __compute_product(psis, vscaling * monomial.numeric * weights)

    # Report elapsed time and number of entries
    toc = time.time() - tic
    num_entries = numpy.prod(numpy.shape(A0))
    debug("%d entries computed in %.3g seconds" % (num_entries, toc), 1)

    return A0

def __init_quadrature(basisfunctions, integral_type):
    "Initialize quadrature for given monomial term."

    # Get shape (check first one, all should be the same)
    shape = basisfunctions[0].element.cell_shape()
    facet_shape = basisfunctions[0].element.facet_shape()

    # Compute number of points to match the degree
    q = __compute_degree(basisfunctions)
    m = (q + 1 + 1) / 2 # integer division gives 2m - 1 >= q
    #debug("Total degree is %d, using %d quadrature point(s) in each dimension" % (q, m), 1)

    # Create quadrature rule and get points and weights
    # FIXME: FIAT ot finiteelement should return shape of facet
    if integral_type == Integral.CELL:
        quadrature = make_quadrature(shape, m)
    elif integral_type == Integral.EXTERIOR_FACET:
        quadrature = make_quadrature(facet_shape, m)
    elif integral_type == Integral.INTERIOR_FACET:
        quadrature = make_quadrature(facet_shape, m)
    points = quadrature.get_points()
    weights = quadrature.get_weights()

    # Compensate for different choice of reference cells in FIAT
    # FIXME: Convince Rob to change his reference elements
    if shape == TRIANGLE:
        if integral_type == Integral.CELL:
            vscaling = 0.25  # Area 1/2 instead of 2
            dscaling = 2.0   # Scaling of derivative
        elif integral_type == Integral.EXTERIOR_FACET:
            vscaling = 0.5   # Length 1 instead of 2
            dscaling = 2.0   # Scaling of derivative        
        elif integral_type == Integral.INTERIOR_FACET:
            vscaling = 0.5   # Length 1 instead of 2
            dscaling = 2.0   # Scaling of derivative        
    elif shape == TETRAHEDRON:
        if integral_type == Integral.CELL:
            vscaling = 0.125 # Volume 1/6 instead of 4/3
            dscaling = 2.0   # Scaling of derivative
        elif integral_type == Integral.EXTERIOR_FACET:
            vscaling = 0.25  # Area 1/2 instead of 2
            dscaling = 2.0   # Scaling of derivative
        elif integral_type == Integral.INTERIOR_FACET:
            vscaling = 0.25  # Area 1/2 instead of 2
            dscaling = 2.0   # Scaling of derivative

    return (points, weights, vscaling, dscaling)

def __init_table(basisfunctions, integral_type, points, facet0, facet1):
    """Initialize table of basis functions and their derivatives at
    the given quadrature points for each element."""

    # Compute maximum number of derivatives for each element
    num_derivatives = {}
    for v in basisfunctions:
        element = v.element
        order = len(v.derivatives)
        if element in num_derivatives:
            num_derivatives[element] = max(order, num_derivatives[element])
        else:
            num_derivatives[element] = order

    # Call FIAT to tabulate the basis functions for each element
    table = {}
    for element in num_derivatives:
        order = num_derivatives[element]
        # Tabulate for different integral types
        if integral_type == Integral.CELL:
            table[(element, None)] = element.tabulate(order, points)
        elif integral_type == Integral.EXTERIOR_FACET:
            table[(element, None)] = element.tabulate(order, points, facet0)
        elif integral_type == Integral.INTERIOR_FACET:
            points0 = reorder_points(points, facet0, element.cell_shape())
            points1 = reorder_points(points, facet1, element.cell_shape())
            table[(element, Restriction.PLUS)]  = element.tabulate(order, points0, facet0)
            table[(element, Restriction.MINUS)] = element.tabulate(order, points1, facet1)

    return table

def __compute_psi(v, table, num_points, dscaling):
    "Compute the table Psi for the given BasisFunction v."

    # We just need to pick the values for Psi from the table, which is
    # somewhat tricky since the table created by tabulate_jet() is a
    # mix of list, dictionary and numpy.array.
    #
    # The dimensions of the resulting table are ordered as follows:
    #
    #     one dimension  corresponding to quadrature points
    #     all dimensions corresponding to auxiliary Indices
    #     all dimensions corresponding to primary   Indices
    #     all dimensions corresponding to secondary Indices
    #
    # All fixed Indices are removed here. The first set of dimensions
    # corresponding to quadrature points and auxiliary Indices are removed
    # later when we sum over these dimensions.

    # Get FiniteElement for v
    element = v.element
    cell_dimension = element.cell_dimension()
    space_dimension = element.space_dimension()

    # Get restriction for v
    restriction = v.restriction

    # Get Indices and shapes for Derivatives
    dindex = [d.index for d in v.derivatives]
    dshape = [cell_dimension for d in v.derivatives]
    dorder = len(dindex)

    # Get Indices and shapes for BasisFunction
    vindex = [v.index]
    vshape = [space_dimension]

    # Get Indices and shapes for components
    if len(v.component) > 1:
        raise RuntimeError, "Can only handle rank 0 or rank 1 tensors."
    if len(v.component) > 0:
        cindex = [v.component[0]]
        cshape = [element.value_dimension(0)]
    else:
        cindex = []
        cshape = []

    # Create list of Indices that label the dimensions of the tensor Psi
    indices = cindex + dindex + vindex
    shapes = cshape + dshape + vshape + [num_points]

    # Initialize tensor Psi: component, derivatives, basis function, points
    Psi = numpy.zeros(shapes, dtype = numpy.float)

    # Iterate over derivative Indices
    dlists = build_indices(dshape) or [[]]
    if len(cindex) > 0:
        etable = table[(element, restriction)]
        for component in range(cshape[0]):
            for dlist in dlists:
                # Translate derivative multiindex to lookup tuple
                dtuple = __multiindex_to_tuple(dlist, cell_dimension)
                # Get values from table
                Psi[component][tuple(dlist)] = etable[component][dorder][dtuple]
    else:
        etable = table[(element, restriction)][dorder]
        for dlist in dlists:
            # Translate derivative multiindex to lookup tuple
            dtuple = __multiindex_to_tuple(dlist, cell_dimension)
            # Get values from table
            Psi[tuple(dlist)] = etable[dtuple]

    # Rearrange Indices as (fixed, auxiliary, primary, secondary)
    (rearrangement, num_indices) = __compute_rearrangement(indices)
    indices = [indices[i] for i in rearrangement]
    Psi = numpy.transpose(Psi, rearrangement + (len(indices),))

    # Remove fixed indices
    for i in range(num_indices[0]):
        Psi = Psi[indices[i].index,...]
    indices = [index for index in indices if not index.type == Index.FIXED]

    # Put quadrature points first
    rank = numpy.rank(Psi)
    Psi = numpy.transpose(Psi, (rank - 1,) + tuple(range(0, rank - 1)))

    # Scale derivatives (FIAT uses different reference element)
    Psi = pow(dscaling, dorder) * Psi

    # Compute auxiliary index positions for current Psi
    bpart = [i.index for i in indices if i.type == Index.AUXILIARY_0]

    return (Psi, indices, bpart)

def __compute_product(psis, weights):
    "Compute special product of list of Psis."

    # The reference tensor is obtained by summing over quadrature
    # points and auxiliary Indices the outer product of all the Psis
    # with the first dimension (corresponding to quadrature points)
    # and all auxiliary dimensions removed.

    # Initialize zero reference tensor (will be rearranged later)
    (shape, indices) = __compute_shape(psis)
    A0 = numpy.zeros(shape, dtype= numpy.float)

    # Initialize list of auxiliary multiindices
    bshape = __compute_auxiliary_shape(psis)
    bindices = build_indices(bshape) or [[]]

    # Sum over quadrature points and auxiliary indices
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
    (rearrangement, num_indices) = __compute_rearrangement(indices)
    A0 = numpy.transpose(A0, rearrangement)

    
    
    return A0

def __compute_degree(basisfunctions):
    "Compute total degree for given monomial term."
    q = 0
    for v in basisfunctions:
        q += v.element.degree()
        for d in v.derivatives:
            q -= 1
    return q

def __compute_rearrangement(indices):
    """Compute rearrangement tuple for given list of Indices, so that
    the tuple reorders the given list of Indices with fixed, primary,
    secondary and auxiliary Indices in rising order."""
    fixed     = __find_indices(indices, Index.FIXED)
    auxiliary = __find_indices(indices, Index.AUXILIARY_0)
    primary   = __find_indices(indices, Index.PRIMARY)
    secondary = __find_indices(indices, Index.SECONDARY)
    assert len(fixed + auxiliary + primary + secondary) == len(indices)
    return (tuple(fixed + auxiliary + primary + secondary), \
            (len(fixed), len(auxiliary), len(primary), len(secondary)))

def __compute_shape(psis):
    "Compute shape of reference tensor from given list of tables."
    shape, indices = [], []
    for (Psi, index, bpart) in psis:
        num_auxiliary = len([0 for i in index if i.type == Index.AUXILIARY_0])
        shape += numpy.shape(Psi)[1 + num_auxiliary:]
        indices += index[num_auxiliary:]
    return (shape, indices)
    
def __compute_auxiliary_shape(psis):
    """Compute shape for auxiliary indices from given list of tables.
    Also compute a list of  mappings from each table to the auxiliary
    dimensions associated with that table."""
    # First find the number of different auxiliary indices (check maximum)
    bs = [b for (Psi, index, bpart) in psis for b in bpart]
    if len(bs) == 0: return []
    bmax = max(bs)
    # Find the dimension for each auxiliary index
    bshape = [0 for i in range(bmax + 1)]
    for (Psi, index, bpart) in psis:
        for i in range(len(bpart)):
            bshape[bpart[i]] = numpy.shape(Psi)[i + 1]
    # Check that we found the shape for each auxiliary index
    if 0 in bshape:
        raise RuntimeError, "Unable to compute the shape for each auxiliary index."
    return bshape

def __find_indices(indices, index_type):
    "Return sorted list of positions for given Index type."
    pos = [i for i in range(len(indices)) if indices[i].type == index_type]
    val = [indices[i].index for i in range(len(indices)) if indices[i].type == index_type]
    return [pos[i] for i in numpy.argsort(val)]

def __multiindex_to_tuple(dindex, cell_dimension):
    """Compute lookup tuple from given derivative
    multiindex. Necessary since the table we get from FIAT is a
    dictionary with the tuples as keys. A derivative tuple specifies
    the number of derivatives in each space dimension, rather than
    listing the space dimensions for the derivatives."""
    dtuple = [0 for i in range(cell_dimension)]
    for d in dindex:
        dtuple[d] += 1
    return tuple(dtuple)
