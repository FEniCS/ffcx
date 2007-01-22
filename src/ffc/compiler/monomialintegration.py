"This module provides efficient integration of monomial forms."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2007-01-18"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Thanks to Robert C. Kirby for suggesting the initial algorithm that
# this implementation is based on.
#
# Modified by Garth N. Wells 2006

# Python modules
import numpy

# FIAT modules
from FIAT.quadrature import *
from FIAT.shapes import *

# FFC common modules
from ffc.common.debug import *
from ffc.common.progress import *

# FFC formlanguage modules
from ffc.formlanguage.index import *
from ffc.formlanguage.algebra import *
from ffc.formlanguage.integral import *
from ffc.formlanguage.multiindex import *

# FFC compiler modules
from pointreordering import *

def integrate(product, facet0, facet1, alignment):
    """Compute the reference tensor for a given monomial term of a
    multilinear form, given as a Product."""

    debug("Pretabulating basis functions at quadrature points")

    # Check for integral type
    type = product.integral.type

    # Initialize quadrature points and weights
    (points, weights, vscaling, dscaling) = __init_quadrature(product.basisfunctions, type)

    # Initialize quadrature table for basis functions
    table = __init_table(product.basisfunctions, type, points, facet0, facet1, alignment)

    # Compute table Psi for each factor
    psis = [__compute_psi(v, table, len(points), dscaling) for v in product.basisfunctions]

    # Compute product of all Psis
    A0 = __compute_product(psis, vscaling * product.numeric * weights)

    return A0

def __init_quadrature(basisfunctions, type):
    "Initialize quadrature for given monomial term."

    debug("Initializing quadrature.", 1)

    # Get shape (check first one, all should be the same)
    shape = basisfunctions[0].element.cell_shape()
    facet_shape = basisfunctions[0].element.facet_shape()

    # Compute number of points to match the degree
    q = __compute_degree(basisfunctions)
    m = (q + 1 + 1) / 2 # integer division gives 2m - 1 >= q
    debug("Total degree is %d, using %d quadrature point(s) in each dimension" % (q, m), 1)

    # Create quadrature rule and get points and weights
    # FIXME: FIAT ot finiteelement should return shape of facet
    if type == Integral.CELL:
        quadrature = make_quadrature(shape, m)
    elif type == Integral.EXTERIOR_FACET:
        quadrature = make_quadrature(facet_shape, m)
    elif type == Integral.INTERIOR_FACET:
        quadrature = make_quadrature(facet_shape, m)
    points = quadrature.get_points()
    weights = quadrature.get_weights()

    # Compensate for different choice of reference cells in FIAT
    # FIXME: Convince Rob to change his reference elements
    if shape == TRIANGLE:
        if type == Integral.CELL:
            vscaling = 0.25  # Area 1/2 instead of 2
            dscaling = 2.0   # Scaling of derivative
        elif type == Integral.EXTERIOR_FACET:
            vscaling = 0.5   # Length 1 instead of 2
            dscaling = 2.0   # Scaling of derivative        
        elif type == Integral.INTERIOR_FACET:
            vscaling = 0.5   # Length 1 instead of 2
            dscaling = 2.0   # Scaling of derivative        
    elif shape == TETRAHEDRON:
        if type == Integral.CELL:
            vscaling = 0.125 # Volume 1/6 instead of 4/3
            dscaling = 2.0   # Scaling of derivative
        elif type == Integral.EXTERIOR_FACET:
            vscaling = 0.25  # Area 1/2 instead of 2
            dscaling = 2.0   # Scaling of derivative
        elif type == Integral.INTERIOR_FACET:
            vscaling = 0.25  # Area 1/2 instead of 2
            dscaling = 2.0   # Scaling of derivative

    return (points, weights, vscaling, dscaling)

def __init_table(basisfunctions, type, points, facet0, facet1, alignment):
    """Initialize table of basis functions and their derivatives at
    the given quadrature points for each element."""

    debug("Precomputing table of basis functions at quadrature points.", 1)
    
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
        if type == Integral.CELL:
            table[(element, None)] = element.tabulate(order, points)
        elif type == Integral.EXTERIOR_FACET:
            table[(element, None)] = element.tabulate(order, points, facet0)
        elif type == Integral.INTERIOR_FACET:
            reordered_points = reorder_points(points, element.facet_shape(), alignment)
            table[(element, Restriction.PLUS)]  = element.tabulate(order, points, facet0)
            table[(element, Restriction.MINUS)] = element.tabulate(order, reordered_points, facet1)

    return table

def __compute_psi(v, table, num_points, dscaling):
    "Compute the table Psi for the given BasisFunction v."

    debug("Computing table for factor v = " + str(v), 1)

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
    shapedim = element.shapedim()
    space_dimension = element.space_dimension()

    # Get restriction for v
    restriction = v.restriction

    # Get Indices and shapes for Derivatives
    dindex = [d.index for d in v.derivatives]
    dshape = [shapedim for d in v.derivatives]
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
                dtuple = __multiindex_to_tuple(dlist, shapedim)
                # Get values from table
                Psi[component][tuple(dlist)] = etable[component][dorder][dtuple]
    else:
        etable = table[(element, restriction)][dorder]
        for dlist in dlists:
            # Translate derivative multiindex to lookup tuple
            dtuple = __multiindex_to_tuple(dlist, shapedim)
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

    debug("Computing product of tables", 1)

    # The reference tensor is obtained by summing over quadrature
    # points and auxiliary Indices the outer product of all the Psis
    # with the first dimension (corresponding to quadrature points)
    # and all auxiliary dimensions removed.

    # Initialize zero reference tensor (will be rearranged later)
    (shape, indices) = __compute_shape(psis)
    A0 = numpy.zeros(shape, dtype= numpy.float)
    debug("Computing the reference tensor (%d entries), this may take some time..." % numpy.size(A0))

    # Initialize list of auxiliary multiindices
    bshape = __compute_auxiliary_shape(psis)
    bindices = build_indices(bshape) or [[]]

    # Sum over quadrature points and auxiliary indices
    num_points = len(weights)
    progress = Progress(num_points * len(bindices))
    for q in range(num_points):
        for b in bindices:            
            # Compute outer products of subtables for current (q, b)
            B = weights[q]
            for (Psi, index, bpart) in psis:
                B = numpy.multiply.outer(B, Psi[ tuple([q] + [b[i] for i in bpart])])

            # Add product to reference tensor
            numpy.add(A0, B, A0)

            # Update progress
            progress += 1

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

def __find_indices(indices, type):
    "Return sorted list of positions for given Index type."
    pos = [i for i in range(len(indices)) if indices[i].type == type]
    val = [indices[i].index for i in range(len(indices)) if indices[i].type == type]
    return [pos[i] for i in numpy.argsort(val)]

def __multiindex_to_tuple(dindex, shapedim):
    """Compute lookup tuple from given derivative
    multiindex. Necessary since the table we get from FIAT is a
    dictionary with the tuples as keys. A derivative tuple specifies
    the number of derivatives in each space dimension, rather than
    listing the space dimensions for the derivatives."""
    dtuple = [0 for i in range(shapedim)]
    for d in dindex:
        dtuple[d] += 1
    return tuple(dtuple)
