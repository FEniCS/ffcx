"This module implements the tabulation of monomial forms, large parts are copied from monomialintegration.py"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-23 -- 2007-06-19"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Python modules
import numpy
import time

# FIAT modules
from FIAT.quadrature import *
from FIAT.shapes import *

# FFC common modules
from ffc.common.debug import *

# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.integral import *
from ffc.compiler.language.restriction import *

# FFC tensor representation modules
from ffc.compiler.representation.tensor.multiindex import *
from ffc.compiler.representation.tensor.pointreordering import *

class Quadrature:
    """This class contains quadrature information

    Attributes:

        points    - quadrture points
        weights   - weights in quadrature_points
    """

    def __init__(self, points, weights):
        "Create quadrture object"

        self.points  = points
        self.weights = weights

def tabulate(monomial, facet0, facet1):
    """Tabulate the element tensor for a given monomial term of a
    multilinear form"""

    tic = time.time()

    # Check for integral type
    integral_type = monomial.integral.type

    # Initialize quadrature points and weights
    (points, weights, vscaling, dscaling) = __init_quadrature(monomial.basisfunctions, integral_type)

    num_quadrature_points = len(weights)

    # Correction of weights by scaling factor (and numeric constant if any. This might be wrong!!!)
    quadrature = Quadrature(points, weights*vscaling*monomial.numeric)

    # Initialize quadrature table for basis functions
    table = __init_table(monomial.basisfunctions, integral_type, points, facet0, facet1)

    # Compute table Psi for each factor
    psis = [__compute_psi(v, table, len(points), dscaling, integral_type) for v in monomial.basisfunctions]

    toc = time.time() - tic

    return (psis, quadrature)
# Use for non-affine mapping?
#    return (map_derivatives, map_element, psis, quadrature)
#    return (Derivatives, num_quadrature_points, psis)

def __init_quadrature(basisfunctions, integral_type):
    "Initialize quadrature for given monomial term."

    # Get shape (check first one, all should be the same)
    shape = basisfunctions[0].element.cell_shape()
    facet_shape = basisfunctions[0].element.facet_shape()

    # Compute number of points to match the degree
    q = __compute_degree(basisfunctions)
    m = (q + 1 + 1) / 2 # integer division gives 2m - 1 >= q
    debug("Total degree is %d, using %d quadrature point(s) in each dimension" % (q, m), 1)

    # Create quadrature rule and get points and weights
    # FIXME: FIAT not finiteelement should return shape of facet
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
            points0 = reorder_points(points, element.cell_shape(), facet0)
            points1 = reorder_points(points, element.cell_shape(), facet1)
            table[(element, Restriction.PLUS)]  = element.tabulate(order, points0, facet0)
            table[(element, Restriction.MINUS)] = element.tabulate(order, points1, facet1)

    return table

def __compute_psi(v, table, num_points, dscaling, integral_type):
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
    # All fixed Indices are removed here.

    # Get cell dimension
    cell_dimension = v.element.cell_dimension()

    # Get indices and shapes for components
    if len(v.component) ==  0:
        cindex = []
        cshape = []
    elif len(v.component) == 1:
        cindex = [v.component[0]]
        cshape = [len(v.component[0].range)]
    else:
        raise RuntimeError, "Can only handle rank 0 or rank 1 tensors."

    # Get indices and shapes for derivatives
    dindex = [d.index for d in v.derivatives]
    dshape = [len(d.index.range) for d in v.derivatives]
    dorder = len(dindex)

    # Get indices and shapes for basis functions
    vindex = [v.index]
    vshape = [len(v.index.range)]

    # Create list of indices that label the dimensions of the tensor Psi
    indices = cindex + dindex + vindex
    shapes = cshape + dshape + vshape + [num_points]
    dimensions = cshape + dshape + vshape

    # Initialize tensor Psi: component, derivatives, basis function, points
    Psi = numpy.zeros(shapes, dtype = numpy.float)

    # Get restriction and handle constants
    restriction = v.restriction
    if restriction == Restriction.CONSTANT:
        if integral_type == Integral.INTERIOR_FACET:
            restriction = Restriction.PLUS
        else:
            restriction = None

    # Iterate over derivative indices
    dlists = build_indices([index.range for index in dindex]) or [[]]
    if len(cindex) > 0:
        etable = table[(v.element, restriction)]
        for component in range(len(cindex[0].range)):
            for dlist in dlists:
                # Translate derivative multiindex to lookup tuple
                dtuple = __multiindex_to_tuple(dlist, cell_dimension)
                # Get values from table
                Psi[component][tuple(dlist)] = etable[cindex[0].range[component]][dorder][dtuple]
    else:
        etable = table[(v.element, restriction)][dorder]
        for dlist in dlists:
            # Translate derivative multiindex to lookup tuple
            dtuple = __multiindex_to_tuple(dlist, cell_dimension)
            # Get values from table
            Psi[tuple(dlist)] = etable[dtuple]


    # Rearrange Indices as (fixed, auxiliary, secondary, primary)
    (rearrangement, num_indices) = __compute_rearrangement(indices)
    indices = [indices[i] for i in rearrangement]

    # Rearrange Psi
    Psi = numpy.transpose(Psi, rearrangement + (len(indices),))

    # Remove fixed indices from Psi and indices
    for i in range(num_indices[0]):
        Psi = Psi[0, ...]
    indices = [index for index in indices if not index.type == Index.FIXED]

    # Scale derivatives (FIAT uses different reference element)
    Psi = pow(dscaling, dorder) * Psi

    # Compute auxiliary index positions for current Psi (I don't use this? problems?)
    bpart = [i.index for i in indices if i.type == Index.AUXILIARY_0]

    # Return Psis the list of indices, and the index of the basis function
    return (Psi, indices, vindex[0])

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
    the tuple reorders the given list of Indices with fixed, auxiliary,
    secondary and primary Indices in rising order."""

    fixed     = __find_indices(indices, Index.FIXED)
    auxiliary = __find_indices(indices, Index.AUXILIARY_0)
    primary   = __find_indices(indices, Index.PRIMARY)
    secondary = __find_indices(indices, Index.SECONDARY)
    assert len(fixed + auxiliary + primary + secondary) == len(indices)
    return (tuple(fixed + auxiliary + secondary + primary), \
            (len(fixed), len(auxiliary), len(secondary), len(primary)))

def __find_indices(indices, index_type):
    "Return unsorted list of positions for given Index type."
#    pos = [i for i in range(len(indices)) if indices[i].type == index_type]
#    val = [indices[i].index for i in range(len(indices)) if indices[i].type == index_type]
#    return [pos[i] for i in numpy.argsort(val)]
    return [i for i in range(len(indices)) if indices[i].type == index_type]

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
