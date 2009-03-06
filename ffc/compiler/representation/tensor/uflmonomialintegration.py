"This module implements efficient integration of monomial forms"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2007-11-26"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Thanks to Robert C. Kirby for suggesting the initial algorithm that
# this implementation is based on.

# Modified by Garth N. Wells 2006
# Modified by Marie E. Rognes (meg@math.uio.no) 2008

# Python modules
import numpy
import time

# UFL modules
from ufl.integral import Measure

# FIAT modules
#from FIAT.quadrature import make_quadrature
#from FIAT.transformedspace import *
#from FIAT.shapes import *

# FFC common modules
from ffc.common.log import debug
#from ffc.common.progress import *

# FFC fem modules
from ffc.fem import create_element
from ffc.fem.quadrature import make_quadrature
from ffc.fem.quadratureelement import QuadratureElement

# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.algebra import *
from ffc.compiler.language.integral import *

# FFC tensor representation modules
from multiindex import *
from pointreordering import *

def integrate(monomial, domain_type, facet0, facet1):
    "Compute reference tensor for given monomial."

    # Start timing (this may take some time)
    tic = time.time()

    # Extract factors and elements
    factors = monomial.factors
    elements = [create_element(factor.element()) for factor in factors]

    # Initialize quadrature points and weights
    (points, weights) = _init_quadrature(factors, elements, domain_type)

    # Initialize quadrature table for basis functions
    table = _init_table(factors, elements, domain_type, points, facet0, facet1)

    print table

    # Compute table Psi for each factor
    psis = [_compute_psi(v, elements[i], table, len(points), domain_type) for (i, v) in enumerate(factors)]

    # Compute product of all Psis
    A0 = _compute_product(psis, monomial.float_value * weights)

    # Report elapsed time and number of entries
    toc = time.time() - tic
    num_entries = numpy.prod(numpy.shape(A0))
    debug("%d entries computed in %.3g seconds" % (num_entries, toc), 1)
    debug("Shape of reference tensor: " + str(numpy.shape(A0)), 1)

    return A0

def _init_quadrature(factors, elements, domain_type):
    "Initialize quadrature for given monomial."

    # Get shapes (check first factor, should be the same for all)
    cell_shape = elements[0].cell_shape()
    facet_shape = elements[0].facet_shape()

    # Compute number of points to match the degree
    degree = _compute_degree(factors, elements)
    num_points = (degree + 2) / 2
    debug("Total degree is %d, using %d quadrature point(s) in each dimension." % (degree, num_points))

    # Check if any basis functions are defined on a quadrature element
    for element in elements:
        if isinstance(element, QuadratureElement):
            num_points = element.num_axis_points()
            debug("Found quadrature element, adjusting number of points to %d.", num_points)
            break

    # Create quadrature rule and get points and weights
    if domain_type == Measure.CELL:
        (points, weights) = make_quadrature(cell_shape, num_points)
    else:
        (points, weights) = make_quadrature(facet_shape, num_points)

    return (points, weights)

def _init_table(factors, elements, domain_type, points, facet0, facet1):
    """Initialize table of basis functions and their derivatives at
    the given quadrature points for each element."""

    # Compute maximum number of derivatives for each element
    num_derivatives = {}
    for (i, v) in enumerate(factors):
        element = elements[i]
        order = len(v.derivative)
        if element in num_derivatives:
            num_derivatives[element] = max(order, num_derivatives[element])
        else:
            num_derivatives[element] = order

    print num_derivatives

    # Call FIAT to tabulate the basis functions for each element
    table = {}
    for element in num_derivatives:
        order = num_derivatives[element]
        # Tabulate for different integral types
        if domain_type == Measure.CELL:
            table[(element, None)] = element.tabulate(order, points)
        elif domain_type == Measure.EXTERIOR_FACET:
            table[(element, None)] = element.tabulate(order, map_to_facet(points, facet0))
        elif domain_type == Measure.INTERIOR_FACET:
            table[(element, Restriction.PLUS)]  = element.tabulate(order, map_to_facet(points, facet0))
            table[(element, Restriction.MINUS)] = element.tabulate(order, map_to_facet(points, facet1))

    return table

def _compute_psi(v, element, table, num_points, domain_type):
    "Compute the table Psi for the given basis function v."

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

    # Get cell dimension
    cell_dimension = element.cell_dimension()

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
    dindex = [d for d in v.derivative]
    dshape = [len(d.index.range) for d in v.derivative]
    dorder = len(dindex)

    # Get indices and shapes for basis functions
    vindex = [v.index]
    vshape = [len(v.index.range)]

    # Create list of indices that label the dimensions of the tensor Psi
    indices = cindex + dindex + vindex
    shapes = cshape + dshape + vshape + [num_points]

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
                dtuple = _multiindex_to_tuple(dlist, cell_dimension)
                # Get values from table
                Psi[component][tuple(dlist)] = etable[cindex[0].range[component]][dorder][dtuple]
    else:
        etable = table[(v.element, restriction)][dorder]
        for dlist in dlists:
            # Translate derivative multiindex to lookup tuple
            dtuple = _multiindex_to_tuple(dlist, cell_dimension)
            # Get values from table
            Psi[tuple(dlist)] = etable[dtuple]

    # Rearrange Indices as (fixed, auxiliary, primary, secondary)
    (rearrangement, num_indices) = _compute_rearrangement(indices)
    indices = [indices[i] for i in rearrangement]
    Psi = numpy.transpose(Psi, rearrangement + (len(indices),))

    # Remove fixed indices
    for i in range(num_indices[0]):
        Psi = Psi[0, ...]
    indices = [index for index in indices if not index.type == Index.FIXED]

    # Put quadrature points first
    rank = numpy.rank(Psi)
    Psi = numpy.transpose(Psi, (rank - 1,) + tuple(range(0, rank - 1)))

    # Compute auxiliary index positions for current Psi
    bpart = [i.index for i in indices if i.type == Index.AUXILIARY_0]

    return (Psi, indices, bpart)

def _compute_product(psis, weights):
    "Compute special product of list of Psis."

    # The reference tensor is obtained by summing over quadrature
    # points and auxiliary Indices the outer product of all the Psis
    # with the first dimension (corresponding to quadrature points)
    # and all auxiliary dimensions removed.

    # Initialize zero reference tensor (will be rearranged later)
    (shape, indices) = _compute_shape(psis)
    A0 = numpy.zeros(shape, dtype= numpy.float)
    # Initialize list of auxiliary multiindices
    bshape = _compute_auxiliary_shape(psis)
    bindices = build_indices([range(b) for b in bshape]) or [[]]

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
    (rearrangement, num_indices) = _compute_rearrangement(indices)
    A0 = numpy.transpose(A0, rearrangement)
    
    return A0

def _compute_degree(factors, elements):
    "Compute total degree for given monomial term."
    q = 0
    for (i, v) in enumerate(factors):
        q += elements[i].degree()
        for d in v.derivative:
            q -= 1
    return q

def _compute_rearrangement(indices):
    """Compute rearrangement tuple for given list of Indices, so that
    the tuple reorders the given list of Indices with fixed, primary,
    secondary and auxiliary Indices in rising order."""
    fixed     = _find_indices(indices, Index.FIXED)
    auxiliary = _find_indices(indices, Index.AUXILIARY_0)
    primary   = _find_indices(indices, Index.PRIMARY)
    secondary = _find_indices(indices, Index.SECONDARY)
    assert len(fixed + auxiliary + primary + secondary) == len(indices)
    return (tuple(fixed + auxiliary + primary + secondary), \
            (len(fixed), len(auxiliary), len(primary), len(secondary)))

def _compute_shape(psis):
    "Compute shape of reference tensor from given list of tables."
    shape, indices = [], []
    for (Psi, index, bpart) in psis:
        num_auxiliary = len([0 for i in index if i.type == Index.AUXILIARY_0])
        shape += numpy.shape(Psi)[1 + num_auxiliary:]
        indices += index[num_auxiliary:]
    return (shape, indices)
    
def _compute_auxiliary_shape(psis):
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

def _find_indices(indices, index_type):
    "Return sorted list of positions for given Index type."
    pos = [i for i in range(len(indices)) if indices[i].type == index_type]
    val = [indices[i].index for i in range(len(indices)) if indices[i].type == index_type]
    return [pos[i] for i in numpy.argsort(val)]

def _multiindex_to_tuple(dindex, cell_dimension):
    """Compute lookup tuple from given derivative
    multiindex. Necessary since the table we get from FIAT is a
    dictionary with the tuples as keys. A derivative tuple specifies
    the number of derivatives in each space dimension, rather than
    listing the space dimensions for the derivatives."""
    dtuple = [0 for i in range(cell_dimension)]
    for d in dindex:
        dtuple[d] += 1
    return tuple(dtuple)
