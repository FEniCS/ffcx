"This module implements the tabulation of monomial forms, large parts are copied from monomialintegration.py"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-23 -- 2007-03-23"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU GPL Version 2"

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
    Derivatives = __derivatives(monomial.basisfunctions, integral_type, points, dscaling, facet0, facet1)

#    print "Derivatives: ", Derivatives
#    quadrature = Quadrature(points, weights)

    # Correction of weights by scaling factor
    quadrature = Quadrature(points, weights*vscaling)
#    points[0] = (-1,-1)
#    print "\n monomial integration, points: \n", points
#    print "\n monomial integration, weights: \n", weights
#    print "\n monomial integration, vscaling: \n", vscaling
#    print "\n monomial integration, dscaling: \n", dscaling
#    print "init table"
    # Initialize quadrature table for basis functions
    table = __init_table(monomial.basisfunctions, integral_type, points, facet0, facet1)

#    print "\n monomial integration, table: \n", table

#    print "\n monomial integration, monomial.basisfunctions: \n", monomial.basisfunctions

#    print "\n monomial integration, monomial.basisfunctions[0].element: \n", \ 
#monomial.basisfunctions[0].element.__doc__

#    print "\n monomial integration, monomial.basisfunctions[0].element.basis(): \n", \
#    monomial.basisfunctions[0].element.basis()

#    print "\n monomial integration, monomial.basisfunctions[0].element.basis()[0]: \n", \
#    monomial.basisfunctions[0].element.basis()[0]

#    print "\n monomial integration, monomial.basisfunctions[0].element.basis()[0](points[0]): \n", \
#    monomial.basisfunctions[0].element.basis()[1]((-1,-1))

#    print "\n monomial integration, monomial.basisfunctions[0].element.basis()[2].deriv(0)(points[0]): \n", \
#    monomial.basisfunctions[0].element.basis()[2].deriv(0)(points[0])

#    print "\n monomial integration, monomial.basisfunctions[0].element.basis()[2].deriv(1)(points[0]): \n", \
#    monomial.basisfunctions[0].element.basis()[2].deriv(1)(points[0])

#    print "\n Computing Psis \n"
    # Compute table Psi for each factor
    psis = [__compute_psi(v, table, len(points), dscaling) for v in monomial.basisfunctions]


#    print "\n monomial integration, psis: \n", psis
#    print "numpy.shape(psis[0][0]): ", numpy.shape(psis[0][0])
#    print "numpy.shape(psis[1][0]): ", numpy.shape(psis[1][0])
#    print "numpy.shape(psis[2][0]): ", numpy.shape(psis[2][0])
#    print "psis[0][1][0].type: ", psis[0][1][0].type

#    print "\n monomial integration, psis[0][1][0]: \n", psis[0][1][0].__doc__

#    print "\n Compute product fo all Psis"
#    print "\n monomial integration, monomial.numeric: \n", monomial.numeric
#    print "\n monomial integration, __compute_product(arg2): \n", vscaling * monomial.numeric * weights

#    print "\n monomial integration, monomial: \n", monomial

    # Compute product of all Psis
#    A0 = __compute_product(psis, vscaling * monomial.numeric * weights)
#    A0 = 1
#    print "A0: ", A0
    # Report elapsed time and number of entries
    toc = time.time() - tic
#    num_entries = numpy.prod(numpy.shape(A0))
#    debug("%d entries computed in %.3g seconds" % (num_entries, toc), 1)

#    A0 = 1
    return (Derivatives, psis, quadrature)
#    return (Derivatives, num_quadrature_points, psis)

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
#        print "element: ", element
#        print "order: ", order
        if element in num_derivatives:
            num_derivatives[element] = max(order, num_derivatives[element])
        else:
            num_derivatives[element] = order

#    print "num_derivatives: ", num_derivatives
    # Call FIAT to tabulate the basis functions for each element
    table = {}
    for element in num_derivatives:
        order = num_derivatives[element]
#        print "order: ", order
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
#    print "v: ", v
#    print "table: ", table
#    print "num_points: ", num_points
#    print "dscaling: ", dscaling

    element = v.element
#    print "element: ", element

    cell_dimension = element.cell_dimension()
#    print "cell_dimension: ", cell_dimension

    space_dimension = element.space_dimension()
#    print "space_dimension: ", space_dimension

    # Get restriction for v
    restriction = v.restriction
#    print "restriction: ", restriction

#    print "v.derivatives: ", v.derivatives

    # Get Indices and shapes for Derivatives
    dindex = [d.index for d in v.derivatives]
#    print "dindex: ", dindex

    dshape = [cell_dimension for d in v.derivatives]
#    print "dshape: ", dshape

    dorder = len(dindex)
#    print "dorder: ", dorder

    # Get Indices and shapes for BasisFunction
    vindex = [v.index]
#    print "vindex: ", vindex

    vshape = [space_dimension]
#    print "vshape: ", vshape

#    print "v.component", v.component

    # Get Indices and shapes for components
    if len(v.component) > 1:
        raise RuntimeError, "Can only handle rank 0 or rank 1 tensors."
    if len(v.component) > 0:
        cindex = [v.component[0]]
        cshape = [element.value_dimension(0)]
    else:
        cindex = []
        cshape = []

#    print "cindex: ", cindex
#    print "cshape: ", cshape

    # Create list of Indices that label the dimensions of the tensor Psi
    indices = cindex + vindex + dindex
#    print "indices: ", indices
    shapes = cshape + dshape + vshape + [num_points]
    dimensions = cshape + vshape + dshape
#    print "shapes: ", shapes
    # Initialize tensor Psi: component, derivatives, basis function, points
    Psi = numpy.zeros(shapes, dtype = numpy.float)
#    print "numpy.shape(Psi): ", numpy.shape(Psi)

#    print "Psi-init:", Psi

    # Iterate over derivative Indices
    dlists = build_indices(dshape) or [[]]
#    print "dlists: ", dlists

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
#        print "etable: ", etable
        for dlist in dlists:
#            print "dlist: ", dlist
            # Translate derivative multiindex to lookup tuple
            dtuple = __multiindex_to_tuple(dlist, cell_dimension)
#            print "dtuple: ", dtuple
            # Get values from table
#            print "tuple(dlist): ", tuple(dlist)
            Psi[tuple(dlist)] = etable[dtuple]
#            print "Psi dic: ", Psi
    # Rearrange Indices as (fixed, auxiliary, primary, secondary)
    (rearrangement, num_indices) = __compute_rearrangement(indices)

#    print "rearrangement, num_indices: ", rearrangement, num_indices

    indices = [indices[i] for i in rearrangement]
#    print "indices  - rearrangement: ", indices
    #Psi = numpy.transpose(Psi, rearrangement + (len(indices),))
#    print "Psi - rearrange: ", Psi
#    print "numpy.shape(Psi): ", numpy.shape(Psi)

    # Remove fixed indices
    for i in range(num_indices[0]):
        Psi = Psi[indices[i].index,...]
    indices = [index for index in indices if not index.type == Index.FIXED]

    # Put quadrature points first
#    rank = numpy.rank(Psi)
#    print "rank: ", rank
#    print "transpose: ", (rank - 1,) + tuple(range(0, rank - 1))
#    print "numpy.shape(Psi): ", numpy.shape(Psi)
#    Psi = numpy.transpose(Psi, (rank - 1,) + tuple(range(0, rank - 1)))

    # Scale derivatives (FIAT uses different reference element)
    Psi = pow(dscaling, dorder) * Psi
#    print "scaling: ", pow(dscaling, dorder)
    # Compute auxiliary index positions for current Psi
    bpart = [i.index for i in indices if i.type == Index.AUXILIARY_0]

#    print "Tabulate.... Psi: ", Psi
#    print "numpy.shape(Psi): ", numpy.shape(Psi)

#    print "Tabulate.... indices: ", indices
#    print "Tabulate.... bpart: ", bpart

#    print "Tabulate.... Psi[0]: ", Psi[0]
#    print "Tabulate.... indices[0]: ", indices[0]
#    print "Tabulate.... indices[1]: ", indices[1]
#    print "Tabulate.... indices[0].index: ", indices[0].index
#    print "Tabulate.... indices[1].index: ", indices[1].index
#    print "Tabulate.... indices[0].type: ", indices[0].type
#    print "Tabulate.... indices[1].type: ", indices[1].type




    return (Psi, indices, dimensions, bpart)

def __derivatives(basisfunctions, integral_type, points, dscaling, facet0, facet1):
    "Compute derivatives at quadrature points for a given monomial term"


    # Get element
    #FIXME: Assuming all elements are the same, picking first
    element = basisfunctions[0].element
    print "basisfunctions: ", basisfunctions
    print "element: ", element
    print "element.num_sub_elements(): ", element.num_sub_elements()

#        if element in num_derivatives:
#            num_derivatives[element] = max(order, num_derivatives[element])
#        else:
#            num_derivatives[element] = order

    # Call FIAT to tabulate the derivatives
#    table = {}
#    for element in elements:
#        order = num_derivatives[element]
#        print "order: ", order
        # Tabulate for different integral types
    if integral_type == Integral.CELL:
#            table[(element, None)] = element.tabulate(order, points)
        derivatives = element.tabulate(1, points)
    elif integral_type == Integral.EXTERIOR_FACET:
        RuntimeError("Not implemented yet!")
#            table[(element, None)] = element.tabulate(order, points, facet0)
    elif integral_type == Integral.INTERIOR_FACET:
        RuntimeError("Not implemented yet!")

#            points0 = reorder_points(points, facet0, element.cell_shape())
#            points1 = reorder_points(points, facet1, element.cell_shape())
#            table[(element, Restriction.PLUS)]  = element.tabulate(order, points0, facet0)
#            table[(element, Restriction.MINUS)] = element.tabulate(order, points1, facet1)
    print "numpy.shape(derivatives): ", numpy.shape(derivatives)
    print "derivatives: ", derivatives
#    print "len(derivatives[1]): ", len(derivatives[1])
#    print "derivatives[1]: ", derivatives[1]

    # FIXME: this is just a quick fix
    directions = [(1,0),(0,1)]
    for d in directions:
        derivatives[1][d] = derivatives[1][d]*dscaling

#    print "derivatives[1]: ", derivatives[1]
#    print "derivatives[1][directions[0]]*dscaling: ", derivatives[1][directions[0]]*dscaling

    # Only return derivatives of the basisfunctions, multiply by scaling factor for derivatives
    return derivatives[1]

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
    return (tuple(fixed + auxiliary + primary + secondary), \
            (len(fixed), len(auxiliary), len(primary), len(secondary)))

def __compute_shape(psis):
    "Compute shape of reference tensor from given list of tables."
    shape, indices = [], []
    for (Psi, index, bpart) in psis:
        num_auxiliary = len([0 for i in index if i.type == Index.AUXILIARY_0])
        print "num_auxiliary: ", num_auxiliary
        print "numpy.shape(Psi)", numpy.shape(Psi)
        print "numpy.shape(Psi)[1:]", numpy.shape(Psi)[1:]
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
