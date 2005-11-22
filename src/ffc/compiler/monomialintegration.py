"This module provides efficient integration of monomial forms."

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-03 -- 2005-11-21"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# Thanks to Robert C. Kirby for suggesting the initial algorithm that
# this implementation is based on.

# Python modules
import Numeric

# FIAT modules
from FIAT.quadrature import *
from FIAT.shapes import *

# FFC common modules
from ffc.common.debug import *

# FFC compiler modules
from algebra import *
from multiindex import *

def integrate(product):
    """Compute the reference tensor for a given monomial term of a
    multilinear form, given as a Product."""

    debug("Computing reference tensor, this may take some time...")

    # Initialize quadrature points and weights
    (points, weights, vscaling, dscaling) = __init_quadrature(product.basisfunctions)

    # Initialize quadrature table for basis functions
    table = __init_table(product.basisfunctions, points)

    # Compute table Psi for each factor
    psis = [__compute_psi(v, table, len(points), dscaling) for v in product.basisfunctions]

    # Compute product of all Psis
    A0 = __compute_product(psis, weights, vscaling * product.numeric)

    return A0

def __init_quadrature(basisfunctions):
    "Initialize quadrature for given monomial term."

    # Get shape (check first one, all should be the same)
    shape = basisfunctions[0].element.shape()

    # Compute number of points to match the degree
    q = __compute_degree(basisfunctions)
    m = (q + 1 + 1) / 2 # integer division gives 2m - 1 >= q
    debug("Total degree is %d, using %d quadrature point(s) in each dimension" % (q, m), 1)

    # Create quadrature rule and get points and weights
    quadrature = make_quadrature(shape, m)
    points = quadrature.get_points()
    weights = quadrature.get_weights()

    # Compensate for different choice of reference cells in FIAT
    # FIXME: Convince Rob to change his reference elements
    if shape == TRIANGLE:
        vscaling = 0.25  # Area 1/2 instead of 2
        dscaling = 2.0   # Scaling of derivative
    elif shape == TETRAHEDRON:
        vscaling = 0.125 # Volume 1/6 instead of 4/3
        dscaling = 2.0   # Scaling of derivative
        
    return (points, weights, vscaling, dscaling)

def __init_table(basisfunctions, points):
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
        table[element] = element.tabulate(order, points)

    return table

def __compute_psi(v, table, num_points, dscaling):
    "Compute the table Psi for the given BasisFunction v."

    # We just need to pick the values for Psi from the table, which is
    # somewhat tricky since the table created by tabulate_jet() is a
    # mix of list, dictionary and Numeric.array.

    # Get FiniteElement for v
    element = v.element
    shapedim = element.shapedim()
    spacedim = element.spacedim()

    # Get Indices and shapes for Derivatives
    dindex = [d.index for d in v.derivatives]
    dshape = [shapedim for d in v.derivatives]
    dorder = len(dindex)

    # Get Indices and shapes for BasisFunction
    vindex = [v.index]
    vshape = [spacedim]

    # Get Indices and shapes for components
    if len(v.component) > 1:
        raise RuntimeError, "Can only handle rank 0 or rank 1 tensors."
    if len(v.component) > 0:
        cindex = [v.component[0]]
        cshape = [element.tensordim(0)]
    else:
        cindex = []
        cshape = []

    # Create list of Indices that label the dimensions of the tensor Psi
    indices = cindex + dindex + vindex
    shapes = cshape + dshape + vshape + [num_points]

    # Initialize tensor Psi: component, derivatives, basis function, points
    Psi = Numeric.zeros(shapes, Numeric.Float)

    # Iterate over derivative Indices
    dlists = build_indices(dshape) or [[]]
    if len(cindex) > 0:
        etable = table[element]
        for component in range(cshape[0]):
            for dlist in dlists:
                # Translate derivative multiindex to lookup tuple
                dtuple = __multiindex_to_tuple(dlist, shapedim)
                # Get values from table
                Psi[component][dlist] = etable[component][dorder][dtuple]
    else:
        etable = table[element][dorder]
        for dlist in dlists:
            # Translate derivative multiindex to lookup tuple
            dtuple = __multiindex_to_tuple(dlist, shapedim)
            # Get values from table
            Psi[dlist] = etable[dtuple]

    # Rearrange Indices as (fixed, primary, secondary, auxiliary)
    (rearrangement, num_indices) = __compute_rearrangement(indices)
    indices = [indices[i] for i in rearrangement]
    Psi = Numeric.transpose(Psi, rearrangement + (len(indices),))

    # Remove fixed indices
    for i in range(num_indices[0]):
        Psi = Psi[indices[i].index,...]
    indices = [index for index in indices if not index.type == "fixed"]

    # Scale derivatives (FIAT uses different reference element)
    Psi = pow(dscaling, dorder) * Psi

    return (Psi, indices)

def __compute_product(psis, weights, vscaling):
    "Compute special product of list of Psis."

    # We want to compute the outer product with respect to all
    # dimensions but the last of all the Psi tensors, and the
    # elementwise product with respect the last dimension
    # which corresponds to quadrature points.
    
    # Start with the quadrature weights
    A0 = vscaling * Numeric.array(weights)

    # Iterate over all Psis
    indices = []
    num_points = len(weights)
    for (Psi, index) in psis:

        # Build Index list so we can rearrange later
        indices += index

        # Initialize new tensor
        newshape = Numeric.shape(A0)[:-1] + Numeric.shape(Psi)[:-1] + (num_points,)
        B = Numeric.zeros(newshape, Numeric.Float)

        # Iterate over quadrature points
        for q in range(num_points):

            # Compute outer product with current Psi
            B[..., q] = Numeric.multiply.outer(A0[..., q], Psi[..., q])

        # Update reference tensor
        A0 = B

    # Sum over quadrature points (tensor contraction on last dimension)
    A0 = Numeric.add.reduce(A0, -1)

    # Rearrange Indices as (primary, secondary, auxiliary)
    (rearrangement, num_indices) = __compute_rearrangement(indices)
    A0 = Numeric.transpose(A0, rearrangement)

    # Sum over auxiliary indices (tensor contraction over dimension pairs)
    num_auxiliary = num_indices[3]
    assert num_auxiliary % 2 == 0
    for i in range(num_auxiliary / 2):
        A0 = Numeric.add.reduce(Numeric.diagonal(A0), -1)
    
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
    fixed     = __find_indices(indices, "fixed")
    primary   = __find_indices(indices, "primary")
    secondary = __find_indices(indices, "secondary")
    auxiliary = __find_indices(indices, "reference tensor auxiliary")
    assert len(fixed + primary + secondary + auxiliary) == len(indices)
    return (tuple(fixed + primary + secondary + auxiliary ), \
            (len(fixed), len(primary), len(secondary), len(auxiliary)))

def __find_indices(indices, type):
    "Return sorted list of positions for given Index type."
    pos = [i for i in range(len(indices)) if indices[i].type == type]
    val = [indices[i].index for i in range(len(indices)) if indices[i].type == type]
    return [pos[i] for i in Numeric.argsort(val)]

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
