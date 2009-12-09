"Reordering of entries in reference tensor for interior facets"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2006-12-01"
__copyright__ = "Copyright (C) 2006-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-09

# Python modules.
import numpy

# FFC tensor representation modules.
from monomialtransformation import MonomialIndex
from multiindex import MultiIndex

def reorder_entries(terms):
    """Reorder entries to compute the reference tensor for an
    interior facet from the the reduced reference tensor"""

    # Iterate over terms
    for term in terms:

        # Compute restrictions corresponding to indices
        (restrictions, idims, adims) = __compute_restrictions(term)
        dims = idims + adims

        # Compute position where to insert
        position = []
        for i in range(len(restrictions)):
            dim = dims[i]
            if restrictions[i] == "+":
                position = position + [slice(0, dim/2)]
            elif restrictions[i] == "-":
                position = position + [slice(dim/2, dim)]
            else:
                position = position + [slice(0, dim)]

        # Initialize empty reference tensor of double size in each dimension
        A0 = numpy.zeros(dims, dtype=numpy.float)

        # Insert reduced reference tensor into reference tensor
        A0[position] = term.A0.A0
        term.A0.A0 = A0

        # Reinitialize indices to new size
        term.A0.primary_multi_index = MultiIndex([range(idim) for idim in idims])
        term.A0.secondary_multi_index = MultiIndex([range(adim) for adim in adims])
        for G in term.GK:
            G.secondary_multi_index = term.A0.secondary_multi_index

def __compute_restrictions(term):
    """Compute restrictions corresponding to indices for given
    term. For indices at basis functions, we need to double the
    size of the reference tensor in the corresponding dimension,
    but for other dimensions corresponding to component indices
    and derivatives, the size remains the same."""

    # Get dimensions for primary and secondary indices
    idims = term.A0.primary_multi_index.dims
    adims = term.A0.secondary_multi_index.dims

    # Get basis functions for term
    arguments = term.monomial.arguments

    # Create empty list of restrictions for indices
    restrictions = [None for i in range(len(idims) + len(adims))]

    # Extract restrictions corresponding to primary indices at basis functions
    for i in range(len(idims)):
        for v in arguments:
            if v.index.index_type == MonomialIndex.PRIMARY and v.index.index_id == i:
                restrictions[i] = v.restriction
                break

    # Extract restrictions corresponding to secondary indices at basis functions
    for i in range(len(adims)):
        for v in arguments:
            if v.index.index_type == MonomialIndex.SECONDARY and v.index.index_id == i:
                restrictions[len(idims) + i] = v.restriction
                break

    # Compute new dimensions
    new_idims = [i for i in idims]
    new_adims = [i for i in adims]
    for i in range(len(new_idims)):
        if not restrictions[i] == None:
            new_idims[i] = 2*new_idims[i]
    for i in range(len(new_adims)):
        if not restrictions[i + len(new_idims)] == None:
            new_adims[i] = 2*new_adims[i]

    return (restrictions, new_idims, new_adims)
