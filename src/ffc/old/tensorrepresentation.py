"""This module implements the computation of the tensor representation
for a given subset of terms in a sum, specified as the terms matching
a given integral type."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-06 -- 2007-01-31"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006

# FFC common modules
from debug import *
from constants import *
from exceptions import *

# FFC compiler.language modules
from indexreordering import *

# FFC compiler modules
from term import *

# FFC optimization modules
from optimization import *

# FIXME: Remove
from declaration import *

def compute_terms(form, type, facet0, facet1, alignment):

    # Check if there are any terms to compute
    num_terms = __terms_to_compile(form, type)
    debug("Number of terms to compile: %d" % num_terms)
    if num_terms == 0:
        return []

    # Reorder indices and compute factorization
    factorization = reorder_indices(form)
#    print "factorization ",factorization

    # Compute terms
    terms = [None for i in range(len(form.monomials))]
    for i in range(len(form.monomials)):
        debug("Compiling term %d" % i, 1)
        p = form.monomials[i]
#        print "p ", p
#        print "i ",i
        if p.integral.type == type:
            # Compute geometry tensor
#            print "ok type"
            G = GeometryTensor(p)
            # Check if reference tensor should be computed
            if factorization[i] == None:
#                print "factor "
                # Compute reference tensor and add term
                A0 = ReferenceTensor(p, facet0, facet1, alignment)
                terms[i] = Term(p, A0, [G])
            else:
                # Add geometry tensor to previous term
                terms[factorization[i]].G += [G]
                A0 = terms[factorization[i]].A0
            print "Shapes: " + str(A0.i.dims) + " " + str(A0.a.dims)

    debug("All terms compiled", 1)

    # Remove terms not computed (factorized)
    [terms.remove(None) for i in range(len(terms)) if None in terms]

    return terms

def compute_reference_tensor(terms, format):
    "Precompute reference tensor according to given format."
    debug("Generating code for reference tensor", 1)
    if not terms or format.format["reference tensor"](0, 0, []) == None: return []
    declarations = []
    for j in range(len(terms)):
        term = terms[j]
        iindices = term.A0.i.indices or [[]]
        aindices = term.A0.a.indices or [[]]
        for i in iindices:
            for a in aindices:
                name = format.format["reference tensor"](j, i, a)
                value = format.format["floating point"](term.A0(i, a))
                declarations += [Declaration(name, value)]

    return declarations

def compute_geometry_tensor(terms, format, g_used, c_used):
    "Precompute geometry tensor according to given format."
    debug("Generating code for geometry tensor", 1)
    if not terms or format.format["geometry tensor"](0, []) == None: return []
    declarations = []
    for j in range(len(terms)):
        # Should be the same, so pick first
        aindices = terms[j].G[0].a.indices
        if not aindices:
            aindices = [[]]
        for a in aindices:
            # Sum factorized values
            name = format.format["geometry tensor"](j, a)
            used = name in g_used
            value = format.format["sum"]([G(a, format, c_used, used) for G in terms[j].G])
            declaration = Declaration(name, value)
            declaration.used = used
            # Add declaration
            declarations += [declaration]

    return declarations

def compute_element_tensor(terms, format, options):
    "Precompute element tensor, including possible optimizations."
    if not terms or format.format["element tensor"]((0,), 0) == None: return []
    if options["optimize"]:
        rank = terms[0].A0.i.rank
        if rank == 2 or rank == 1:
            return __compute_element_tensor_optimized(terms, format)
        else:
            debug("Only rank 2 tensors can currently be optimized with FErari, generating default code")
            return __compute_element_tensor_default(terms, format)
    else:
        return __compute_element_tensor_default(terms, format)

def check_used(terms, format, g_used):
    """Check which declarations of g are actually used, i.e,
    which entries of the geometry tensor that get multiplied with
    nonzero entries of the reference tensor."""
    if not terms or format.format["geometry tensor"](0, []) == None: return []
    iindices = terms[0].A0.i.indices or [[]] # All primary ranks are equal
    for i in iindices:
        for j in range(len(terms)):
            A0 = terms[j].A0
            if A0.a.indices: aindices = A0.a.indices
            else: aindices = [[]]
            for a in aindices:
                a0 = A0(tuple(i), tuple(a))
                gk = format.format["geometry tensor"](j, a)
                if abs(a0) > FFC_EPSILON:
                    g_used.add(gk)

def __compute_element_tensor_default(terms, format):
    """Precompute element tensor without optimizations except for
    dropping multiplication with zero."""
    debug("Generating code for element tensor", 1)         
    declarations = []
    iindices = terms[0].A0.i.indices or [[]] # All primary ranks are equal

    # Prefetch formats to speed up code generation
    format_element_tensor  = format.format["element tensor"]
    format_geometry_tensor = format.format["geometry tensor"]
    format_sum             = format.format["sum"]
    format_subtract        = format.format["subtract"]
    format_multiply        = format.format["multiplication"]
    format_floating_point  = format.format["floating point"]

    # Generate code for geometry tensor elements
    gk_tensor = [ ( [(format_geometry_tensor(j, a), a) \
                     for a in __aindices(terms[j].A0) ], j) \
                  for j in range(len(terms)) ]

    # Generate code for computing the element tensor
    k = 0
    num_dropped = 0
    num_ops = 0
    zero = format_floating_point(0.0)
    for i in iindices:
        value = None
        name = format_element_tensor(i, k)
        for (gka, j) in gk_tensor:
            A0 = terms[j].A0
            for (gk, a) in gka:
                a0 = A0.A0[tuple(i + a)]
                if abs(a0) > FFC_EPSILON:
                    if value and a0 < 0.0:
                        value = format_subtract([value, format_multiply([format_floating_point(-a0), gk])])
                    elif value:
                        value = format_sum([value, format_multiply([format_floating_point(a0), gk])])
                    else:
                        value = format_multiply([format_floating_point(a0), gk])
                    num_ops += 1
                else:
                    num_dropped += 1
        value = value or zero
        declarations += [Declaration(name, value)]
        k += 1
    debug("Number of zeros dropped from reference tensor: " + str(num_dropped), 1)
    return (declarations, num_ops)

def __compute_element_tensor_optimized(terms, format):
    "Precompute element tensor with FErari optimizations."
    debug("Generating optimized code for element tensor", 1)
    # Call FErari to do optimizations
    return optimize(terms, format)

def __terms_to_compile(form, type):
    "Count the number of terms to be computed."
    count = 0
    for p in form.monomials:
        if p.integral.type == type:
            count += 1
    return count

def __aindices(A0):
    if A0.a.indices: aindices = A0.a.indices
    else: aindices = [[]]
    return aindices
