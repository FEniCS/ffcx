__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-06 -- 2006-10-17"
__copyright__ = "Copyright (C) 2004-2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006

# Python modules
from sets import Set # (Could use built-in set with Python 2.4)

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *
from ffc.common.exceptions import *

# FFC compiler modules
from term import *
from reorder import *
from declaration import *
from optimization import *

class ElementTensor:
    """An ElementTensor represents the element tensor of a
    multi-linear form and consist of a list of Terms, each containing
    a pair of a ReferenceTensor and a GeometryTensor.

    Attributes:

        terms   - a list of Terms (products A0 * GK)
        a0      - a list of precomputed reference tensor declarations
        gK      - a list of precomputed geometry tensor declarations
        aK      - a list of precomputed element tensor declarations
        facet   - number of the associated facet (None for interior)
        num_ops - number of operations in computation of element tensor
    """

    def __init__(self, sum, type, format, cK_used, gK_used, options, facet):
        "Create ElementTensor."

        # Check that all terms have integrals
        self.__check_integrals(sum)

        # Reset number of operations
        self.num_ops = 0

        # Check if there are any terms to compute
        num_terms = self.__terms_to_compile(sum, type)
        debug("Number of terms to compile: %d" % num_terms)
        if num_terms == 0:
            self.terms = []
            return

        # Reorder indices and compute factorization
        factorization = reorder_indices(sum)

        # Compute terms
        self.terms = [None for i in range(len(sum.products))]
        for i in range(len(sum.products)):
            debug("Compiling term %d" % i, 1)
            p = sum.products[i]
            if p.integral.type == type:
                # Compute geometry tensor
                GK = GeometryTensor(p)
                # Check if reference tensor should be computed
                if factorization[i] == None:
                    # Compute reference tensor and add term
                    A0 = ReferenceTensor(p, facet)
                    self.terms[i] = Term(p, A0, [GK])
                else:
                    # Add geometry tensor to previous term
                    self.terms[factorization[i]].GKs += [GK]
        debug("All terms compiled", 1)

        # Remove terms not computed (factorized)
        [self.terms.remove(None) for i in range(len(self.terms)) if None in self.terms]
    
        # Compute reference tensor declarations (just need to pick the values)
        self.a0 = self.__compute_reference_tensor(format)

        # Compute element tensor declarations
        self.aK = self.__compute_element_tensor(format, options)

        # Compute geometry tensor declarations
        self.__check_used(format, gK_used)
        self.gK = self.__compute_geometry_tensor(format, gK_used, cK_used)

        # Save facet
        self.facet = facet

        return

    def __compute_reference_tensor(self, format):
        "Precompute reference tensor according to given format."
        debug("Generating code for reference tensor", 1)
        if not self.terms or format.format["reference tensor"](0, 0, []) == None: return []
        declarations = []
        for j in range(len(self.terms)):
            term = self.terms[j]
            iindices = term.A0.i.indices or [[]]
            aindices = term.A0.a.indices or [[]]
            for i in iindices:
                for a in aindices:
                    name = format.format["reference tensor"](j, i, a)
                    value = format.format["floating point"](term.A0(i, a))
                    declarations += [Declaration(name, value)]
        return declarations

    def __compute_geometry_tensor(self, format, gK_used, cK_used):
        "Precompute geometry tensor according to given format."
        debug("Generating code for geometry tensor", 1)
        if not self.terms or format.format["geometry tensor"](0, []) == None: return []
        declarations = []
        for j in range(len(self.terms)):
            # Should be the same, so pick first
            aindices = self.terms[j].GKs[0].a.indices
            if not aindices:
                aindices = [[]]
            for a in aindices:
                # Sum factorized values
                name = format.format["geometry tensor"](j, a)
                used = name in gK_used
                value = format.format["sum"]([GK(a, format, cK_used, used) for GK in self.terms[j].GKs])
                declaration = Declaration(name, value)
                declaration.used = used
                # Add declaration
                declarations += [declaration]

        return declarations

    def __compute_element_tensor(self, format, options):
        "Precompute element tensor, including possible optimizations."
        if not self.terms or format.format["element tensor"]((0,), 0) == None: return []
        if options["optimize"]:
            rank = self.terms[0].A0.i.rank
            if rank == 2 or rank == 1:
            #if rank == 2:
                return self.__compute_element_tensor_optimized(format)
            else:
                debug("Only rank 2 tensors can currently be optimized with FErari, generating default code")
                return self.__compute_element_tensor_default(format)
        else:
            return self.__compute_element_tensor_default(format)

    def __compute_element_tensor_default(self, format):
        """Precompute element tensor without optimizations except for
        dropping multiplication with zero."""
        debug("Generating code for element tensor", 1)         
        declarations = []
        iindices = self.terms[0].A0.i.indices or [[]] # All primary ranks are equal

        # Prefetch formats to speed up code generation
        format_element_tensor  = format.format["element tensor"]
        format_geometry_tensor = format.format["geometry tensor"]
        format_sum             = format.format["sum"]
        format_subtract        = format.format["subtract"]
        format_multiply        = format.format["multiplication"]
        format_floating_point  = format.format["floating point"]

        # Generate code for geometry tensor elements
        gk_tensor = [ ( [(format_geometry_tensor(j, a), a) \
                         for a in self.__aindices(j) ], j) \
                         for j in range(len(self.terms)) ]

        # Generate code for computing the element tensor
        k = 0
        num_dropped = 0
        num_ops = 0
        zero = format_floating_point(0.0)
        for i in iindices:
            value = None
            name = format_element_tensor(i, k)
            for (gka, j) in gk_tensor:
                A0 = self.terms[j].A0
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
        self.num_ops = num_ops
        return declarations

    def __compute_element_tensor_optimized(self, format):
        "Precompute element tensor with FErari optimizations."
        debug("Generating optimized code for element tensor", 1)
        # Call FErari to do optimizations
        (declarations, self.num_ops) = optimize(self.terms, format)
        return declarations

    def __check_used(self, format, gK_used):
        """Check which declarations of gK are actually used, i.e,
        which entries of the geometry tensor that get multiplied with
        nonzero entries of the reference tensor."""
        if not self.terms or format.format["geometry tensor"](0, []) == None: return []
        iindices = self.terms[0].A0.i.indices or [[]] # All primary ranks are equal
        for i in iindices:
            for j in range(len(self.terms)):
                A0 = self.terms[j].A0
                if A0.a.indices: aindices = A0.a.indices
                else: aindices = [[]]
                for a in aindices:
                    a0 = A0(tuple(i), tuple(a))
                    gk = format.format["geometry tensor"](j, a)
                    if abs(a0) > FFC_EPSILON:
                        gK_used.add(gk)
        return

    def __check_integrals(self, sum):
        "Check that all terms have integrals."
        for p in sum.products:
            if not p.integral:
                raise FormError, (p, "Missing integral in term.")

    def __terms_to_compile(self, sum, type):
        "Count the number of terms to be computed."
        count = 0
        for p in sum.products:
            if p.integral.type == type:
                count += 1
        return count

    def __aindices(self, j):
        A0 = self.terms[j].A0
        if A0.a.indices: aindices = A0.a.indices
        else: aindices = [[]]
        return aindices
