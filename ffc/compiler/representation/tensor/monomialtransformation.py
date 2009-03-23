"Transformation of monomial representations of UFL forms."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-03-06 -- 2009-03-06"
__copyright__ = "Copyright (C) 2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# UFL modules
from ufl import BasisFunction, Function
from ufl.indexing import FixedIndex

# FFC common modules
from ffc.common.log import ffc_assert

# FFC fem modules
from ffc.fem import create_element
from ffc.fem.finiteelement import AFFINE, CONTRAVARIANT_PIOLA, COVARIANT_PIOLA

# FFC tensor representation modules
from monomialextraction import MonomialForm

# Index counters
_current_secondary_index = 0
_current_internal_index = 0
_current_external_index = 0

def next_secondary_index():
    global _current_secondary_index
    _current_secondary_index += 1
    return _current_secondary_index - 1

def next_internal_index():
    global _current_internal_index
    _current_internal_index += 1
    return _current_internal_index - 1

def next_external_index():
    global _current_external_index
    _current_external_index += 1
    return _current_external_index - 1

def reset_indices():
    global _current_secondary_index
    global _current_internal_index
    global _current_external_index
    _current_secondary_index = 0
    _current_internal_index = 0
    _current_external_index = 0

class MonomialIndex:

    FIXED      = "fixed"      # Integer index
    PRIMARY    = "primary"    # Argument basis function index
    SECONDARY  = "secondary"  # Index appearing both inside and outside integral
    INTERNAL   = "internal"   # Index appearing only inside integral
    EXTERNAL   = "external"   # Index appearing only outside integral

    def __init__(self, index_type=None, index_range=None, index_id=None):
        self.index_type = index_type
        self.index_range = index_range
        self.index_id = index_id

    def __lt__(self, other):
        return self.index_id < other.index_id

    def __call__(self, primary=None, secondary=None, internal=None, external=None):
        "Evaluate index at current index list."

        if self.index_type == MonomialIndex.FIXED:
            return self.index_id
        elif self.index_type == MonomialIndex.PRIMARY:
            if not primary:
                raise RuntimeError, "Missing index values for primary indices."
            return primary[self.index_id]
        elif self.index_type == MonomialIndex.SECONDARY:
            if not secondary:
                raise RuntimeError, "Missing index values for secondary indices."
            return secondary[self.index_id]
        elif self.index_type == MonomialIndex.INTERNAL:
            if not internal:
                raise RuntimeError, "Missing index values for internal auxiliary indices."
            return internal[self.index_id]
        elif self.index_type == MonomialIndex.EXTERNAL:
            if not external:
                raise RuntimeError, "Missing index values for external auxiliary indices."
            return external[self.index_id]
        else:
            raise RuntimeError, "Unknown index type " + str(self.type)
        return

    def __str__(self):
        if self.index_type == MonomialIndex.FIXED:
            return str(self.index_id)
        elif self.index_type == MonomialIndex.PRIMARY:
            return "i_" + str(self.index_id)
        elif self.index_type == MonomialIndex.SECONDARY:
            return "a_" + str(self.index_id)
        elif self.index_type == MonomialIndex.INTERNAL:
            return "g_" + str(self.index_id)
        elif self.index_type == MonomialIndex.EXTERNAL:
            return "b_" + str(self.index_id)
        else:
            return "?"

class MonomialDeterminant:

    # FIXME: Handle restrictions for determinants

    def __init__(self):
        self.power = 0
        self.restriction = None

    def __str__(self):
        if self.power == 0:
            return "|det F'|"
        elif self.power == 1:
            return "|det F'| (det F')"
        else:
            return "|det F'| (det F')^%s" % str(self.power)

class MonomialCoefficient:

    def __init__(self, index, number):
        self.index = index
        self.number = number

    def __str__(self):
        return "c_" + str(self.index)

class MonomialTransform:

    J = "J"
    JINV = "JINV"

    def __init__(self, index0, index1, transform_type, restriction):
        self.index0 = index0
        self.index1 = index1
        self.transform_type = transform_type
        self.restriction = restriction

    def __str__(self):
        if self.restriction is None:
            r = ""
        else:
            r = "(%s)" % str(self.restriction)
        if self.transform_type == "J":
            return "dx_%s/dX_%s%s" % (str(self.index0), str(self.index1), r)
        else:
            return "dX_%s/dx_%s%s" % (str(self.index0), str(self.index1), r)

class MonomialBasisFunction:

    def __init__(self, element, index, components, derivatives, restriction):
        self.element = element
        self.index = index
        self.components = components
        self.derivatives = derivatives
        self.restriction = restriction

    def __str__(self):
        if len(self.components) == 0:
            c = ""
        else:
            c = "[%s]" % ", ".join(str(c) for c in self.components)
        if len(self.derivatives) == 0:
            d0 = ""
            d1 = ""
        else:
            d0 = "(" + " ".join("d/dX_%s" % str(d) for d in self.derivatives) + " "
            d1 = ")"
        if self.restriction is None:
            r = ""
        else:
            r = "(%s)" % str(self.restriction)
        v = "V_" + str(self.index)
        return d0 + v + r + c + d1

class TransformedMonomial:

    def __init__(self, monomial):

        # Reset monomial data
        self.float_value = monomial.float_value
        self.determinant = MonomialDeterminant()
        self.coefficients = []
        self.transforms = []
        self.basis_functions = []

        # Reset index counters
        reset_indices()

        # Initialize index map
        index_map = {}

        # Iterate over factors
        for f in monomial.factors:

            print "Factor:", f

            # Extract element and dimensions
            element = create_element(f.element())
            vdim = element.space_dimension()
            gdim = element.geometric_dimension()
            cdim = element.num_sub_elements()

            # Extract basis function index and coefficients
            if isinstance(f.function, BasisFunction):
                vindex = MonomialIndex(MonomialIndex.PRIMARY,
                                       index_range=range(vdim),
                                       index_id=f.function.count())
            elif isinstance(f.function, Function):
                vindex = MonomialIndex(index_range=range(vdim))
                coefficient = MonomialCoefficient(vindex, f.function.count())
                self.coefficients.append(coefficient)

            # Extract components
            components = []
            for c in f.components:
                if c in index_map:
                    index = index_map[c]
                elif isinstance(c, FixedIndex):
                    index = MonomialIndex(index_type=MonomialIndex.FIXED,
                                          index_range=[int(c)],
                                          index_id=int(c))
                else:
                    index = MonomialIndex(index_range=range(cdim))
                index_map[c] = index
                components.append(index)
            if len(components) > 1:
                raise MonomialException, "Can only handle rank 0 or rank 1 tensors."
            
            # Handle non-affine mappings (Piola)
            if len(components) > 0:

                # Get type of mapping
                mapping = element.mapping(components[0].index_range)

                # FIXME: Work in progress!

                # Add transforms where appropriate
                if mapping == CONTRAVARIANT_PIOLA:
                    # phi(x) = (det J)^{-1} J Phi(X) 
                    index0 = components[0]
                    index1 = MonomialIndex(index_range=[])
                    transform = MonomialTransform(index0, index1, MonomialTransform.J, f.restriction)
                    self.determinant.power -= 1
                elif mapping == COVARIANT_PIOLA:
                    # phi(x) = J^{-T} Phi(X)
                    index1 = components[0]                    
                    index0 = MonomialIndex(index_range=[])
                    transform = MonomialTransform(index1, index0, MonomialTransform.JINV, f.restriction)

            # Extract derivatives / transforms
            derivatives = []
            for d in f.derivatives:
                index0 = MonomialIndex(index_range=range(gdim))
                if d in index_map:
                    index1 = index_map[d]
                elif isinstance(d, FixedIndex):
                    index1 = MonomialIndex(index_type=MonomialIndex.FIXED,
                                           index_range=[int(d)],
                                           index_id=int(d))
                else:
                    index1 = MonomialIndex(index_range=range(gdim))
                index_map[d] = index1
                transform = MonomialTransform(index0, index1, MonomialTransform.JINV, f.restriction)
                self.transforms.append(transform)
                derivatives.append(index0)

            # Extract restriction
            restriction = f.restriction

            print "derivatives =", derivatives

            # Create basis function
            v = MonomialBasisFunction(element, vindex, components, derivatives, restriction)
            self.basis_functions.append(v)

        # Figure out secondary and auxiliary indices
        internal_indices = self.extract_internal_indices(None)
        external_indices = self.extract_external_indices(None)        
        for i in internal_indices + external_indices:

            # Skip already visited indices
            if not i.index_type is None:
                continue

            # Set index type and id
            num_internal = len([j for j in internal_indices if j == i])
            num_external = len([j for j in external_indices if j == i])
            if num_internal == 1 and num_external == 1:
                i.index_type = MonomialIndex.SECONDARY
                i.index_id   = next_secondary_index()
            elif num_internal == 2:
                i.index_type = MonomialIndex.INTERNAL
                i.index_id   = next_internal_index()
                print "Found internal:", i
            elif num_external == 2:
                i.index_type = MonomialIndex.EXTERNAL
                i.index_id   = next_external_index()
                print "Found external:", i, i.index_type
            else:
                raise RuntimeError, "Summation index does not appear exactly twice: " + str(i)

    def extract_internal_indices(self, index_type=None):
        "Return list of indices appearing inside integral."
        indices = []
        for v in self.basis_functions:
            indices += [v.index] + v.components + v.derivatives
        return [i for i in indices if i.index_type == index_type]

    def extract_external_indices(self, index_type=None):
        "Return list of indices appearing outside integral."
        indices = [c.index for c in self.coefficients] + \
                  [t.index0 for t in self.transforms]  + \
                  [t.index1 for t in self.transforms]
        return [i for i in indices if i.index_type == index_type]        

    def extract_indices(self, index_type=None):
        "Return all indices for monomial."
        return self.extract_internal_indices(index_type) + \
               self.extract_external_indices(index_type)

    def extract_unique_indices(self, index_type=None):
        "Return all unique indices for monomial."
        indices = []
        for index in self.extract_indices(index_type):
            if not index in indices:
                indices.append(index)
        return indices

    def __str__(self):
        factors = []
        if not self.float_value == 1.0:
            factors.append(self.float_value)
        factors.append(self.determinant)
        factors += self.coefficients
        factors += self.transforms
        return " * ".join([str(f) for f in factors]) + " | " + " * ".join([str(v) for v in self.basis_functions])
    
def transform_monomial_form(monomial_form):
    "Transform monomial form to reference element."

    # Check that we get a Form
    ffc_assert(isinstance(monomial_form, MonomialForm), "Expecting a MonomialForm.")

    # Transform each monomial
    for (integrand, measure) in monomial_form:
        for (i, monomial) in enumerate(integrand.monomials):
            integrand.monomials[i] = TransformedMonomial(monomial)
