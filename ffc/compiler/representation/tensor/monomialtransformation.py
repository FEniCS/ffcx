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

# FFC tensor representation modules
from monomialextraction import MonomialForm

# Index counters
_current_secondary_index = 0
_current_auxiliary_index = 0

def next_secondary_index():
    global _current_secondary_index
    _current_secondary_index += 1
    return _current_secondary_index - 1

def next_auxiliary_index():
    global _current_auxiliary_index
    _current_auxiliary_index += 1
    return _current_auxiliary_index - 1

def reset_indices():
    global _current_secondary_index
    global _current_auxiliary_index
    _current_secondary_index = 0
    _current_auxiliary_index = 0

class MonomialIndex:

    PRIMARY = "primary"
    SECONDARY = "secondary"
    AUXILIARY = "auxiliary"
    FIXED = "fixed"

    def __init__(self, index_type, index_range, index_id=None):
        if index_type == MonomialIndex.SECONDARY and index_id is None:
            index_id = next_secondary_index()
        elif index_type == MonomialIndex.AUXILIARY and index_id is None:
            index_id = next_auxiliary_index()
        self.index_type = index_type
        self.index_range = index_range
        self.index_id = index_id

    def __lt__(self, other):
        return self.index_id < other.index_id

    def __str__(self):
        if self.index_type == MonomialIndex.PRIMARY:
            return "i_" + str(self.index_id)
        elif self.index_type == MonomialIndex.SECONDARY:
            return "a_" + str(self.index_id)
        elif self.index_type == MonomialIndex.AUXILIARY:
            return "b_" + str(self.index_id)
        elif self.index_type == MonomialIndex.FIXED:
            return str(self.index_id)

class MonomialRestriction:

    PLUS = "plus"
    MINUS = "minus"
    CONSTANT = "constant"

class MonomialDeterminant:

    def __init__(self):
        self.power = 0
        self.restriction = None

    def __str__(self):
        if self.power == 0:
            return ""
        elif self.power == 1:
            return "(det F')"
        else:
            return "(det F')^%s" + str(self.power)

class MonomialCoefficient:

    def __init__(self, index, number):
        self.index = index
        self.number = number

    def __str__(self):
        return "c_" + str(self.index)

class MonomialTransform:

    J = "J"
    JINV = "JINV"

    def __init__(self, index0, index1):
        self.index0 = index0
        self.index1 = index1
        self.transform_type = MonomialTransform.JINV

    def __str__(self):
        return "dX_%s/dx_%s" % (str(self.index0), str(self.index1))

class MonomialBasisFunction:

    def __init__(self, element, index, components, derivatives):
        self.element = element
        self.index = index
        self.components = components
        self.derivatives = derivatives

        # FIXME: Handle restriction
        self.restriction = None

    def __str__(self):
        c = ""
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
        v = "V_" + str(self.index)
        return d0 + v + c + d1

class TransformedMonomial:

    def __init__(self, monomial):

        # Reset monomial data
        self.float_value = 1.0
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

            # Extract element and dimensions
            element = create_element(f.element())
            vdim = element.space_dimension()
            gdim = element.geometric_dimension()
            cdim = element.num_sub_elements()

            # Extract basis function index and coefficients
            if isinstance(f.function, BasisFunction):
                vindex = MonomialIndex(MonomialIndex.PRIMARY, range(vdim), f.function.count())
            elif isinstance(f.function, Function):
                vindex = MonomialIndex(MonomialIndex.SECONDARY, range(vdim))
                coefficient = MonomialCoefficient(vindex, f.function.count())
                self.coefficients.append(coefficient)

            # Extract components
            components = []
            for c in f.components:
                if c in index_map:
                    index = index_map[c]
                elif isinstance(c, FixedIndex):
                    index = MonomialIndex(MonomialIndex.FIXED, [int(c)], int(c))
                else:
                    index = MonomialIndex(MonomialIndex.AUXILIARY, range(cdim))
                index_map[c] = index
                components.append(index)

            # Extract derivatives / transforms
            derivatives = []
            for d in f.derivatives:
                index0 = MonomialIndex(MonomialIndex.SECONDARY, range(gdim))
                if d in index_map:
                    index1 = index_map[d]
                elif isinstance(d, FixedIndex):
                    index1 = MonomialIndex(MonomialIndex.FIXED, [int(d)], int(d))
                else:
                    index1 = MonomialIndex(MonomialIndex.AUXILIARY, range(gdim))
                index_map[d] = index1
                transform = MonomialTransform(index0, index1)
                self.transforms.append(transform)
                derivatives.append(index0)

            print "derivatives =", derivatives

            # Create basis function
            v = MonomialBasisFunction(element, vindex, components, derivatives)
            self.basis_functions.append(v)

    def indices(self):
        "Return all indices for monomial."
        return [c.index for c in self.coefficients] + \
               [t.index0 for t in self.transforms] + \
               [t.index1 for t in self.transforms] + \
               [v.index for v in self.basis_functions]

    def __str__(self):
        factors = []
        if not self.float_value == 1.0:
            factors.append(float_value)
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
