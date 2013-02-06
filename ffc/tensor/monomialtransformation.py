"Transformation of monomial representations of UFL forms."

# Copyright (C) 2009 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Kristian B. Oelgaard, 2009
# Modified by Marie E. Rognes, 2010
#
# First added:  2009-03-06
# Last changed: 2010-02-17

# UFL modules
from ufl.classes import Argument
from ufl.classes import Coefficient
from ufl.classes import FixedIndex
from ufl.permutation import build_component_numbering

# FFC modules
from ffc.log import info, error, ffc_assert
from ffc.fiatinterface import create_element
from ffc.utils import all_equal
from ffc.representationutils import transform_component

# FFC tensor representation modules
from ffc.tensor.monomialextraction import MonomialForm
from ffc.tensor.monomialextraction import MonomialException

def transform_monomial_form(monomial_form):
    "Transform monomial form to reference element."

    info("Transforming monomial form to reference element")

    # Check that we get a monomial form
    ffc_assert(isinstance(monomial_form, MonomialForm),
               "Expecting a MonomialForm.")

    # Note that we check if each monomial has been transformed before
    # and if so we leave it untouched. This is to prevent repeated
    # transformation (which fails) which may sometimes happen as a
    # result of extracted integrands being cached by the monomial
    # extraction.

    # Transform each integral
    for (integrand, measure) in monomial_form:
        for (i, monomial) in enumerate(integrand.monomials):
            if not isinstance(monomial, TransformedMonomial):
                integrand.monomials[i] = TransformedMonomial(monomial)

class MonomialIndex:
    """
    This class represents a monomial index. Each index has a type,
    a range and a unique id. Valid index types are listed below.
    """

    FIXED     = "fixed"      # Integer index
    PRIMARY   = "primary"    # Argument basis function index
    SECONDARY = "secondary"  # Index appearing both inside and outside integral
    INTERNAL  = "internal"   # Index appearing only inside integral
    EXTERNAL  = "external"   # Index appearing only outside integral

    def __init__(self, index=None, index_type=None, index_range=None, index_id=None):
        "Create index with given type, range and id."
        if isinstance(index, MonomialIndex):
            self.index_type = index.index_type
            self.index_range = [i for i in index.index_range]
            self.index_id = index.index_id
        else:
            self.index_type = index_type
            self.index_range = index_range
            self.index_id = index_id

    def __lt__(self, other):
        "Comparison operator."
        return self.index_id < other.index_id

    def __call__(self, primary=None, secondary=None, internal=None, external=None):
        "Evaluate index at current index list."

        if self.index_type == MonomialIndex.FIXED:
            return self.index_range[0]
        elif self.index_type == MonomialIndex.PRIMARY:
            if not primary:
                error("Missing index values for primary indices.")
            return primary[self.index_id]
        elif self.index_type == MonomialIndex.SECONDARY:
            if not secondary:
                error("Missing index values for secondary indices.")
            return secondary[self.index_id]
        elif self.index_type == MonomialIndex.INTERNAL:
            if not internal:
                error("Missing index values for internal auxiliary indices.")
            return internal[self.index_id]
        elif self.index_type == MonomialIndex.EXTERNAL:
            if not external:
                error("Missing index values for external auxiliary indices.")
            return external[self.index_id]
        else:
            error("Unknown index type " + str(self.type))

    def __add__(self, offset):
        "Add offset to index range."
        index = MonomialIndex(self)
        index.index_range = [offset + i for i in index.index_range]
        return index

    def __sub__(self, offset):
        "Subtract offset from index range."
        return self + (-offset)

    def __str__(self):
        "Return informal string representation (pretty-print)."
        if self.index_type == MonomialIndex.FIXED:
            return str(self.index_range[0])
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
    "This class representes a determinant factor in a monomial."

    # FIXME: Handle restrictions for determinants

    def __init__(self):
        "Create empty monomial determinant."
        self.power = 0
        self.restriction = None

    def __str__(self):
        "Return informal string representation (pretty-print)."
        if self.power == 0:
            return "|det F'|"
        elif self.power == 1:
            return "|det F'| (det F')"
        else:
            return "|det F'| (det F')^%s" % str(self.power)

class MonomialCoefficient:
    "This class represents a coefficient in a monomial."

    def __init__(self, index, number):
        "Create monomial coefficient for given index and number."
        self.index = index
        self.number = number

    def __str__(self):
        "Return informal string representation (pretty-print)."
        return "c_" + str(self.index)

class MonomialTransform:
    "This class represents a transform (mapping derivative) in a form."

    J = "J"
    JINV = "JINV"

    def __init__(self, index0, index1, transform_type, restriction, offset):
        "Create monomial transform."

        # Set data
        self.index0 = index0
        self.index1 = index1
        self.transform_type = transform_type
        self.restriction = restriction
        self.offset = offset

        # Subtract offset for fixed indices. Note that the index subtraction
        # creates a new index instance. This is ok here since a fixed index
        # does not need to match any other index (being the same instance)
        # in index summation and index extraction.
        if index0.index_type is MonomialIndex.FIXED:
            self.index0 = index0 - offset
        if index1.index_type is MonomialIndex.FIXED:
            self.index1 = index1 - offset

    def __str__(self):
        "Return informal string representation (pretty-print)."
        if self.restriction is None:
            r = ""
        else:
            r = "(%s)" % str(self.restriction)
        if self.transform_type == "J":
            return "dx_%s/dX_%s%s" % (str(self.index0), str(self.index1), r)
        else:
            return "dX_%s/dx_%s%s" % (str(self.index0), str(self.index1), r)

class MonomialArgument:
    """
    This class represents a monomial argument, that is, a derivative of
    a scalar component of a basis function on the reference element.
    """

    def __init__(self, element, index, components, derivatives, restriction):
        "Create monomial argument."
        self.element = element
        self.index = index
        self.components = components
        self.derivatives = derivatives
        self.restriction = restriction

    def __str__(self):
        "Return informal string representation (pretty-print)."
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
    """
    This class represents a monomial form after transformation to the
    reference element.
    """

    def __init__(self, monomial):
        "Create transformed monomial from given monomial."

        # Reset monomial data
        self.float_value = monomial.float_value
        self.determinant = MonomialDeterminant()
        self.coefficients = []
        self.transforms = []
        self.arguments = []

        # Reset index counters
        _reset_indices()

        # Initialize index map
        index_map = {}

        # Iterate over factors
        for f in monomial.factors:

            # Create FIAT element
            ufl_element = f.element()
            fiat_element = create_element(ufl_element)

            # Note nifty aspect here: when gdim != tdim, it might be
            # (in particular, is for H(div)/H(curl), that the value
            # dimension is different for the physical and reference
            # elements.

            # Get number of components
            # FIXME: Can't handle tensor-valued elements: vdim = shape[0]
            shape = ufl_element.value_shape()
            assert(len(shape) <= 1), \
                "MonomialTransformation does not handle tensor-valued elements"
            if len(shape) == 0:
                vdim = 1
            else:
                vdim = shape[0]

            # Extract dimensions
            sdim = fiat_element.space_dimension()
            gdim = ufl_element.cell().geometric_dimension()
            tdim = ufl_element.cell().topological_dimension()

            # Extract basis function index and coefficients
            if isinstance(f.function, Argument):
                vindex = MonomialIndex(index_type=MonomialIndex.PRIMARY,
                                       index_range=range(sdim),
                                       index_id=f.count())

            elif isinstance(f.function, Coefficient):
                vindex = MonomialIndex(index_range=range(sdim))
                coefficient = MonomialCoefficient(vindex, f.count())
                self.coefficients.append(coefficient)

            # Extract components
            components = self._extract_components(f, index_map, vdim)

            if len(components) > 1:
                raise MonomialException, "Can only handle rank 0 or rank 1 tensors."

            # Handle non-affine mappings (Piola)
            if len(components) > 0:

                # We can only handle rank 1 elements for now
                component = components[0]

                # Get mapping (all need to be equal)
                mappings = []
                for i in component.index_range:
                    (offset, ufl_sub_element) = ufl_element.extract_component(i)
                    fiat_sub_element = create_element(ufl_sub_element)
                    mappings.extend(fiat_sub_element.mapping())
                if not all_equal(mappings):
                    raise MonomialException, ("Mappings differ: " + str(mappings))
                mapping = mappings[0]

                # Get component index relative to its sub element and its sub element
                (component_index, sub_element) = ufl_element.extract_component(component.index_range[0])

                # Get offset
                if len(component_index) == 0:
                    offset = 0
                else:
                    offset = component.index_range[0] - component_index[0]

                # MER: Need to handle mappings in special ways if gdim
                # != tdim and some Piolas are present. This could
                # probably be merged with the offset code above, but I
                # was not able to wrap my head around the offsets
                # always referring to component.index_range[0].
                if (gdim != tdim):
                    assert len(component.index_range) == 1, \
                        "Component transform not implemented for this case. Please request this feature."
                    component, offset = transform_component(component.index_range[0], offset, ufl_element)
                    component = MonomialIndex(index_type=MonomialIndex.FIXED,
                                              index_range=[component], index_id=None)
                    components = [component, ]

                # Add transforms where appropriate
                if mapping == "contravariant piola":
                    # phi(x) = (det J)^{-1} J Phi(X)
                    index0 = component
                    index1 = MonomialIndex(index_range=range(tdim)) + offset
                    transform = MonomialTransform(index0, index1, MonomialTransform.J,
                                                  f.restriction, offset)
                    self.transforms.append(transform)
                    self.determinant.power -= 1
                    components[0] = index1
                elif mapping == "covariant piola":
                    # phi(x) = J^{-T} Phi(X)
                    index0 = MonomialIndex(index_range=range(tdim)) + offset
                    index1 = component
                    transform = MonomialTransform(index0, index1, MonomialTransform.JINV,
                                                  f.restriction, offset)
                    self.transforms.append(transform)
                    components[0] = index0

            # Extract derivatives / transforms
            derivatives = []
            for d in f.derivatives:
                index0 = MonomialIndex(index_range=range(tdim))
                if d in index_map:
                    index1 = index_map[d]
                elif isinstance(d, FixedIndex):
                    index1 = MonomialIndex(index_type=MonomialIndex.FIXED,
                                           index_range=[int(d)],
                                           index_id=int(d))
                else:
                    index1 = MonomialIndex(index_range=range(gdim))
                index_map[d] = index1
                transform = MonomialTransform(index0, index1, MonomialTransform.JINV, f.restriction, 0)

                self.transforms.append(transform)
                derivatives.append(index0)

            # Extract restriction
            restriction = f.restriction

            # Create basis function
            v = MonomialArgument(ufl_element, vindex, components, derivatives, restriction)
            self.arguments.append(v)

        # Figure out secondary and auxiliary indices
        internal_indices = self._extract_internal_indices(None)
        external_indices = self._extract_external_indices(None)
        for i in internal_indices + external_indices:

            # Skip already visited indices
            if not i.index_type is None:
                continue

            # Set index type and id
            num_internal = len([j for j in internal_indices if j == i])
            num_external = len([j for j in external_indices if j == i])

            if num_internal == 1 and num_external == 1:
                i.index_type = MonomialIndex.SECONDARY
                i.index_id   = _next_secondary_index()
            elif num_internal == 2 and num_external == 0:
                i.index_type = MonomialIndex.INTERNAL
                i.index_id   = _next_internal_index()
            elif num_internal == 0 and num_external == 2:
                i.index_type = MonomialIndex.EXTERNAL
                i.index_id   = _next_external_index()
            else:
                raise Exception("Summation index does not appear exactly twice: %s" % str(i))

    def extract_unique_indices(self, index_type=None):
        "Return all unique indices for monomial w.r.t. type and id (not range)."
        indices = []
        for index in self._extract_indices(index_type):
            if not index in indices:
                indices.append(index)
        return indices

    def _extract_components(self, f, index_map, vdim):
        "Return list of components."
        components = []
        for c in f.components:
            if c in index_map:
                index = index_map[c]
            elif isinstance(c, FixedIndex):
                # Map component using component map from UFL.
                # KBO: Is this the right place to add, and do we only have
                # scalar components in the tensor representation at this stage
                # in the representation?
                comp_map, comp_num = build_component_numbering(f.element().value_shape(), f.element().symmetry())
                comp = comp_map[(int(c),)]
                index = MonomialIndex(index_type=MonomialIndex.FIXED,
                                      index_range=[comp],
                                      index_id=None)
            else:
                index = MonomialIndex(index_range=range(vdim))
            index_map[c] = index
            components.append(index)
        return components

    def _extract_internal_indices(self, index_type=None):
        "Return list of indices appearing inside integral."
        indices = []
        for v in self.arguments:
            indices += [v.index] + v.components + v.derivatives
        return [i for i in indices if i.index_type == index_type]

    def _extract_external_indices(self, index_type=None):
        "Return list of indices appearing outside integral."
        indices = [c.index for c in self.coefficients] + \
                  [t.index0 for t in self.transforms]  + \
                  [t.index1 for t in self.transforms]
        return [i for i in indices if i.index_type == index_type]

    def _extract_indices(self, index_type=None):
        "Return all indices for monomial."
        return self._extract_internal_indices(index_type) + \
               self._extract_external_indices(index_type)

    def __str__(self):
        "Return informal string representation (pretty-print)."
        factors = []
        if not self.float_value == 1.0:
            factors.append(self.float_value)
        factors.append(self.determinant)
        factors += self.coefficients
        factors += self.transforms
        return " * ".join([str(f) for f in factors]) + " | " + " * ".join([str(v) for v in self.arguments])

# Index counters
_current_secondary_index = 0
_current_internal_index = 0
_current_external_index = 0

def _next_secondary_index():
    "Return next available secondary index."
    global _current_secondary_index
    _current_secondary_index += 1
    return _current_secondary_index - 1

def _next_internal_index():
    "Return next available internal index."
    global _current_internal_index
    _current_internal_index += 1
    return _current_internal_index - 1

def _next_external_index():
    "Return next available external index."
    global _current_external_index
    _current_external_index += 1
    return _current_external_index - 1

def _reset_indices():
    "Reset all index counters."
    global _current_secondary_index
    global _current_internal_index
    global _current_external_index
    _current_secondary_index = 0
    _current_internal_index = 0
    _current_external_index = 0
