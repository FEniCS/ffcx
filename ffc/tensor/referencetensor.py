# Copyright (C) 2004-2009 Anders Logg
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
# Modified by Garth N. Wells 2006
# Modified by Kristian B. Oelgaard, 2009.
#
# First added:  2004-11-03
# Last changed: 2011-11-28

# FFC modules.
from ffc.log import debug

# FFC tensor representation modules.
from monomialintegration import integrate
from monomialtransformation import MonomialIndex
from multiindex import create_multiindex

class ReferenceTensor:
    """
    This class represents the reference tensor for a monomial term of
    a multilinear form.
    """

    def __init__(self,
                 monomial,
                 domain_type,
                 facet0, facet1,
                 quadrature_order,
                 quadrature_rule,
                 cellname,
                 facet_cellname):
        "Create reference tensor for given monomial."

        # Compute reference tensor
        self.A0 = integrate(monomial,
                            domain_type,
                            facet0, facet1,
                            quadrature_order,
                            quadrature_rule,
                            cellname,
                            facet_cellname)

        # Extract indices
        primary_indices   = monomial.extract_unique_indices(MonomialIndex.PRIMARY)
        secondary_indices = monomial.extract_unique_indices(MonomialIndex.SECONDARY)
        internal_indices  = monomial.extract_unique_indices(MonomialIndex.INTERNAL)

        # Create multiindices
        self.primary_multi_index   = create_multiindex(primary_indices)
        self.secondary_multi_index = create_multiindex(secondary_indices)
        self.internal_multi_index  = create_multiindex(internal_indices)

        # Store monomial
        self.monomial = monomial

        debug("Primary multi index:   " + str(self.primary_multi_index))
        debug("Secondary multi index: " + str(self.secondary_multi_index))
        debug("Internal multi index:  " + str(self.internal_multi_index))
