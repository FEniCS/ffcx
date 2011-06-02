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
# Modified by Marie E. Rognes, 2007
# Modified by Kristian B. Oelgaard, 2009
#
# First added:  2004-11-03
# Last changed: 2009-12-21

# FFC modules.
from ffc.log import debug

# FFC tensor representation modules.
from monomialtransformation import MonomialIndex
from multiindex import create_multiindex

class GeometryTensor:
    """
    This class represents the geometry tensor for a monomial term of a
    multilinear form.
    """

    def __init__(self, monomial):
        "Create geometry tensor for given monomial."

        # Save monomial data
        self.determinant = monomial.determinant
        self.coefficients = monomial.coefficients
        self.transforms = monomial.transforms

        # Extract indices
        secondary_indices = monomial.extract_unique_indices(MonomialIndex.SECONDARY)
        external_indices  = monomial.extract_unique_indices(MonomialIndex.EXTERNAL)

        # Create multiindices
        self.secondary_multi_index = create_multiindex(secondary_indices)
        self.external_multi_index  = create_multiindex(external_indices)

        debug("Secondary multi index: " + str(self.secondary_multi_index))
        debug("External multi index:  " + str(self.external_multi_index))
