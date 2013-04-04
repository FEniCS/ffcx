# Copyright (C) 2012 Benjamin Kehlet
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
# Modified by Marie E. Rognes, 2012
#
# First added:  2012-08-15
# Last changed: 2012-09-07


from FIAT import finite_element, polynomial_set, dual_set, functional, reference_element
import time_elements_ext as ext
import numpy

class TimeElementDualSet(dual_set.DualSet):
  """. """
  def __init__(self, family, degree):
    assert(family == "Lobatto" or family == "Radau"), \
        "Unknown time element '%s'" % family
    if family == "Lobatto" :
        assert (degree > 0), "Lobatto not defined for degree < 1!"
    else :
        assert(degree >= 0), "Degree must be >= 0"

    # ids is a map from mesh entity (numbers) to dof numbers
    ids = {}

    # dofs is a list of functionals
    dofs = []

    # Only defined in 1D (on an inteval)
    cell = reference_element.UFCInterval()

    self.coords = (ext.compute_lobatto_points(degree) if family == "Lobatto"
                   else ext.compute_radau_points(degree))
    points = [(c,) for c in self.coords]

    # Create dofs from points
    dofs = [functional.PointEvaluation(cell, point)
            for point in points]

    # Create ids
    if family == "Lobatto":
      ids[0] = {0: [0], 1: [len(points)-1]}
      ids[1] = {0: range(1, len(points)-1)}
    elif family == "Radau":
      ids[0] = {0: [], 1: []}
      ids[1] = {0: range(len(points))} # Treat all Radau points as internal
    else:
      error("Undefined family: %s" % family)

    # Initialize dual set
    dual_set.DualSet.__init__(self, dofs, cell, ids)

class TimeElement(finite_element.FiniteElement):
  """."""
  def __init__(self, family, degree):
    "Create time element with given (polynomial degree)."

    # Only defined in 1D (on an inteval)
    cell = reference_element.UFCInterval()

    # Initialize polynomial space of degree 'degree'
    polynomial_space = polynomial_set.ONPolynomialSet(cell, degree)

    # Create dual (degrees of freedom)
    dual = TimeElementDualSet(family, degree)

    # Initialize super class
    finite_element.FiniteElement.__init__(self,
                                          polynomial_space,
                                          dual,
                                          degree
                                          )

  def compute_quadrature_weights(self) :
    """Compute the quadrature weights by solving a linear system of equations
    for exact integration of polynomials. We compute the integrals over
    [-1,1] of the Legendre polynomials of degree <= n - 1; These integrals
    are all zero, except for the integral of P0 which is 2.

    This requires that the n-point quadrature rule is exact at least for
    polynomials of degree n-1."""

    n = len(self.dual.coords)

    # Special case n = 0
    if n == 0 :
      weights[0] = 2.0;
      return weights

    # Initialize linear system
    A = ext.compute_legendre_coeffs(self.dual.coords)

    b = numpy.zeros(n)
    b[0] = 2.0

    weights = numpy.linalg.solve(A, b)

    # Weights are computed on interval [-1, 1]. Scale to reference interval
    return weights/2.0


class LobattoElement(TimeElement):
  """."""
  def __init__(self, degree):
      "Create Lobatto element with given (polynomial degree)."
      TimeElement.__init__(self, "Lobatto", degree)

class RadauElement(TimeElement):
  """."""
  def __init__(self, degree):
      "Create Radau element with given (polynomial degree)."
      TimeElement.__init__(self, "Radau", degree)

