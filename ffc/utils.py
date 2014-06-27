# Copyright (C) 2005-2010 Anders Logg
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
#
# First added:  2005-02-04
# Last changed: 2014-04-02

# Python modules.
import operator
import functools
import itertools

# FFC modules.
from log import error

def product(sequence):
    "Return the product of all elements in a sequence."
    # Copied from UFL
    return functools.reduce(operator.__mul__, sequence, 1)

def all_equal(sequence):
    "Check that all items in list are equal."
    return sequence[:-1] == sequence[1:]

def pick_first(sequence):
    "Check that all values are equal and return the value."
    if not all_equal(sequence):
        error("Values differ: " + str(sequence))
    return sequence[0]

def listcopy(sequence):
    """Create a copy of the list, calling the copy constructor on each
    object in the list (problems when using copy.deepcopy)."""
    if not sequence:
        return []
    else:
        return [object.__class__(object) for object in sequence]

def compute_permutations(k, n, skip = []):
   """Compute all permutations of k elements from (0, n) in rising order.
   Any elements that are contained in the list skip are not included."""
   if k == 1:
       return [(i,) for i in range(n) if not i in skip]
   pp = compute_permutations(k - 1, n, skip)
   permutations = []
   for i in range(n):
       if i in skip:
           continue
       for p in pp:
           if i < p[0]:
               permutations += [(i, ) + p]
   return permutations

def compute_derivative_tuples(n, gdim):
    """Compute the list of all derivative tuples for derivatives of
    given total order n and given geometric dimension gdim. This
    function returns two lists. The first is a list of tuples, where
    each tuple of length n specifies the coordinate directions of the
    n derivatives. The second is a corresponding list of tuples, where
    each tuple of length gdim specifies the number of derivatives in
    each direction. Both lists have length gdim^n and are ordered as
    expected by the UFC function tabulate_basis_derivatives.

    Example: If n = 2 and gdim = 3, then the nice tuples are

      (0, 0)  <-->  (2, 0, 0)  <-->  d^2/dxdx
      (0, 1)  <-->  (1, 1, 0)  <-->  d^2/dxdy
      (0, 2)  <-->  (1, 0, 1)  <-->  d^2/dxdz
      (1, 0)  <-->  (1, 1, 0)  <-->  d^2/dydx
      (1, 1)  <-->  (0, 2, 0)  <-->  d^2/dydy
      (1, 2)  <-->  (0, 1, 1)  <-->  d^2/dydz
      (2, 0)  <-->  (1, 0, 1)  <-->  d^2/dzdx
      (2, 1)  <-->  (0, 1, 1)  <-->  d^2/dzdy
      (2, 2)  <-->  (0, 0, 2)  <-->  d^2/dzdz
    """

    # Create list of derivatives (note that we have d^n derivatives)
    deriv_tuples = [d for d in itertools.product(*(n*[list(range(0, gdim))]))]

    # Translate from list of derivative tuples to list of tuples
    # expressing the number of derivatives in each dimension...
    _deriv_tuples = [tuple(len([_d for _d in d if _d == i]) for i in range(gdim))
                     for d in deriv_tuples]

    return deriv_tuples, _deriv_tuples
