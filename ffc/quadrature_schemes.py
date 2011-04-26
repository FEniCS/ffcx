"""Quadrature schemes on cells

This module generates quadrature schemes on reference cells that integrate
exactly a polynomial of a given degree using a specified scheme. The
UFC definition of a reference cell is used.

Scheme options are:

  scheme="default"

  scheme="canonical" (collapsed Gauss scheme supplied by FIAT)

Background on the schemes:

  Keast rules for tetrahedra:
    Keast, P. Moderate-degree tetrahedral quadrature formulas, Computer
    Methods in Applied Mechanics and Engineering 55(3):339-348, 1986.
    http://dx.doi.org/10.1016/0045-7825(86)90059-9
"""
__author__ = "Garth N. Wells (gnw20@cam.ac.uk)"
__date__ = "2011-04-19"
__copyright__ = "Copyright (C) 2011 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# NumPy
from numpy import array, arange, float64

# UFL
import ufl

# FFC modules
from ffc.log import debug, error
from ffc.fiatinterface import reference_cell
from ffc.fiatinterface import create_quadrature as fiat_create_quadrature

# Dictionary mapping from domain (cell) to dimension
from ufl.geometry import domain2dim

def create_quadrature(shape, degree, scheme="default"):
    """
    Generate quadrature rule (points, weights) for given shape
    that will integrate an polynomial of order 'degree' exactly.
    """

    # FIXME: KBO: Can this be handled more elegantly?
    # Handle point case
    if isinstance(shape, int) and shape == 0 or domain2dim[shape] == 0:
        return ([()], array([1.0,]))

    if scheme == "default":
        if shape == "tetrahedron":
            return _tetrahedron_scheme(degree)
        elif shape == "triangle":
            return _triangle_scheme(degree)
        else:
          return _fiat_scheme(shape, degree)
    elif scheme == "canonical":
        return _fiat_scheme(shape, degree)
    else:
        error("Unkown quadrature scheme: %s." % scheme)

def _fiat_scheme(shape, degree):
    """Get quadrature scheme from FIAT interface"""

    # Number of points per axis for exact integration
    num_points_per_axis = (degree + 1 + 1) / 2

    # Create and return FIAT quadrature rulet
    return fiat_create_quadrature(shape, num_points_per_axis)


def _triangle_scheme(degree):
    """Return a quadrature scheme on a triangle of specified order. Falls
    back on canonical rule for higher orders."""

    if degree == 0 or degree == 1:
        # Scheme from Zienkiewicz and Taylor, 1 point, degree of precision 1
        x = array([ [1.0/3.0, 1.0/3.0] ])
        w = array([0.5])
    elif degree == 2:
        # Scheme from Strang and Fix, 3 points, degree of precision 2
        x = array([ [1.0/6.0, 1.0/6.0],
                          [1.0/6.0, 2.0/3.0],
                          [2.0/3.0, 1.0/6.0] ])
        w = arange(3, dtype=float64)
        w[:] = 1.0/6.0
    elif degree == 3:
        # Scheme from Strang and Fix, 6 points, degree of precision 3
        x = array([ [0.659027622374092, 0.231933368553031],
                    [0.659027622374092, 0.109039009072877],
                    [0.231933368553031, 0.659027622374092],
                    [0.231933368553031, 0.109039009072877],
                    [0.109039009072877, 0.659027622374092],
                    [0.109039009072877, 0.231933368553031] ])
        w = arange(6, dtype=float64)
        w[:] = 1.0/12.0
    elif degree == 4:
        # Scheme from Strang and Fix, 6 points, degree of precision 4
        x = array([ [0.816847572980459, 0.091576213509771],
                    [0.091576213509771, 0.816847572980459],
                    [0.091576213509771, 0.091576213509771],
                    [0.108103018168070, 0.445948490915965],
                    [0.445948490915965, 0.108103018168070],
                    [0.445948490915965, 0.445948490915965] ])
        w = arange(6, dtype=float64)
        w[0:3] = 0.109951743655322
        w[3:6] = 0.223381589678011
        w = w/2.0
    elif degree == 5:
        # Scheme from Strang and Fix, 7 points, degree of precision 5
        x = array([ [0.33333333333333333, 0.33333333333333333],
                    [0.79742698535308720, 0.10128650732345633],
                    [0.10128650732345633, 0.79742698535308720],
                    [0.10128650732345633, 0.10128650732345633],
                    [0.05971587178976981, 0.47014206410511505],
                    [0.47014206410511505, 0.05971587178976981],
                    [0.47014206410511505, 0.47014206410511505] ])
        w = arange(7, dtype=float64)
        w[0]   = 0.22500000000000000
        w[1:4] = 0.12593918054482717
        w[4:7] = 0.13239415278850616
        w = w/2.0
    elif degree == 5:
        # Scheme from Strang and Fix, 12 points, degree of precision 6
        x = array([ [0.873821971016996, 0.063089014491502],
                    [0.063089014491502, 0.873821971016996],
                    [0.063089014491502, 0.063089014491502],
                    [0.501426509658179, 0.249286745170910],
                    [0.249286745170910, 0.501426509658179],
                    [0.249286745170910, 0.249286745170910],
                    [0.636502499121399, 0.310352451033785],
                    [0.636502499121399, 0.053145049844816],
                    [0.310352451033785, 0.636502499121399],
                    [0.310352451033785, 0.053145049844816],
                    [0.053145049844816, 0.636502499121399],
                    [0.053145049844816, 0.310352451033785] ])
        w = arange(12, dtype=float64)
        w[0:3]  = 0.050844906370207
        w[3:6]  = 0.116786275726379
        w[6:12] = 0.082851075618374
        w = w/2.0
    else:
        # Get canonical scheme
        x, w = _fiat_scheme("triangle", degree)

    # Return scheme
    return x, w


def _tetrahedron_scheme(degree):
    """Return a quadrature scheme on a tetrahedron of specified
    degree. Falls back on canonical rule for higher orders"""

    if degree == 0 or degree == 1:
        # Scheme from Zienkiewicz and Taylor, 1 point, degree of precision 1
        x = array([ [1.0/4.0, 1.0/4.0, 1.0/4.0] ])
        w = array([1.0/6.0])
    elif degree == 2:
        # Scheme from Zienkiewicz and Taylor, 4 points, degree of precision 2
        a, b = 0.585410196624969, 0.138196601125011
        x = array([ [a, b, b],
                    [b, a, b],
                    [b, b, a],
                    [b, b, b] ])
        w = arange(4, dtype=float64)
        w[:] = 1.0/24.0
    elif degree == 3:
        # Scheme from Zienkiewicz and Taylor, 5 points, degree of precision 3
        # Note: this scheme has a negative weight
        x = array([ [0.2500000000000000, 0.2500000000000000, 0.2500000000000000],
                    [0.5000000000000000, 0.1666666666666666, 0.1666666666666666],
                    [0.1666666666666666, 0.5000000000000000, 0.1666666666666666],
                    [0.1666666666666666, 0.1666666666666666, 0.5000000000000000],
                    [0.1666666666666666, 0.1666666666666666, 0.1666666666666666] ])
        w = arange(5, dtype=float64)
        w[0]   = -0.8
        w[1:5] =  0.45
        w = w/6.0
    elif degree == 4:
        # Keast rule, 14 points, degree of precision 4
        # Values taken from http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
        # (KEAST5)
        x = array([ [0.0000000000000000, 0.5000000000000000, 0.5000000000000000],
                    [0.5000000000000000, 0.0000000000000000, 0.5000000000000000],
                    [0.5000000000000000, 0.5000000000000000, 0.0000000000000000],
                    [0.5000000000000000, 0.0000000000000000, 0.0000000000000000],
                    [0.0000000000000000, 0.5000000000000000, 0.0000000000000000],
                    [0.0000000000000000, 0.0000000000000000, 0.5000000000000000],
                    [0.6984197043243866, 0.1005267652252045, 0.1005267652252045],
                    [0.1005267652252045, 0.1005267652252045, 0.1005267652252045],
                    [0.1005267652252045, 0.1005267652252045, 0.6984197043243866],
                    [0.1005267652252045, 0.6984197043243866, 0.1005267652252045],
                    [0.0568813795204234, 0.3143728734931922, 0.3143728734931922],
                    [0.3143728734931922, 0.3143728734931922, 0.3143728734931922],
                    [0.3143728734931922, 0.3143728734931922, 0.0568813795204234],
                    [0.3143728734931922, 0.0568813795204234, 0.3143728734931922] ])
        w = arange(14, dtype=float64)
        w[0:6]   = 0.0190476190476190
        w[6:10]  = 0.0885898247429807
        w[10:14] = 0.1328387466855907
        w = w/6.0
    elif degree == 5:
        # Keast rule, 15 points, degree of precision 5
        # Values taken from http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
        # (KEAST6)
        x = array([ [0.2500000000000000, 0.2500000000000000, 0.2500000000000000],
                    [0.0000000000000000, 0.3333333333333333, 0.3333333333333333],
                    [0.3333333333333333, 0.3333333333333333, 0.3333333333333333],
                    [0.3333333333333333, 0.3333333333333333, 0.0000000000000000],
                    [0.3333333333333333, 0.0000000000000000, 0.3333333333333333],
                    [0.7272727272727273, 0.0909090909090909, 0.0909090909090909],
                    [0.0909090909090909, 0.0909090909090909, 0.0909090909090909],
                    [0.0909090909090909, 0.0909090909090909, 0.7272727272727273],
                    [0.0909090909090909, 0.7272727272727273, 0.0909090909090909],
                    [0.4334498464263357, 0.0665501535736643, 0.0665501535736643],
                    [0.0665501535736643, 0.4334498464263357, 0.0665501535736643],
                    [0.0665501535736643, 0.0665501535736643, 0.4334498464263357],
                    [0.0665501535736643, 0.4334498464263357, 0.4334498464263357],
                    [0.4334498464263357, 0.0665501535736643, 0.4334498464263357],
                    [0.4334498464263357, 0.4334498464263357, 0.0665501535736643] ])
        w = arange(15, dtype=float64)
        w[0]    = 0.1817020685825351
        w[1:5]  = 0.0361607142857143
        w[5:9]  = 0.0698714945161738
        w[9:15] = 0.0656948493683187
        w = w/6.0
    elif degree == 6:
        # Keast rule, 24 points, degree of precision 6
        # Values taken from http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
        # (KEAST7)
        x = array([ [0.3561913862225449, 0.2146028712591517, 0.2146028712591517],
                    [0.2146028712591517, 0.2146028712591517, 0.2146028712591517],
                    [0.2146028712591517, 0.2146028712591517, 0.3561913862225449],
                    [0.2146028712591517, 0.3561913862225449, 0.2146028712591517],
                    [0.8779781243961660, 0.0406739585346113, 0.0406739585346113],
                    [0.0406739585346113, 0.0406739585346113, 0.0406739585346113],
                    [0.0406739585346113, 0.0406739585346113, 0.8779781243961660],
                    [0.0406739585346113, 0.8779781243961660, 0.0406739585346113],
                    [0.0329863295731731, 0.3223378901422757, 0.3223378901422757],
                    [0.3223378901422757, 0.3223378901422757, 0.3223378901422757],
                    [0.3223378901422757, 0.3223378901422757, 0.0329863295731731],
                    [0.3223378901422757, 0.0329863295731731, 0.3223378901422757],
                    [0.2696723314583159, 0.0636610018750175, 0.0636610018750175],
                    [0.0636610018750175, 0.2696723314583159, 0.0636610018750175],
                    [0.0636610018750175, 0.0636610018750175, 0.2696723314583159],
                    [0.6030056647916491, 0.0636610018750175, 0.0636610018750175],
                    [0.0636610018750175, 0.6030056647916491, 0.0636610018750175],
                    [0.0636610018750175, 0.0636610018750175, 0.6030056647916491],
                    [0.0636610018750175, 0.2696723314583159, 0.6030056647916491],
                    [0.2696723314583159, 0.6030056647916491, 0.0636610018750175],
                    [0.6030056647916491, 0.0636610018750175, 0.2696723314583159],
                    [0.0636610018750175, 0.6030056647916491, 0.2696723314583159],
                    [0.2696723314583159, 0.0636610018750175, 0.6030056647916491],
                    [0.6030056647916491, 0.2696723314583159, 0.0636610018750175] ])
        w = arange(24, dtype=float64)
        w[0:4]   = 0.0399227502581679
        w[4:8]   = 0.0100772110553207
        w[8:12]  = 0.0553571815436544
        w[12:24] = 0.0482142857142857
        w = w/6.0
    else:
        # Get canonical scheme
        x, w =  _fiat_scheme("tetrahedron", degree)

    # Return scheme
    return x, w
