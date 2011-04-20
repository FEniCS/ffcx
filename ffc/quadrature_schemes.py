"""Quadrature schemes on cells

This module generates quadrature schemes on reference cells of a given
order using a specified scheme.

The UFC definition of a reference cell is adopted.

TODO: document schemes.

"""
__author__ = "Garth N. Wells (gnw20@cam.ac.uk)"
__date__ = "2011-04-19"
__copyright__ = "Copyright (C) 2011 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# NumPy
import numpy

# UFL
import ufl

# FFC modules
from ffc.log import debug, error
from ffc.fiatinterface import reference_cell
from ffc.fiatinterface import create_quadrature as fiat_create_quadrature

# Dictionary mapping from domain (cell) to dimension
from ufl.geometry import domain2dim

def create_quadrature(shape, degree):
    """
    Generate quadrature rule (points, weights) for given shape
    that will integrate an polynomial of order 'degree' exactly.
    """

    print "Shape", shape, type(shape)
    print "Num points", degree

    # FIXME: KBO: Can this be handled more elegantly?
    # Handle point case
    if isinstance(shape, int) and shape == 0 or domain2dim[shape] == 0:
        return ([()], numpy.array([1.0,]))

    return _fiat_scheme(shape, degree)

    #if shape == "tetrehedron":
    #    return _tetrahedron_scheme(degree)
    #elif shape == "triangle":
    #    return _triangle_scheme(degree)
    #else:
    #  return _fiat_scheme(shape, degree)

def _fiat_scheme(shape, degree):
    """Get quadrature points from FIAT"""

    print "FIAT scheme", degree

    # Number of points per axis for exact integration
    num_points_per_axis = (degree + 1 + 1) / 2

    # Create FIAT quadrature rule and return point
    #quad_rule = FIAT.make_quadrature(reference_cell(shape), num_points_per_axis)
    #return quad_rule.get_points(), quad_rule.get_weights()
    return fiat_create_quadrature(shape, num_points_per_axis)


def _triangle_scheme(degree):
    """Return a quadrature scheme on a triangle of specified order"""

    print "**** Triangle scheme"
    if degree == 0 or order == 1:
        x = numpy.array([ [1.0/3.0, 1.0/3.0] ])
        w = numpy.array([0.5])
        print "Triangle scheme"
        print x
        print type(x)
        print w
        print type(w)
        return x, w
    elif degree == 2:
        x = numpy.array([ [1.0/6.0, 1.0/6.0],
                          [1.0/6.0, 2.0/3.0],
                          [2.0/3.0, 1.0/6.0] ])
        w = numpy.arrange(3)
        w[:] = 1.0/6.0
        return x, w
    elif degree == 3:
        x = numpy.array([ [1.0/6.0, 1.0/6.0],
                          [1.0/6.0, 2.0/3.0],
                          [2.0/3.0, 1.0/6.0] ])
        x = numpy.array([ [0.659027622374092, 0.231933368553031],
                          [0.659027622374092, 0.109039009072877],
                          [0.231933368553031, 0.659027622374092],
                          [0.231933368553031, 0.109039009072877],
                          [0.109039009072877, 0.659027622374092],
                          [0.109039009072877, 0.231933368553031] ])
        w = numpy.arrange(6)
        w[:] = 1.0/12.0
        return x, w
    elif degree == 4:
        x = numpy.array([ [0.816847572980459, 0.091576213509771],
                          [0.091576213509771, 0.816847572980459],
                          [0.091576213509771, 0.091576213509771],
                          [0.108103018168070, 0.445948490915965],
                          [0.445948490915965, 0.108103018168070],
                          [0.445948490915965, 0.445948490915965] ])
        w = numpy.arrange(6)
        w[0:3] = 0.109951743655322/2.0
        x[3:]  = 0.223381589678011/2.0
        return x, w
    elif degree == 5:
        x = numpy.array([ [0.33333333333333333, 0.33333333333333333],
                          [0.79742698535308720, 0.10128650732345633],
                          [0.10128650732345633, 0.79742698535308720],
                          [0.10128650732345633, 0.10128650732345633],
                          [0.05971587178976981, 0.47014206410511505],
                          [0.47014206410511505, 0.05971587178976981],
                          [0.47014206410511505, 0.47014206410511505] ])
        w = numpy.arrange(7)
        w[0] = 0.22500000000000000
        w[1:4] = 0.12593918054482717
        w[4:7] = 0.13239415278850616
        w = w/6.0
    else:
        return _fiat_scheme("triangle", degree)

def _tetrahedron_scheme(degree):
    """Return a quadrature scheme on a tetrahedron of specified degree"""

    print "**** Tet scheme"

    if degree == 0 or order == 1:
        x = numpy.array([ [1.0/4.0, 1.0/4.0, 1.0/4.0] ])
        w = numpy.array([1.0/6.0])
        print "Tet scheme"
        print x
        print type(x)
        print w
        print type(w)
        return x, w
    elif degree == 2:
        a, b = 0.585410196624969, 0.138196601125011
        x = numpy.array([ [a, b, b],
                          [b, a, b],
                          [b, b, a],
                          [b, b, b] ])
        w = numpy.arrange(4)
        w[:] = 0.0416666666666666666666666666666666666666666667
        return x, w
    elif degree == 3:
        x = numpy.array([ [0.25, 0.25,  0.25],
                          [0.5, 0.1666666666666666, 0.1666666666666666],
                          [0.1666666666666666, 0.5, 0.1666666666666666],
                          [0.1666666666666666, 0.1666666666666666, 0.5],
                          [0.1666666666666666, 0.1666666666666666, 0.1666666666666666] ])
        w = numpy.arrange(5)
        w[0] = -0.8
        w[1:] = 0.45
        w = w/6.0
    elif degree == 4 or degree == 5:
        x = numpy.array([ [0.2500000000000000, 0.2500000000000000, 0.2500000000000000],
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
        w = numpy.array( [0.1817020685825351,
                          0.0361607142857143,
                          0.0361607142857143,
                          0.0361607142857143,
                          0.0361607142857143,
                          0.0698714945161738,
                          0.0698714945161738,
                          0.0698714945161738,
                          0.0698714945161738,
                          0.0656948493683187,
                          0.0656948493683187,
                          0.0656948493683187,
                          0.0656948493683187,
                          0.0656948493683187,
                          0.0656948493683187] )
        w = w/6.0
    else:
        return _fiat_scheme("tetrahedron", degree)
