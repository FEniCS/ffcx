#!/usr/bin/env py.test
# -*- coding: utf-8 -*-

"""
Tests generating code for the different ufc integral types,
checking that the variations between them are correctly
reflected in the code.
"""

from ufl import *
#from ffc.uflacs.backends.? import ?

def test_cell_integral_body(gtest):
    """
      class cell_integral
        virtual void tabulate_tensor(double* A,
                                     const double * const * w,
                                     const double* coordinate_dofs,
                                     int cell_orientation) const = 0;
    """

    pre = """
    """

    post = """
    """

    code = """
    """

    gtest.add(pre + code + post)


def test_exterior_facet_integral_body(gtest):
    """
      class exterior_facet_integral
        virtual void tabulate_tensor(double* A,
                                     const double * const * w,
                                     const double* coordinate_dofs,
                                     std::size_t facet) const = 0;
    """

    pre = """
    """

    post = """
    """

    code = """
    """

    gtest.add(pre + code + post)


def test_interior_facet_integral_body(gtest):
    """
      class interior_facet_integral
        virtual void tabulate_tensor(double* A,
                                     const double * const * w,
                                     const double* coordinate_dofs_0,
                                     const double* coordinate_dofs_1,
                                     std::size_t facet_0,
                                     std::size_t facet_1) const = 0;
    """

    pre = """
    """

    post = """
    """

    code = """
    """

    gtest.add(pre + code + post)


def test_vertex_integral_body(gtest):
    """
      class vertex_integral
        virtual void tabulate_tensor(double* A,
                                     const double * const * w,
                                     const double* coordinate_dofs,
                                     std::size_t vertex) const = 0;
    """

    pre = """
    """

    post = """
    """

    code = """
    """

    gtest.add(pre + code + post)


def test_custom_integral_body(gtest):
    """
      class custom_integral
        virtual void tabulate_tensor(double* A,
                                     const double * const * w,
                                     const double* coordinate_dofs,
                                     std::size_t num_quadrature_points,
                                     const double* quadrature_points,
                                     const double* quadrature_weights) const = 0;
    """

    pre = """
    """

    post = """
    """

    code = """
    """

    gtest.add(pre + code + post)
