# Copyright (C) 2016 Miklos Homolya, Jan Blechta
#
# This file is part of FFC and contains snippets originally from tsfc.
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

from __future__ import absolute_import

import numpy
import six

import coffee.base as coffee

import gem
from gem.node import traversal

from tsfc.fiatinterface import create_element
from tsfc.mixedelement import MixedElement
from tsfc.coffee import SCALAR_TYPE


class Kernel(object):
    __slots__ = ("ast", "integral_type", "oriented", "subdomain_id",
                 "coefficient_numbers", "__weakref__")
    """A compiled Kernel object.

    :kwarg ast: The COFFEE ast for the kernel.
    :kwarg integral_type: The type of integral.
    :kwarg oriented: Does the kernel require cell_orientations.
    :kwarg subdomain_id: What is the subdomain id for this kernel.
    :kwarg coefficient_numbers: A list of which coefficients from the
        form the kernel needs.
    """
    def __init__(self, ast=None, integral_type=None, oriented=False,
                 subdomain_id=None, coefficient_numbers=()):
        # Defaults
        self.ast = ast
        self.integral_type = integral_type
        self.oriented = oriented
        self.subdomain_id = subdomain_id
        self.coefficient_numbers = coefficient_numbers
        super(Kernel, self).__init__()


class KernelBuilderBase(object):
    """Helper class for building local assembly kernels."""

    def __init__(self, interior_facet=False):
        """Initialise a kernel builder.

        :arg interior_facet: kernel accesses two cells
        """
        assert isinstance(interior_facet, bool)
        self.interior_facet = interior_facet

        self.prepare = []
        self.finalise = []

        self.coefficient_map = {}

    def apply_glue(self, prepare=None, finalise=None):
        """Append glue code for operations that are not handled in the
        GEM abstraction.

        Current uses: mixed interior facet mess

        :arg prepare: code snippets to be prepended to the kernel
        :arg finalise: code snippets to be appended to the kernel
        """
        if prepare is not None:
            self.prepare.extend(prepare)
        if finalise is not None:
            self.finalise.extend(finalise)

    def construct_kernel(self, name, args, body):
        """Construct a COFFEE function declaration with the
        accumulated glue code.

        :arg name: function name
        :arg args: function argument list
        :arg body: function body (:class:`coffee.Block` node)
        :returns: :class:`coffee.FunDecl` object
        """
        assert isinstance(body, coffee.Block)
        body_ = coffee.Block(self.prepare + body.children + self.finalise)
        return coffee.FunDecl("void", name, args, body_, pred=["static", "inline"])

    @property
    def coefficient_mapper(self):
        """A function that maps :class:`ufl.Coefficient`s to GEM
        expressions."""
        return lambda coefficient: self.coefficient_map[coefficient]

    def arguments(self, arguments, indices):
        """Prepare arguments. Adds glue code for the arguments.

        :arg arguments: :class:`ufl.Argument`s
        :arg indices: GEM argument indices
        :returns: COFFEE function argument and GEM expression
                  representing the argument tensor
        """
        funarg, prepare, expressions, finalise = prepare_arguments(
            arguments, indices, interior_facet=self.interior_facet)
        self.apply_glue(prepare, finalise)
        return funarg, expressions

    def coefficients(self, coefficients, name, mode=None):
        """Prepare a coefficient. Adds glue code for the coefficient
        and adds the coefficient to the coefficient map.

        :arg coefficient: iterable of :class:`ufl.Coefficient`s
        :arg name: coefficient name
        :arg mode: see :func:`prepare_coefficient`
        :returns: COFFEE function argument for the coefficient
        """
        funarg, prepare, expressions = prepare_coefficients(
            coefficients, name, mode=mode,
            interior_facet=self.interior_facet)
        self.apply_glue(prepare)
        for i, coefficient in enumerate(coefficients):
            self.coefficient_map[coefficient] = expressions[i]
        return funarg

    def coordinates(self, coefficient, name, mode=None):
        """Prepare a coordinates. Adds glue code for the coefficient
        and adds the coefficient to the coefficient map.

        :arg coefficient: :class:`ufl.Coefficient`
        :arg name: coefficient name
        :arg mode: see :func:`prepare_coefficient`
        :returns: COFFEE function argument for the coefficient
        """
        funarg, prepare, expression = prepare_coordinates(
            coefficient, name, mode=mode,
            interior_facet=self.interior_facet)
        self.apply_glue(prepare)
        self.coefficient_map[coefficient] = expression
        return funarg


class KernelBuilder(KernelBuilderBase):
    """Helper class for building a :class:`Kernel` object."""

    def __init__(self, integral_type, subdomain_id):
        """Initialise a kernel builder."""
        super(KernelBuilder, self).__init__(integral_type.startswith("interior_facet"))

        self.kernel = Kernel(integral_type=integral_type, subdomain_id=subdomain_id)
        self.local_tensor = None
        self.coordinates_arg = None
        self.coefficient_args = []
        self.coefficient_split = {}

    def set_arguments(self, arguments, indices):
        """Process arguments.

        :arg arguments: :class:`ufl.Argument`s
        :arg indices: GEM argument indices
        :returns: GEM expression representing the return variable
        """
        self.local_tensor, expressions = self.arguments(arguments, indices)
        return expressions

    def set_coordinates(self, coefficient, name, mode=None):
        """Prepare the coordinate field.

        :arg coefficient: :class:`ufl.Coefficient`
        :arg name: coordinate coefficient name
        :arg mode: see :func:`prepare_coefficient`
        """
        self.coordinates_arg = self.coordinates(coefficient, name, mode)

    def set_coefficients(self, integral_data, form_data):
        """Prepare the coefficients of the form.

        :arg integral_data: UFL integral data
        :arg form_data: UFL form data
        """
        from ufl import Coefficient, MixedElement as ufl_MixedElement, FunctionSpace
        coefficients = []
        coefficient_numbers = []
        # enabled_coefficients is a boolean array that indicates which
        # of reduced_coefficients the integral requires.
        for i in range(len(integral_data.enabled_coefficients)):
            if integral_data.enabled_coefficients[i]:
                coefficient = form_data.reduced_coefficients[i]
                if type(coefficient.ufl_element()) == ufl_MixedElement:
                    split = [Coefficient(FunctionSpace(coefficient.ufl_domain(), element))
                             for element in coefficient.ufl_element().sub_elements()]
                    coefficients.extend(split)
                    self.coefficient_split[coefficient] = split
                else:
                    coefficients.append(coefficient)
                # This is which coefficient in the original form the
                # current coefficient is.
                # Consider f*v*dx + g*v*ds, the full form contains two
                # coefficients, but each integral only requires one.
                coefficient_numbers.append(form_data.original_coefficient_positions[i])
        self.coefficient_args.append(self.coefficients(coefficients, "w"))
        self.kernel.coefficient_numbers = tuple(coefficient_numbers)

    def require_cell_orientations(self):
        """Set that the kernel requires cell orientations."""
        self.kernel.oriented = True

    def construct_kernel(self, name, body):
        """Construct a fully built :class:`Kernel`.

        This function contains the logic for building the argument
        list for assembly kernels.

        :arg name: function name
        :arg body: function body (:class:`coffee.Block` node)
        :returns: :class:`Kernel` object
        """
        args = [self.local_tensor]
        args.extend(self.coefficient_args)
        args.append(self.coordinates_arg)
        if self.kernel.oriented:
            args.append(cell_orientations_coffee_arg)
        if self.kernel.integral_type in ["exterior_facet", "exterior_facet_vert"]:
            args.append(coffee.Decl("unsigned int",
                                    coffee.Symbol("facet", rank=(1,)),
                                    qualifiers=["const"]))
        elif self.kernel.integral_type in ["interior_facet", "interior_facet_vert"]:
            args.append(coffee.Decl("unsigned int",
                                    coffee.Symbol("facet", rank=(2,)),
                                    qualifiers=["const"]))

        self.kernel.ast = KernelBuilderBase.construct_kernel(self, name, args, body)
        return self.kernel


def prepare_coefficients(coefficients, name, mode=None, interior_facet=False):
    """Bridges the kernel interface and the GEM abstraction for
    Coefficients.  Mixed element Coefficients are rearranged here for
    interior facet integrals.

    :arg coefficient: iterable of UFL Coefficients
    :arg name: unique name to refer to the Coefficient in the kernel
    :arg mode: 'manual_loop' or 'list_tensor'; two ways to deal with
               interior facet integrals on mixed elements
    :arg interior_facet: interior facet integral?
    :returns: (funarg, prepare, expressions)
         funarg     - :class:`coffee.Decl` function argument
         prepare    - list of COFFEE nodes to be prepended to the
                      kernel body
         expressions- GEM expressions referring to the Coefficient
                      values
    """
    num_coefficients = len(coefficients)
    funarg = coffee.Decl(SCALAR_TYPE, coffee.Symbol(name),
                         pointers=[("const",), ()],
                         qualifiers=["const"])

    expressions = []
    i = gem.Index()
    for j, coefficient in enumerate(coefficients):
        fiat_element = create_element(coefficient.ufl_element())
        shape = (fiat_element.space_dimension(),)
        expression = gem.ComponentTensor(
            gem.Indexed(gem.Variable(name, (num_coefficients,) + shape),
                        (j, i)),
            (i,))
        expressions.append(expression)

    return funarg, [], expressions


def prepare_coordinates(coefficient, name, mode=None, interior_facet=False):
    """Bridges the kernel interface and the GEM abstraction for
    coordinates.

    :arg coefficient: UFL Coefficient
    :arg name: unique name to refer to the Coefficient in the kernel
    :arg mode: 'manual_loop' or 'list_tensor'; two ways to deal with
               interior facet integrals on mixed elements
    :arg interior_facet: interior facet integral?
    :returns: (funarg, prepare, expression)
         funarg     - :class:`coffee.Decl` function argument
         prepare    - list of COFFEE nodes to be prepended to the
                      kernel body
         expression - GEM expression referring to the Coefficient
                      values
    """
    fiat_element = create_element(coefficient.ufl_element())

    shape = (fiat_element.space_dimension(),)
    funarg = coffee.Decl(SCALAR_TYPE, coffee.Symbol(name),
                         pointers=[("",)],
                         qualifiers=["const"])

    i = gem.Index()
    expression = gem.ComponentTensor(
        gem.Indexed(gem.Variable(name+"_dofs", shape), (i,)),
        (i,))

    # Translate UFC/DOLFIN coords into P1 dofs
    # FIXME: Implement higher-orders, see strange scheme in
    #            dolfin/mesh/Cell.h:get_coordinate_dofs,
    #        maybe fetching P2 dofs from ufc_element will be needed!
    # FIXME: Can do it without copying? Rather offsetting expression indices?
    #        Something like
    #            expression = gem.ComponentTensor(
    #                gem.Indexed(gem.Variable(name, shape), (i+j*gdim,),),
    #                   (i*num_nodes+j,))
    #        Or using integer division/remainder on single index?
    new_name = expression.children[0].children[0].name
    gdim = coefficient.ufl_element().cell().geometric_dimension()
    assert len(shape) == 1 and shape[0] % gdim == 0
    num_nodes = shape[0] / gdim
    i, j = gem.Index(), gem.Index()
    prepare = [coffee.Decl(SCALAR_TYPE, coffee.Symbol(new_name, rank=shape)),
               coffee_for(i, gdim,
                   coffee_for(j, num_nodes,
                       coffee.Assign(
                           coffee.Symbol(new_name, rank=(i,), offset=((num_nodes, j),)),
                           coffee.Symbol(name,     rank=(j,), offset=((gdim,      i),))
                       )
                   )
               )
              ]

    return funarg, prepare, expression


def prepare_arguments(arguments, indices, interior_facet=False):
    """Bridges the kernel interface and the GEM abstraction for
    Arguments.  Vector Arguments are rearranged here for interior
    facet integrals.

    :arg arguments: UFL Arguments
    :arg indices: Argument indices
    :arg interior_facet: interior facet integral?
    :returns: (funarg, prepare, expression, finalise)
         funarg      - :class:`coffee.Decl` function argument
         prepare     - list of COFFEE nodes to be prepended to the
                       kernel body
         expressions - GEM expressions referring to the argument
                       tensor
         finalise    - list of COFFEE nodes to be appended to the
                       kernel body
    """
    elements = tuple(create_element(arg.ufl_element()) for arg in arguments)

    shape = tuple(element.space_dimension() for element in elements)
    funarg = coffee.Decl(SCALAR_TYPE, coffee.Symbol("A"), pointers=[()])

    if len(arguments) == 0:
        shape = (1,)
        indices = (0,)
    expression = gem.Indexed(gem.Variable("AA", shape), indices)

    reshape = coffee.Decl(SCALAR_TYPE,
                          coffee.Symbol("(&%s)" % expression.children[0].name,
                                        rank=shape),
                          init="*reinterpret_cast<%s (*)%s>(%s)" %
                              (SCALAR_TYPE,
                               "".join("[%s]"%i for i in shape),
                               funarg.sym.gencode()
                              )
                          )
    zero = coffee.FlatBlock("memset(%s, 0, %d * sizeof(*%s));"%
        (funarg.sym.gencode(), numpy.product(shape), funarg.sym.gencode()))
    prepare = [reshape, zero]

    return funarg, prepare, [expression], []


def coffee_for(index, extent, body):
    """Helper function to make a COFFEE loop.

    :arg index: :class:`coffee.Symbol` loop index
    :arg extent: loop extent (integer)
    :arg body: loop body (COFFEE node)
    :returns: COFFEE loop
    """
    return coffee.For(coffee.Decl("int", index, init=0),
                      coffee.Less(index, extent),
                      coffee.Incr(index, 1),
                      body)


def needs_cell_orientations(ir):
    """Does a multi-root GEM expression DAG references cell
    orientations?"""
    # cell orientation always needed in UFC
    return True


cell_orientations_coffee_arg = coffee.Decl("int", coffee.Symbol("cell_orientation"))
"""COFFEE function argument for cell orientations"""
