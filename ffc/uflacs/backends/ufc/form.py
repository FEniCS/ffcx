# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

# Note: Most of the code in this file is a direct translation from the old implementation in FFC

from ffc.classname import make_integral_classname
from ffc.uflacs.backends.ufc.generator import ufc_generator, integral_name_templates, ufc_integral_types
from ffc.uflacs.backends.ufc.utils import generate_return_new, generate_return_new_switch


def create_delegate(integral_type, declname, impl):
    def _delegate(self, L, ir, parameters):
        return impl(self, L, ir, parameters, integral_type, declname)
    _delegate.__doc__ = impl.__doc__ % {"declname": declname, "integral_type": integral_type}
    return _delegate


def add_ufc_form_integral_methods(cls):
    """This function generates methods on the class it decorates,
    for each integral name template and for each integral type.

    This allows implementing e.g. create_###_integrals once in the
    decorated class as '_create_foo_integrals', and this function will
    expand that implementation into 'create_cell_integrals',
    'create_exterior_facet_integrals', etc.

    Name templates are taken from 'integral_name_templates' and 'ufc_integral_types'.
    """
    # The dummy name "foo" is chosen for familiarity for ffc developers
    dummy_integral_type = "foo"

    for template in integral_name_templates:
        implname = "_" + (template % (dummy_integral_type,))
        impl = getattr(cls, implname)
        for integral_type in ufc_integral_types:
            declname = template % (integral_type,)
            _delegate = create_delegate(integral_type, declname, impl)
            setattr(cls, declname, _delegate)
    return cls


@add_ufc_form_integral_methods
class ufc_form(ufc_generator):
    """Each function maps to a keyword in the template. See documentation of ufc_generator.

    The exceptions are functions on the form
        def _*_foo_*(self, L, ir, parameters, integral_type, declname)
    which add_ufc_form_integral_methods will duplicate for foo = each integral type.
    """
    def __init__(self):
        ufc_generator.__init__(self, "form")

    def num_coefficients(self, L, num_coefficients):
        return L.Return(num_coefficients)

    def rank(self, L, rank):
        return L.Return(rank)

    def original_coefficient_position(self, L, ir):
        i = L.Symbol("i")
        positions = ir["original_coefficient_position"]

        # Check argument
        msg = "Invalid original coefficient index."

        if positions:
            code = [
            L.If(L.GE(i, len(positions)),
                L.Throw("std::runtime_error", msg))
            ]

            # Throwing a lot into the 'typename' string here but
            # no plans for building a full C++ type system
            typename = "static const std::vector<std::size_t>"
            position = L.Symbol("position")
            initializer_list = L.VerbatimExpr("{" + ", ".join(str(i) for i in positions) + "}")
            code += [
                L.VariableDecl(typename, position, value=initializer_list),
                L.Return(position[i]),
                ]
            return code
        else:
            code = [ L.Throw("std::runtime_error", msg),
                     L.Return(i) ]
        return code

    def create_coordinate_finite_element(self, L, ir):
        classnames = ir["create_coordinate_finite_element"]
        assert len(classnames) == 1
        return generate_return_new(L, classnames[0], factory=ir["jit"])

    def create_coordinate_dofmap(self, L, ir):
        classnames = ir["create_coordinate_dofmap"]
        assert len(classnames) == 1
        return generate_return_new(L, classnames[0], factory=ir["jit"])

    def create_coordinate_mapping(self, L, ir):
        classnames = ir["create_coordinate_mapping"]
        assert len(classnames) == 1  # list of length 1 until we support multiple domains
        return generate_return_new(L, classnames[0], factory=ir["jit"])

    def create_finite_element(self, L, ir):
        i = L.Symbol("i")
        classnames = ir["create_finite_element"]
        return generate_return_new_switch(L, i, classnames, factory=ir["jit"])

    def create_dofmap(self, L, ir):
        i = L.Symbol("i")
        classnames = ir["create_dofmap"]
        return generate_return_new_switch(L, i, classnames, factory=ir["jit"])


    # This group of functions are repeated for each
    # foo_integral by add_ufc_form_integral_methods:

    def _max_foo_subdomain_id(self, L, ir, parameters, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        # e.g. max_subdomain_id = ir["max_cell_subdomain_id"]
        max_subdomain_id = ir[declname]
        return L.Return(int(max_subdomain_id))

    def _has_foo_integrals(self, L, ir, parameters, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        # e.g. has_integrals = ir["has_cell_integrals"]
        has_integrals = ir[declname]
        return L.Return(bool(has_integrals))

    def _create_foo_integral(self, L, ir, parameters, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        # e.g. subdomain_ids, classnames = ir["create_cell_integral"]
        subdomain_ids, classnames = ir[declname]
        subdomain_id = L.Symbol("subdomain_id")
        return generate_return_new_switch(L, subdomain_id, classnames,
                                          subdomain_ids, factory=ir["jit"])

    def _create_default_foo_integral(self, L, ir, parameters, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        # e.g. classname = ir["create_default_cell_integral"]
        classname = ir[declname]
        if classname is None:
            return L.Return(L.Null())
        else:
            return generate_return_new(L, classname, factory=ir["jit"])
