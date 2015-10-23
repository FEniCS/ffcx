# -*- coding: utf-8 -*-
# Copyright (C) 2009-2015 Anders Logg and Martin Sandve Alnæs
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

from uflacs.backends.ufc.generator import ufc_generator, integral_name_templates, ufc_integral_types
from uflacs.backends.ufc.utils import generate_return_new_switch

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

            # Binding variables explicitly because Python closures don't
            # capture the value of integral_type for each iteration here
            def _delegate(self, L, ir, integral_type=integral_type, declname=declname, impl=impl):
                return impl(self, L, ir, integral_type, declname)
            _delegate.__doc__ = impl.__doc__ % {"declname": declname, "integral_type": integral_type}

            setattr(cls, declname, _delegate)
    return cls

@add_ufc_form_integral_methods
class ufc_form(ufc_generator):
    def __init__(self):
        ufc_generator.__init__(self, "form")

    def topological_dimension(self, L, ir):
        "Default implementation of returning topological dimension fetched from ir."
        tdim = ir["topological_dimension"]
        return L.Return(L.LiteralInt(tdim))

    def geometric_dimension(self, L, ir):
        "Default implementation of returning geometric dimension fetched from ir."
        gdim = ir["geometric_dimension"]
        return L.Return(L.LiteralInt(gdim))

    def num_coefficients(self, L, ir):
        value = ir["num_coefficients"]
        return L.Return(L.LiteralInt(value))

    def rank(self, L, ir):
        value = ir["rank"]
        return L.Return(L.LiteralInt(value))

    def original_coefficient_position(self, L, ir):
        i = L.Symbol("i")

        positions = ir["original_coefficient_position"]

        position = L.Symbol("position")

        # Throwing a lot into the 'typename' string here but no plans for building a full C++ type system
        typename = "static const std::vector<std::size_t>"
        initializer_list = L.VerbatimExpr("{" + ", ".join(str(i) for i in positions) + "}")
        code = L.StatementList([
            L.VariableDecl(typename, position, value=initializer_list),
            L.Return(position[i]),
            ])
        return code

    def create_coordinate_finite_element(self, L, ir):
        classname = ir["create_coordinate_finite_element"] # FIXME: ffc provides element id, not classname
        return L.Return(L.New(classname))
        # TODO: Use factory functions instead, here and in all create_* functions:
        #classname = ir["coordinate_finite_element_classname"] # Not in FFC
        #factoryname = make_factory_function_name(classname)
        #return L.Return(L.Call(factoryname))

    def create_coordinate_dofmap(self, L, ir):
        classname = ir["create_coordinate_dofmap"] # FIXME: ffc provides element id, not classname
        return L.Return(L.New(classname))

    def create_finite_element(self, L, ir):
        i = L.Symbol("i")
        classnames = ir["create_finite_element"] # FIXME: ffc provides element id, not classname
        return generate_return_new_switch(L, i, classnames)

    def create_dofmap(self, L, ir):
        i = L.Symbol("i")
        classnames = ir["create_dofmap"] # FIXME: ffc provides element id, not classname
        return generate_return_new_switch(L, i, classnames)

    def create_coordinate_mapping(self, L, ir):
        classname = ir["create_coordinate_mapping"] # FIXME: ffc provides element id, not classname
        return L.Return(L.New(classname))

    def _max_foo_subdomain_id(self, L, ir, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        # e.g. max_subdomain_id = ir["max_cell_subdomain_id"]
        max_subdomain_id = ir[declname]
        return L.Return(L.LiteralInt(max_subdomain_id))

    def _has_foo_integrals(self, L, ir, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        # e.g. has_integrals = ir["has_cell_integrals"]
        has_integrals = ir[declname]
        return L.Return(L.LiteralBool(has_integrals))

    def _create_foo_integral(self, L, ir, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        subdomain_id = L.Symbol("subdomain_id")
        # e.g. classname = ir["create_cell_integral"]
        classnames = ir[declname] # FIXME: ffc provides element id, not classname
        return generate_return_new_switch(L, subdomain_id, classnames)

    def _create_default_foo_integral(self, L, ir, integral_type, declname):
        "Return implementation of ufc::form::%(declname)s()."
        # e.g. classname = ir["create_default_cell_integral"]
        classname = ir[declname] # FIXME: ffc provides element id, not classname
        if classname:
            return L.Return(L.New(classname))
        else:
            return L.Return(L.Null())
