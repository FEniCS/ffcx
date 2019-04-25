# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

from ffc.codegeneration import form_template as ufc_form
from ffc.codegeneration.utils import (generate_return_new,
                                      generate_return_new_switch)
from ffc.ir.representation import ufc_integral_types

# These are the method names in ufc_form that are specialized for each
# integral type
integral_name_templates = ("get_{}_integral_ids", "create_{}_integral")


def create_delegate(integral_type, declname, impl):
    def _delegate(self, L, ir, parameters):
        return impl(self, L, ir, parameters, integral_type, declname)

    _delegate.__doc__ = impl.__doc__ % {"declname": declname, "integral_type": integral_type}
    return _delegate


def add_ufc_form_integral_methods(cls):
    """Generate methods on the class decorated by this function,
    for each integral name template and for each integral type.

    This allows implementing e.g. create_###_integrals once in the
    decorated class as '_create_foo_integrals', and this function will
    expand that implementation into 'create_cell_integrals',
    'create_exterior_facet_integrals', etc.

    Name templates are taken from 'integral_name_templates' and
    'ufc_integral_types'.

    """
    # The dummy name "foo" is chosen for familiarity for ffc developers
    dummy_integral_type = "foo"

    for template in integral_name_templates:
        implname = "_" + (template.format(dummy_integral_type))
        impl = getattr(cls, implname)
        for integral_type in ufc_integral_types:
            declname = template.format(integral_type)
            _delegate = create_delegate(integral_type, declname, impl)
            setattr(cls, declname, _delegate)
    return cls


@add_ufc_form_integral_methods
class UFCForm:
    """Each function maps to a keyword in the template.

    The exceptions are functions on the form

        def _*_foo_*(self, L, ir, parameters, integral_type, declname)

    which add_ufc_form_integral_methods will duplicate for foo = each integral type.
    """

    def original_coefficient_position(self, L, ir):
        i = L.Symbol("i")
        positions = ir.original_coefficient_position

        # Check argument
        msg = "Invalid original coefficient index."
        if positions:
            code = [L.If(L.GE(i, len(positions)), [L.Comment(msg), L.Return(-1)])]
            position = L.Symbol("position")
            code += [
                L.ArrayDecl("static const int64_t", position, len(positions), positions),
                L.Return(position[i]),
            ]
            return code
        else:
            code = [L.Comment(msg), L.Return(-1)]
        return code

    def generate_coefficient_position_to_name_map(self, L, ir):
        """Generate code that maps name to number."""
        cnames = ir.coefficient_names
        assert ir.num_coefficients == len(cnames)
        names = L.Symbol("names")
        if (len(cnames) == 0):
            code = [L.Return(L.Null())]
        else:
            code = [L.ArrayDecl("static const char*", names, len(cnames), cnames)]
            code += [L.Return(names)]
        return L.StatementList(code)

    def create_coordinate_finite_element(self, L, ir):
        classnames = ir.create_coordinate_finite_element
        assert len(classnames) == 1
        return generate_return_new(L, classnames[0])

    def coordinate_finite_element_declaration(self, L, ir):
        classname = ir.create_coordinate_finite_element
        code = "ufc_finite_element* create_{name}(void);\n".format(name=classname[0])
        return code

    def create_coordinate_dofmap(self, L, ir):
        classnames = ir.create_coordinate_dofmap
        assert len(classnames) == 1
        return generate_return_new(L, classnames[0])

    def coordinate_dofmap_declaration(self, L, ir):
        classname = ir.create_coordinate_dofmap
        code = "ufc_dofmap* create_{name}(void);\n".format(name=classname[0])
        return code

    def create_coordinate_mapping(self, L, ir):
        classnames = ir.create_coordinate_mapping
        # list of length 1 until we support multiple domains
        assert len(classnames) == 1
        return generate_return_new(L, classnames[0])

    def coordinate_mapping_declaration(self, L, ir):
        classname = ir.create_coordinate_mapping
        code = "ufc_coordinate_mapping* create_{name}(void);\n".format(name=classname[0])
        return code

    def create_finite_element(self, L, ir):
        i = L.Symbol("i")
        classnames = ir.create_finite_element
        return generate_return_new_switch(L, i, classnames)

    def finite_element_declaration(self, L, ir):
        classnames = set(ir.create_finite_element)
        code = ""
        for name in classnames:
            code += "ufc_finite_element* create_{name}(void);\n".format(name=name)
        return code

    def create_dofmap(self, L, ir):
        i = L.Symbol("i")
        classnames = ir.create_dofmap
        return generate_return_new_switch(L, i, classnames)

    def dofmap_declaration(self, L, ir):
        classnames = set(ir.create_dofmap)
        code = ""
        for name in classnames:
            code += "ufc_dofmap* create_{name}(void);\n".format(name=name)
        return code

    # This group of functions are repeated for each foo_integral by
    # add_ufc_form_integral_methods:

    def _get_foo_integral_ids(self, L, ir, parameters, integral_type, declname):
        """Return implementation of ufc::form::%(declname)s()."""
        code = []
        ids = L.Symbol("ids")
        for i, v in enumerate(getattr(ir, declname)[0]):
            code += [L.Assign(ids[i], v)]
        code += [L.Return()]
        return L.StatementList(code)

    def _create_foo_integral(self, L, ir, parameters, integral_type, declname):
        """Return implementation of ufc::form::%(declname)s()."""
        # e.g. subdomain_ids, classnames = ir.create_cell_integral
        subdomain_ids, classnames = getattr(ir, declname)
        subdomain_id = L.Symbol("subdomain_id")
        return generate_return_new_switch(L, subdomain_id, classnames, subdomain_ids)


def generator(ir, parameters):
    """Generate UFC code for a form"""

    factory_name = ir.classname

    d = {}
    d["factory_name"] = factory_name
    d["signature"] = "\"{}\"".format(ir.signature)
    d["rank"] = ir.rank
    d["num_coefficients"] = ir.num_coefficients

    d["num_cell_integrals"] = len(ir.create_cell_integral[0])
    d["num_exterior_facet_integrals"] = len(ir.create_exterior_facet_integral[0])
    d["num_interior_facet_integrals"] = len(ir.create_interior_facet_integral[0])
    d["num_vertex_integrals"] = len(ir.create_vertex_integral[0])
    d["num_custom_integrals"] = len(ir.create_custom_integral[0])

    import ffc.codegeneration.C.cnodes as L
    generator = UFCForm()

    statements = generator.original_coefficient_position(L, ir)
    d["original_coefficient_position"] = L.StatementList(statements)

    d["coefficient_name_map"] = generator.generate_coefficient_position_to_name_map(L, ir)

    d["create_coordinate_finite_element"] = generator.create_coordinate_finite_element(L, ir)
    d["coordinate_finite_element_declaration"] = generator.coordinate_finite_element_declaration(L, ir)
    d["create_coordinate_dofmap"] = generator.create_coordinate_dofmap(L, ir)
    d["coordinate_dofmap_declaration"] = generator.coordinate_dofmap_declaration(L, ir)
    d["create_coordinate_mapping"] = generator.create_coordinate_mapping(L, ir)
    d["coordinate_mapping_declaration"] = generator.coordinate_mapping_declaration(L, ir)
    d["create_finite_element"] = generator.create_finite_element(L, ir)
    d["finite_element_declaration"] = generator.finite_element_declaration(L, ir)
    d["create_dofmap"] = generator.create_dofmap(L, ir)
    d["dofmap_declaration"] = generator.dofmap_declaration(L, ir)

    d["get_cell_integral_ids"] = generator.get_cell_integral_ids(L, ir, parameters)
    d["get_exterior_facet_integral_ids"] = generator.get_exterior_facet_integral_ids(L, ir, parameters)
    d["get_interior_facet_integral_ids"] = generator.get_interior_facet_integral_ids(L, ir, parameters)
    d["get_vertex_integral_ids"] = generator.get_vertex_integral_ids(L, ir, parameters)
    d["get_custom_integral_ids"] = generator.get_custom_integral_ids(L, ir, parameters)

    d["create_cell_integral"] = generator.create_cell_integral(L, ir, parameters)
    d["create_interior_facet_integral"] = generator.create_interior_facet_integral(
        L, ir, parameters)
    d["create_exterior_facet_integral"] = generator.create_exterior_facet_integral(
        L, ir, parameters)
    d["create_vertex_integral"] = generator.create_vertex_integral(L, ir, parameters)
    d["create_custom_integral"] = generator.create_custom_integral(L, ir, parameters)

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fields = [fname for _, fname, _, _ in Formatter().parse(ufc_form.factory) if fname]
    assert set(fields) == set(d.keys()), "Mismatch between keys in template and in formattting dict"

    # Format implementation code
    implementation = ufc_form.factory.format_map(d)

    # Format declaration
    declaration = ufc_form.declaration.format(factory_name=factory_name)

    return declaration, implementation
