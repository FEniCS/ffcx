# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

# Note: Most of the code in this file is a direct translation from the
# old implementation in FFC

import logging

from ffcx.codegeneration import form_template

logger = logging.getLogger("ffcx")


def generator(ir, parameters):
    """Generate UFC code for a form."""
    logger.info("Generating code for form:")
    logger.info(f"--- rank: {ir.rank}")
    logger.info(f"--- name: {ir.name}")

    import ffcx.codegeneration.C.cnodes as L

    d = {}
    d["factory_name"] = ir.name
    d["name_from_uflfile"] = ir.name_from_uflfile
    d["signature"] = f"\"{ir.signature}\""
    d["rank"] = ir.rank
    d["num_coefficients"] = ir.num_coefficients
    d["num_constants"] = ir.num_constants

    code = []
    cases = []
    for itg_type in ("cell", "interior_facet", "exterior_facet"):
        cases += [(L.Symbol(itg_type), L.Return(len(ir.subdomain_ids[itg_type])))]
    code += [L.Switch("integral_type", cases, default=L.Return(0))]
    d["num_integrals"] = L.StatementList(code)

    if len(ir.original_coefficient_position) > 0:
        d["original_coefficient_position_init"] = L.ArrayDecl(
            "int", f"original_coefficient_position_{ir.name}",
            values=ir.original_coefficient_position, sizes=len(ir.original_coefficient_position))
        d["original_coefficient_position"] = f"original_coefficient_position_{ir.name}"
    else:
        d["original_coefficient_position_init"] = ""
        d["original_coefficient_position"] = L.Null()

    cnames = ir.coefficient_names
    assert ir.num_coefficients == len(cnames)
    names = L.Symbol("names")
    if (len(cnames) == 0):
        code = [L.Return(L.Null())]
    else:
        code = [L.ArrayDecl("static const char*", names, len(cnames), cnames)]
        code += [L.Return(names)]
    d["coefficient_name_map"] = L.StatementList(code)

    cstnames = ir.constant_names
    names = L.Symbol("names")
    if len(cstnames) == 0:
        code = [L.Return(L.Null())]
    else:
        code = [L.ArrayDecl("static const char*", names, len(cstnames), cstnames)]
        code += [L.Return(names)]
    d["constant_name_map"] = L.StatementList(code)

    if len(ir.finite_elements) > 0:
        d["finite_elements"] = f"finite_elements_{ir.name}"
        d["finite_elements_init"] = L.ArrayDecl("ufcx_finite_element*", f"finite_elements_{ir.name}", values=[
                                                L.AddressOf(L.Symbol(el)) for el in ir.finite_elements],
                                                sizes=len(ir.finite_elements))
    else:
        d["finite_elements"] = L.Null()
        d["finite_elements_init"] = ""

    if len(ir.dofmaps) > 0:
        d["dofmaps"] = f"dofmaps_{ir.name}"
        d["dofmaps_init"] = L.ArrayDecl("ufcx_dofmap*", f"dofmaps_{ir.name}", values=[
            L.AddressOf(L.Symbol(dofmap)) for dofmap in ir.dofmaps], sizes=len(ir.dofmaps))
    else:
        d["dofmaps"] = L.Null()
        d["dofmaps_init"] = ""

    code = []
    cases = []
    code_ids = []
    cases_ids = []
    for itg_type in ("cell", "interior_facet", "exterior_facet"):
        if len(ir.integral_names[itg_type]) > 0:
            code += [L.ArrayDecl(
                "static ufcx_integral*", f"integrals_{itg_type}_{ir.name}",
                values=[L.AddressOf(L.Symbol(itg)) for itg in ir.integral_names[itg_type]],
                sizes=len(ir.integral_names[itg_type]))]
            cases.append((L.Symbol(itg_type), L.Return(L.Symbol(f"integrals_{itg_type}_{ir.name}"))))

            code_ids += [L.ArrayDecl(
                "static int", f"integral_ids_{itg_type}_{ir.name}",
                values=ir.subdomain_ids[itg_type], sizes=len(ir.subdomain_ids[itg_type]))]
            cases_ids.append((L.Symbol(itg_type), L.Return(L.Symbol(f"integral_ids_{itg_type}_{ir.name}"))))

    code += [L.Switch("integral_type", cases, default=L.Return(L.Null()))]
    code_ids += [L.Switch("integral_type", cases_ids, default=L.Return(L.Null()))]
    d["integrals"] = L.StatementList(code)

    d["integral_ids"] = L.StatementList(code_ids)

    code = []
    function_name = L.Symbol("function_name")

    # FIXME: Should be handled differently, revise how
    # ufcx_function_space is generated
    for (name, (element, dofmap, cmap_family, cmap_degree, cmap_celltype, cmap_variant)) in ir.function_spaces.items():
        code += [f"static ufcx_function_space functionspace_{name} ="]
        code += ["{"]
        code += [f".finite_element = &{element},"]
        code += [f".dofmap = &{dofmap},"]
        code += [f".geometry_family = \"{cmap_family}\","]
        code += [f".geometry_degree = {cmap_degree},"]
        code += [f".geometry_basix_cell = {int(cmap_celltype)},"]
        code += [f".geometry_basix_variant = {int(cmap_variant)}"]
        code += ["};"]

    _if = L.If
    for name in ir.function_spaces.keys():
        condition = L.EQ(L.Call("strcmp", (function_name, L.LiteralString(name))), 0)
        code += [_if(condition, L.Return(L.Symbol(f"&functionspace_{name}")))]
        _if = L.ElseIf

    code += ["return NULL;\n"]

    d["functionspace"] = L.StatementList(code)

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fields = [fname for _, fname, _, _ in Formatter().parse(form_template.factory) if fname]
    assert set(fields) == set(d.keys()), "Mismatch between keys in template and in formattting dict"

    # Format implementation code
    implementation = form_template.factory.format_map(d)

    # Format declaration
    declaration = form_template.declaration.format(factory_name=d["factory_name"],
                                                   name_from_uflfile=d["name_from_uflfile"])

    return declaration, implementation
