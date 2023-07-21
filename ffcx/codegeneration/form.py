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


def generator(ir, options):
    """Generate UFC code for a form."""
    logger.info("Generating code for form:")
    logger.info(f"--- rank: {ir.rank}")
    logger.info(f"--- name: {ir.name}")

    d = {}
    d["factory_name"] = ir.name
    d["name_from_uflfile"] = ir.name_from_uflfile
    d["signature"] = f'"{ir.signature}"'
    d["rank"] = ir.rank
    d["num_coefficients"] = ir.num_coefficients
    d["num_constants"] = ir.num_constants

    cases = "\n".join(
        f"case {itg_type}:\n return {len(ir.subdomain_ids[itg_type])};"
        for itg_type in ("cell", "interior_facet", "exterior_facet")
    )
    code = f"""switch (integral_type)
    {{
    {cases}
    }}"""
    d["num_integrals"] = code

    if len(ir.original_coefficient_position) > 0:
        n = len(ir.original_coefficient_position)
        vals = ", ".join(str(i) for i in ir.original_coefficient_position)
        d[
            "original_coefficient_position_init"
        ] = f"int original_coefficient_position_{ir.name}[{n}] = {{{vals}}};\n"

        d["original_coefficient_position"] = f"original_coefficient_position_{ir.name}"
    else:
        d["original_coefficient_position_init"] = ""
        d["original_coefficient_position"] = "NULL"

    cnames = ir.coefficient_names
    assert ir.num_coefficients == len(cnames)
    if len(cnames) == 0:
        code = "return NULL;"
    else:
        cnames_str = ", ".join('"{name}"' for name in cnames)
        code = f"static const char* names[{len(cnames)}]= {{{cnames_str}}};\n"
        code += "return names;\n"
    d["coefficient_name_map"] = code

    cstnames = ir.constant_names
    if len(cstnames) == 0:
        code = "return NULL;\n"
    else:
        cstnames_str = ", ".join('"{name}"' for name in cstnames)
        code = f"static const char* names[{len(cstnames)}] = {{{cstnames_str}}};\n"
        code += "return names;\n"
    d["constant_name_map"] = code

    if len(ir.finite_elements) > 0:
        d["finite_elements"] = f"finite_elements_{ir.name}"
        elems = ", ".join(f"&{el}" for el in ir.finite_elements)
        n = len(ir.finite_elements)
        d[
            "finite_elements_init"
        ] = f"ufcx_finite_element* finite_elements_{ir.name}[{n}] = {{{elems}}};"

    else:
        d["finite_elements"] = "NULL"
        d["finite_elements_init"] = ""

    if len(ir.dofmaps) > 0:
        d["dofmaps"] = f"dofmaps_{ir.name}"
        n = len(ir.dofmaps)
        vals = ", ".join(f"&{dm}" for dm in ir.dofmaps)
        d["dofmaps_init"] = f"ufcx_dofmap* dofmaps_{ir.name}[{n}] = {{{vals}}};"
    else:
        d["dofmaps"] = "NULL"
        d["dofmaps_init"] = ""

    code = []
    cases = []
    code_ids = []
    cases_ids = []
    for itg_type in ("cell", "interior_facet", "exterior_facet"):
        if len(ir.integral_names[itg_type]) > 0:
            vals = ", ".join(f"&{itg}" for itg in ir.integral_names[itg_type])
            n = len(ir.integral_names[itg_type])
            code += [
                f"static ufcx_integral* integrals_{itg_type}_{ir.name}[{n}] = {{{vals}}};"
            ]
            cases.append(f"case {itg_type}:")
            cases.append(f"return integrals_{itg_type}_{ir.name};")

            vals = ", ".join(str(i) for i in ir.subdomain_ids[itg_type])
            n = len(ir.subdomain_ids[itg_type])
            code_ids += [
                f"static int integral_ids_{itg_type}_{ir.name}[{n}] = {{{vals}}};"
            ]
            cases_ids.append(f"case {itg_type}:")
            cases_ids.append(f"return integral_ids_{itg_type}_{ir.name};")

    code += ["switch (integral_type){"]
    code += cases
    code += ["default:"]
    code += ["return NULL;}"]
    code_ids += ["switch(integral_type){"]
    code_ids += cases_ids
    code_ids += ["default:"]
    code_ids += ["return NULL;}"]

    d["integrals"] = "\n".join(code)
    d["integral_ids"] = "\n".join(code_ids)

    code = []
    function_name = "function_name"

    # FIXME: Should be handled differently, revise how
    # ufcx_function_space is generated
    for name, (
        element,
        dofmap,
        cmap_family,
        cmap_degree,
        cmap_celltype,
        cmap_variant,
    ) in ir.function_spaces.items():
        code += [f"static ufcx_function_space functionspace_{name} ="]
        code += ["{"]
        code += [f".finite_element = &{element},"]
        code += [f".dofmap = &{dofmap},"]
        code += [f'.geometry_family = "{cmap_family}",']
        code += [f".geometry_degree = {cmap_degree},"]
        code += [f".geometry_basix_cell = {int(cmap_celltype)},"]
        code += [f".geometry_basix_variant = {int(cmap_variant)}"]
        code += ["};"]

    for name in ir.function_spaces.keys():
        code += [
            f'if (strcmp({function_name}, "{name}") == 0) return(&functionspace_{name});'
        ]
    code += ["return NULL;\n"]

    d["functionspace"] = "\n".join(code)

    # Check that no keys are redundant or have been missed
    from string import Formatter

    fields = [
        fname for _, fname, _, _ in Formatter().parse(form_template.factory) if fname
    ]
    assert set(fields) == set(
        d.keys()
    ), "Mismatch between keys in template and in formatting dict"

    # Format implementation code
    implementation = form_template.factory.format_map(d)

    # Format declaration
    declaration = form_template.declaration.format(
        factory_name=d["factory_name"], name_from_uflfile=d["name_from_uflfile"]
    )

    return declaration, implementation
