# -*- coding: utf-8 -*-
"Code generator for quadrature representation."

# Copyright (C) 2009-2014 Kristian B. Oelgaard
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
# Modified by Mehdi Nikbakht 2010
# Modified by Anders Logg 2013-2014
# Modified by Martin Sandve Alnæs 2013-2014

# Python modules
import functools
import numpy

# UFL modules
from ufl.utils.sorting import sorted_by_key
from ufl.utils.derivativetuples import compute_derivative_tuples
from ufl import custom_integral_types

# FFC modules
from ffc.log import ffc_assert, error, warning
from ffc.quadrature.cpp import format, remove_unused, indent

from ffc.representationutils import initialize_integral_code

# Utility and optimization functions for quadraturegenerator
from ffc.quadrature.symbolics import generate_aux_constants


def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."

    # Prefetch formatting to speedup code generation (well...)
    ret = format["return"]

    # Generate code
    code = initialize_integral_code(ir, prefix, parameters)
    code["num_cells"] = indent(ret(ir["num_cells"]), 4)
    code["tabulate_tensor"] = indent(_tabulate_tensor(ir, prefix, parameters), 4)
    code["additional_includes_set"] = ir["additional_includes_set"]

    precision = ir["integrals_metadata"].get("precision")
    if precision is not None and precision != parameters["precision"]:
        warning("Ignoring precision in integral metadata compiled "
                "using quadrature representation. Not implemented.")

    return code


def _tabulate_tensor(ir, prefix, parameters):
    "Generate code for a single integral (tabulate_tensor())."

    # Prefetch formatting to speedup code generation
    f_comment = format["comment"]
    f_G = format["geometry constant"]
    f_const_double = format["assign"]
    f_switch = format["switch"]
    f_float = format["float"]
    f_assign = format["assign"]
    f_A = format["element tensor"]
    f_r = format["free indices"][0]
    f_loop = format["generate loop"]
    f_int = format["int"]
    f_facet = format["facet"]

    # Get data
    opt_par = ir["optimise_parameters"]
    integral_type = ir["integral_type"]
    gdim = ir["geometric_dimension"]
    tdim = ir["topological_dimension"]
    num_facets = ir["num_facets"]
    num_vertices = ir["num_vertices"]
    prim_idims = ir["prim_idims"]
    integrals = ir["trans_integrals"]
    geo_consts = ir["geo_consts"]
    oriented = ir["needs_oriented"]
    element_data = ir["element_data"]
    num_cells = ir["num_cells"]

    # Create sets of used variables
    used_weights = set()
    used_psi_tables = set()
    used_nzcs = set()
    trans_set = set()
    sets = [used_weights, used_psi_tables, used_nzcs, trans_set]

    affine_tables = {}  # TODO: This is not populated anywhere, remove?
    quadrature_rules = ir["quadrature_rules"]

    operations = []
    if integral_type == "cell":

        # Generate code for computing element tensor
        tensor_code, mem_code, num_ops = _generate_element_tensor(integrals,
                                                                  sets,
                                                                  opt_par,
                                                                  gdim,
                                                                  tdim)
        tensor_code = "\n".join(tensor_code)

        # Set operations equal to num_ops (for printing info on operations).
        operations.append([num_ops])

        # Generate code for basic geometric quantities
        jacobi_code = ""
        jacobi_code += format["compute_jacobian"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += format["compute_jacobian_inverse"](tdim, gdim)
        if oriented:
            jacobi_code += format["orientation"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += format["scale factor snippet"]

        # Generate code for cell volume and circumradius
        jacobi_code += "\n\n" + format["generate cell volume"](tdim, gdim,
                                                               integral_type)
        jacobi_code += "\n\n" + format["generate circumradius"](tdim, gdim,
                                                                integral_type)

    elif integral_type == "exterior_facet":

        # Iterate over facets
        cases = [None for i in range(num_facets)]
        for i in range(num_facets):
            # Update transformer with facets and generate case code +
            # set of used geometry terms.
            c, mem_code, ops = _generate_element_tensor(integrals[i], sets,
                                                        opt_par, gdim, tdim)
            case = [f_comment("Total number of operations to compute element tensor (from this point): %d" % ops)]
            case += c
            cases[i] = "\n".join(case)

            # Save number of operations (for printing info on
            # operations).
            operations.append([i, ops])

        # Generate tensor code for all cases using a switch.
        tensor_code = f_switch(f_facet(None), cases)

        # Generate code for basic geometric quantities
        jacobi_code = ""
        jacobi_code += format["compute_jacobian"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += format["compute_jacobian_inverse"](tdim, gdim)
        if oriented:
            jacobi_code += format["orientation"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += "\n\n" + format["facet determinant"](tdim, gdim)
        jacobi_code += "\n\n" + format["generate normal"](tdim, gdim,
                                                          integral_type)
        jacobi_code += "\n\n" + format["generate facet area"](tdim, gdim)
        if tdim == 3:
            jacobi_code += "\n\n" + format["generate min facet edge length"](tdim, gdim)
            jacobi_code += "\n\n" + format["generate max facet edge length"](tdim, gdim)

        # Generate code for cell volume and circumradius
        jacobi_code += "\n\n" + format["generate cell volume"](tdim, gdim,
                                                               integral_type)
        jacobi_code += "\n\n" + format["generate circumradius"](tdim, gdim,
                                                                integral_type)

    elif integral_type == "interior_facet":

        # Modify the dimensions of the primary indices because we have
        # a macro element
        prim_idims = [d * 2 for d in prim_idims]

        # Iterate over combinations of facets
        cases = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                # Update transformer with facets and generate case
                # code + set of used geometry terms.
                c, mem_code, ops = _generate_element_tensor(integrals[i][j],
                                                            sets, opt_par,
                                                            gdim, tdim)
                case = [f_comment("Total number of operations to compute element tensor (from this point): %d" % ops)]
                case += c
                cases[i][j] = "\n".join(case)

                # Save number of operations (for printing info on
                # operations).
                operations.append([i, j, ops])

        # Generate tensor code for all cases using a switch.
        tensor_code = f_switch(f_facet("+"),
                               [f_switch(f_facet("-"),
                                         cases[i]) for i in range(len(cases))])

        # Generate code for basic geometric quantities
        jacobi_code = ""
        for _r in ["+", "-"]:
            jacobi_code += format["compute_jacobian"](tdim, gdim, r=_r)
            jacobi_code += "\n"
            jacobi_code += format["compute_jacobian_inverse"](tdim, gdim, r=_r)
            if oriented:
                jacobi_code += format["orientation"](tdim, gdim, r=_r)
            jacobi_code += "\n"
        jacobi_code += "\n\n" + format["facet determinant"](tdim, gdim, r="+")
        jacobi_code += "\n\n" + format["generate normal"](tdim, gdim,
                                                          integral_type)
        jacobi_code += "\n\n" + format["generate facet area"](tdim, gdim)
        if tdim == 3:
            jacobi_code += "\n\n" + format["generate min facet edge length"](tdim, gdim, r="+")
            jacobi_code += "\n\n" + format["generate max facet edge length"](tdim, gdim, r="+")

        # Generate code for cell volume and circumradius
        jacobi_code += "\n\n" + format["generate cell volume"](tdim, gdim,
                                                               integral_type)
        jacobi_code += "\n\n" + format["generate circumradius"](tdim, gdim,
                                                                integral_type)

    elif integral_type == "vertex":

        # Iterate over vertices
        cases = [None for i in range(num_vertices)]
        for i in range(num_vertices):
            # Update transformer with vertices and generate case code
            # + set of used geometry terms.
            c, mem_code, ops = _generate_element_tensor(integrals[i],
                                                        sets,
                                                        opt_par,
                                                        gdim,
                                                        tdim)
            case = [f_comment("Total number of operations to compute element tensor (from this point): %d" % ops)]
            case += c
            cases[i] = "\n".join(case)

            # Save number of operations (for printing info on
            # operations).
            operations.append([i, ops])

        # Generate tensor code for all cases using a switch.
        tensor_code = f_switch(format["vertex"], cases)

        # Generate code for basic geometric quantities
        jacobi_code = ""
        jacobi_code += format["compute_jacobian"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += format["compute_jacobian_inverse"](tdim, gdim)
        if oriented:
            jacobi_code += format["orientation"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += "\n\n" + format["facet determinant"](tdim, gdim)  # FIXME: This is not defined in a point???

    elif integral_type in custom_integral_types:

        # Set number of cells
        if integral_type == "cutcell":
            num_cells = 1
        elif integral_type == "interface":
            num_cells = 2
        elif integral_type == "overlap":
            num_cells = 2

        # Warning that more than two cells in only partly supported.
        # The missing piece is to couple multiple cells to
        # restrictions other than '+' and '-'.
        if num_cells > 2:
            warning("Custom integrals with more than two cells only partly supported.")

        # Modify the dimensions of the primary indices because we have a macro element
        if num_cells == 2:
            prim_idims = [d * 2 for d in prim_idims]

        # Check whether we need to generate facet normals
        generate_custom_facet_normal = num_cells == 2

        # Generate code for computing element tensor
        tensor_code, mem_code, num_ops = _generate_element_tensor(integrals,
                                                                  sets,
                                                                  opt_par,
                                                                  gdim,
                                                                  tdim,
                                                                  generate_custom_facet_normal)

        tensor_code = "\n".join(tensor_code)

        # Set operations equal to num_ops (for printing info on
        # operations).
        operations.append([num_ops])

        # FIXME: Jacobi code is only needed when we use cell volume or
        # circumradius.
        # FIXME: Does not seem to be removed by removed_unused.

        # Generate code for basic geometric quantities
        jacobi_code = ""
        for i in range(num_cells):
            r = i if num_cells > 1 else None
            jacobi_code += "\n"
            jacobi_code += f_comment("--- Compute geometric quantities on cell %d ---" % i)
            jacobi_code += "\n\n"
            if num_cells > 1:
                jacobi_code += f_comment("Extract vertex coordinates\n")
                jacobi_code += format["extract_cell_coordinates"]((tdim + 1) * gdim * i, r=i)
                jacobi_code += "\n\n"
            jacobi_code += format["compute_jacobian"](tdim, gdim, r=r)
            jacobi_code += "\n"
            jacobi_code += format["compute_jacobian_inverse"](tdim, gdim, r=r)
            jacobi_code += "\n"
            jacobi_code += format["generate cell volume"](tdim, gdim,
                                                          integral_type,
                                                          r=r if num_cells > 1 else None)
            jacobi_code += "\n"
            jacobi_code += format["generate circumradius"](tdim, gdim,
                                                           integral_type,
                                                           r=r if num_cells > 1 else None)
            jacobi_code += "\n"

    else:
        error("Unhandled integral type: " + str(integral_type))

    # After we have generated the element code for all facets we can
    # remove the unused transformations.
    common = [remove_unused(jacobi_code, trans_set)]

    # FIXME: After introduction of custom integrals, the common code
    # here is not really common anymore. Think about how to
    # restructure this function.

    # Add common code except for custom integrals
    if integral_type not in custom_integral_types:
        common += _tabulate_weights([quadrature_rules[p] for p in sorted(used_weights)])

        # Add common code for updating tables
        name_map = ir["name_map"]
        tables = ir["unique_tables"]
        tables.update(affine_tables)  # TODO: This is not populated anywhere, remove?
        common += _tabulate_psis(tables, used_psi_tables, name_map, used_nzcs,
                                 opt_par, integral_type, gdim)

    # Add special tabulation code for custom integral
    else:
        common += _evaluate_basis_at_quadrature_points(used_psi_tables,
                                                       gdim,
                                                       element_data,
                                                       prefix,
                                                       num_vertices,
                                                       num_cells)

    # Reset the element tensor (array 'A' given as argument to
    # tabulate_tensor() by assembler)
    # Handle functionals.
    common += [f_comment("Reset values in the element tensor.")]
    if prim_idims == []:
        common += [f_assign(f_A(f_int(0)), f_float(0))]
    else:
        dim = functools.reduce(lambda v, u: v * u, prim_idims)
        common += f_loop([f_assign(f_A(f_r), f_float(0))], [(f_r, 0, dim)])

    # Create the constant geometry declarations (only generated if
    # simplify expressions are enabled).
    geo_ops, geo_code = generate_aux_constants(geo_consts, f_G, f_const_double)
    if geo_code:
        common += [f_comment("Number of operations to compute geometry constants: %d." % geo_ops)]
        common += [format["declaration"](format["float declaration"], f_G(len(geo_consts)))]
        common += geo_code

    # Add comments.
    common += ["", f_comment("Compute element tensor using UFL quadrature representation")]
    common += [f_comment("Optimisations: %s" % ", ".join([str((k, opt_par[k]))
                                                          for k in sorted(opt_par.keys())]))]

    for ops in operations:
        # Add geo ops count to integral ops count for writing info.
        if isinstance(ops[-1], int):
            ops[-1] += geo_ops

    return "\n".join(common) + "\n" + tensor_code


def _generate_element_tensor(integrals, sets, optimise_parameters, gdim, tdim,
                             generate_custom_facet_normal=False):
    "Construct quadrature code for element tensors."

    # Prefetch formats to speed up code generation.
    f_comment = format["comment"]
    f_ip = format["integration points"]
    f_I = format["ip constant"]
    f_loop = format["generate loop"]
    f_ip_coords = format["generate ip coordinates"]
    f_coords = format["coordinate_dofs"]
    f_double = format["float declaration"]
    f_decl = format["declaration"]
    f_X = format["ip coordinates"]
    f_C = format["conditional"]

    # Initialise return values.
    element_code = []
    tensor_ops_count = 0

    # TODO: KBO: The members_code was used when I generated the
    # load_table.h file which could load tables of basisfunction. This
    # feature has not been reimplemented. However, with the new design
    # where we only tabulate unique tables (and only non-zero entries)
    # it doesn't seem to be necessary. Should it be deleted?
    members_code = ""

    # We receive a dictionary {num_points: form,}.
    # Loop points and forms.
    for points, terms, functions, ip_consts, coordinate, conditionals in integrals:

        element_code += ["", f_comment("Loop quadrature points for integral.")]

        ip_code = []
        num_ops = 0

        # Generate code to compute coordinates if used.
        if coordinate:
            name, gdim, ip, r = coordinate
            element_code += ["", f_comment("Declare array to hold physical coordinate of quadrature point.")]
            element_code += [f_decl(f_double, f_X(points, gdim))]
            ops, coord_code = f_ip_coords(gdim, tdim, points, name, ip, r)
            ip_code += ["", f_comment("Compute physical coordinate of quadrature point, operations: %d." % ops)]
            ip_code += [coord_code]
            num_ops += ops
            # Update used psi tables and transformation set.
            sets[1].add(name)
            sets[3].add(f_coords(r))

        # Generate code to compute function values.
        if functions:
            func_code, ops = _generate_functions(functions, sets)
            ip_code += func_code
            num_ops += ops

        # Generate code to compute conditionals (might depend on coordinates
        # and function values so put here).
        # TODO: Some conditionals might only depend on geometry so they
        # should be moved outside if possible.
        if conditionals:
            ip_code += [f_decl(f_double, f_C(len(conditionals)))]
            # Sort conditionals (need to in case of nested conditionals).
            reversed_conds = dict([(n, (o, e)) for e, (t, o, n) in conditionals.items()])
            for num in range(len(conditionals)):
                name = format["conditional"](num)
                ops, expr = reversed_conds[num]
                ip_code += [f_comment("Compute conditional, operations: %d." % ops)]
                ip_code += [format["assign"](name, expr)]
                num_ops += ops

        # Generate code for ip constant declarations.
        ip_const_ops, ip_const_code = generate_aux_constants(ip_consts, f_I,
                                                             format["assign"], True)
        num_ops += ip_const_ops
        if ip_const_code:
            ip_code += ["", f_comment("Number of operations to compute ip constants: %d" % ip_const_ops)]
            ip_code += [format["declaration"](format["float declaration"], f_I(len(ip_consts)))]
            ip_code += ip_const_code

        # Generate code to evaluate the element tensor.
        integral_code, ops = _generate_integral_code(points, terms, sets,
                                                     optimise_parameters)
        num_ops += ops
        if points is None:
            quadrature_ops = "unknown"
            tensor_ops_count = "unknown"
        else:
            quadrature_ops = num_ops * points
            tensor_ops_count += quadrature_ops
        ip_code += integral_code
        element_code.append(f_comment
                            ("Number of operations to compute element tensor for following IP loop = %s" % str(quadrature_ops)))

        # Generate code for custom facet normal if necessary
        if generate_custom_facet_normal:
            for line in ip_code:
                if "n_00" in line:
                    ip_code = [format["facet_normal_custom"](gdim)] + ip_code
                    break

        # Loop code over all IPs.
        if points == 0:
            element_code.append(f_comment("Only 1 integration point, omitting IP loop."))
            element_code += ip_code
        elif points is None:
            num_points = "num_quadrature_points"
            element_code += f_loop(ip_code, [(f_ip, 0, num_points)])
        else:
            element_code += f_loop(ip_code, [(f_ip, 0, points)])

    return (element_code, members_code, tensor_ops_count)


def _generate_functions(functions, sets):
    "Generate declarations for functions and code to compute values."

    f_comment = format["comment"]
    f_double = format["float declaration"]
    f_F = format["function value"]
    f_float = format["floating point"]
    f_decl = format["declaration"]
    f_r = format["free indices"][0]
    f_iadd = format["iadd"]
    f_loop = format["generate loop"]

    # Create the function declarations.
    code = ["", f_comment("Coefficient declarations.")]
    code += [f_decl(f_double, f_F(n), f_float(0)) for n in range(len(functions))]

    # Get sets.
    used_psi_tables = sets[1]
    used_nzcs = sets[2]

    # Sort functions after loop ranges.
    function_list = {}
    for key, val in functions.items():
        if val[1] in function_list:
            function_list[val[1]].append(key)
        else:
            function_list[val[1]] = [key]

    total_ops = 0
    # Loop ranges and get list of functions.
    for loop_range, list_of_functions in sorted(function_list.items()):
        function_expr = {}
        # Loop functions.
        func_ops = 0
        for function in list_of_functions:
            # Get name and number.
            number, range_i, ops, psi_name, u_nzcs, ufl_element = functions[function]

            # Add name to used psi names and non zeros name to used_nzcs.
            used_psi_tables.add(psi_name)
            used_nzcs.update(u_nzcs)

            # TODO: This check can be removed for speed later.
            ffc_assert(number not in function_expr, "This is definitely not supposed to happen!")

            # Convert function (that might be a symbol) to a simple
            # string and save.
            function = str(function)
            function_expr[number] = function

            # Get number of operations to compute entry and add to
            # function operations count.
            func_ops += (ops + 1) * range_i

        # Add function operations to total count
        total_ops += func_ops
        code += ["", f_comment("Total number of operations to compute function values = %d" % func_ops)]

        # Sort the functions according to name and create loop to
        # compute the function values.
        lines = [f_iadd(f_F(n), function_expr[n]) for n in sorted(function_expr.keys())]
        code += f_loop(lines, [(f_r, 0, loop_range)])  # TODO: If loop_range == 1, this loop may be unneccessary. Not sure if it's safe to just skip it.

    return code, total_ops


def _generate_integral_code(points, terms, sets, optimise_parameters):
    "Generate code to evaluate the element tensor."

    # Prefetch formats to speed up code generation.
    f_comment = format["comment"]
    f_iadd = format["iadd"]
    f_A = format["element tensor"]
    f_loop = format["generate loop"]
    f_B = format["basis constant"]

    # Initialise return values.
    code = []
    num_ops = 0
    loops = {}

    # Extract sets.
    used_weights, used_psi_tables, used_nzcs, trans_set = sets

    # Loop terms and create code.
    for loop, (data, entry_vals) in sorted(terms.items()):

        # If we don't have any entry values, there's no need to
        # generate the loop.
        if not entry_vals:
            continue

        # Get data.
        t_set, u_weights, u_psi_tables, u_nzcs, basis_consts = data

        # If we have a value, then we also need to update the sets of
        # used variables.
        trans_set.update(t_set)
        used_weights.update(u_weights)
        used_psi_tables.update(u_psi_tables)
        used_nzcs.update(u_nzcs)

        # Generate code for basis constant declarations.
        basis_const_ops, basis_const_code = generate_aux_constants(basis_consts,
                                                                   f_B,
                                                                   format["assign"], True)
        decl_code = []
        if basis_consts:
            decl_code = [format["declaration"](format["float declaration"],
                                               f_B(len(basis_consts)))]
        loops[loop] = [basis_const_ops, decl_code + basis_const_code]

        for entry, value, ops in entry_vals:
            # Compute number of operations to compute entry (add 1
            # because of += in assignment).
            entry_ops = ops + 1

            # Create comment for number of operations
            entry_ops_comment = f_comment("Number of operations to compute entry: %d" % entry_ops)
            entry_code = f_iadd(f_A(entry), value)
            loops[loop][0] += entry_ops
            loops[loop][1] += [entry_ops_comment, entry_code]

    # Write all the loops of basis functions.
    for loop, ops_lines in sorted(loops.items()):
        ops, lines = ops_lines
        prim_ops = functools.reduce(lambda i, j: i * j, [ops] + [l[2] for l in loop])
        # Add number of operations for current loop to total count.
        num_ops += prim_ops
        code += ["", f_comment("Number of operations for primary indices: %d" % prim_ops)]
        code += f_loop(lines, loop)

    return code, num_ops


def _tabulate_weights(quadrature_rules):
    "Generate table of quadrature weights."

    # Prefetch formats to speed up code generation.
    f_float = format["floating point"]
    f_table = format["static const float declaration"]
    f_sep = format["list separator"]
    f_weight = format["weight"]
    f_component = format["component"]
    f_group = format["grouping"]
    f_decl = format["declaration"]
    f_tensor = format["tabulate tensor"]
    f_comment = format["comment"]
    f_int = format["int"]

    code = ["", f_comment("Array of quadrature weights.")]

    # Loop tables of weights and create code.
    for weights, points in quadrature_rules:
        # FIXME: For now, raise error if we don't have weights.
        # We might want to change this later.
        ffc_assert(weights.any(), "No weights.")

        # Create name and value for weight.
        num_points = len(points)
        name = f_weight(num_points)
        value = f_float(weights[0])
        if len(weights) > 1:
            name += f_component("", f_int(num_points))
            value = f_tensor(weights)
        code += [f_decl(f_table, name, value)]

        # Tabulate the quadrature points (uncomment for different
        # parameters).  1) Tabulate the points as: p0, p1, p2, with p0
        # = (x0, y0, z0) etc.  Use f_float to format the value (enable
        # variable precision).
        formatted_points = [f_group(f_sep.join([f_float(val) for val in point]))
                            for point in points]

        # Create comment.
        comment = "Quadrature points on the UFC reference element: " \
                  + f_sep.join(formatted_points)
        code += [f_comment(comment)]

        # 2) Tabulate the coordinates of the points p0, p1, p2 etc.
        #    X: x0, x1, x2
        #    Y: y0, y1, y2
        #    Z: z0, z1, z2
#            comment = "Quadrature coordinates on the UFC reference element: "
#            code += [format["comment"](comment)]

# All points have the same number of coordinates.
#            num_coord = len(points[0])

# All points have x-coordinates.
#            xs = [f_float(p[0]) for p in points]
#            comment = "X: " + f_sep.join(xs)
#            code += [format["comment"](comment)]

#            ys = []
#            zs = []
# Tabulate y-coordinate if we have 2 or more coordinates.
#            if num_coord >= 2:
#                ys = [f_float(p[1]) for p in points]
#                comment = "Y: " + f_sep.join(ys)
#                code += [format["comment"](comment)]
# Only tabulate z-coordinate if we have 3 coordinates.
#            if num_coord == 3:
#                zs = [f_float(p[2]) for p in points]
#                comment = "Z: " + f_sep.join(zs)
#                code += [format["comment"](comment)]

        code += [""]

    return code


def _tabulate_psis(tables, used_psi_tables, inv_name_map, used_nzcs,
                   optimise_parameters, integral_type, gdim):
    "Tabulate values of basis functions and their derivatives at quadrature points."

    # Prefetch formats to speed up code generation.
    f_comment = format["comment"]
    f_table = format["static const float declaration"]
    f_component = format["component"]
    f_const_uint = format["static const uint declaration"]
    f_nzcolumns = format["nonzero columns"]
    f_list = format["list"]
    f_decl = format["declaration"]
    f_tensor = format["tabulate tensor"]
    f_new_line = format["new line"]
    f_int = format["int"]

    # FIXME: Check if we can simplify the tabulation
    code = []
    code += [f_comment("Values of basis functions at quadrature points.")]

    # Get list of non zero columns, if we ignore ones, ignore columns
    # with one component.
    if optimise_parameters["ignore ones"]:
        nzcs = []
        for key, val in sorted(inv_name_map.items()):
            # Check if we have a table of ones or if number of
            # non-zero columns is larger than one.
            if val[1] and len(val[1][1]) > 1 or not val[3]:
                nzcs.append(val[1])
    else:
        nzcs = [val[1] for key, val in sorted(inv_name_map.items())
                if val[1]]

    # TODO: Do we get arrays that are not unique?
    new_nzcs = []
    for nz in nzcs:
        # Only get unique arrays.
        if nz not in new_nzcs:
            new_nzcs.append(nz)

    # Construct name map.
    name_map = {}
    if inv_name_map:
        for name in sorted(inv_name_map):
            if inv_name_map[name][0] in name_map:
                name_map[inv_name_map[name][0]].append(name)
            else:
                name_map[inv_name_map[name][0]] = [name]

    # Loop items in table and tabulate.
    for name in sorted(used_psi_tables):

        # Only proceed if values are still used (if they're not
        # remapped).
        vals = tables[name]

        if vals is not None:

            # Add declaration to name.
            ip, dofs = numpy.shape(vals)
            decl_name = f_component(name, [ip, dofs])

            # Generate array of values.
            value = f_tensor(vals)
            code += [f_decl(f_table, decl_name, f_new_line + value), ""]

        # Tabulate non-zero indices.
        if optimise_parameters["eliminate zeros"]:
            if name in sorted(name_map):
                for n in name_map[name]:
                    if inv_name_map[n][1] and inv_name_map[n][1] in new_nzcs:
                        i, cols = inv_name_map[n][1]
                        if i not in used_nzcs:
                            continue
                        code += [f_comment("Array of non-zero columns")]
                        value = f_list([f_int(c) for c in list(cols)])
                        name_col = f_component(f_nzcolumns(i), len(cols))
                        code += [f_decl(f_const_uint, name_col, value), ""]

                        # Remove from list of columns.
                        new_nzcs.remove(inv_name_map[n][1])
    return code


def _evaluate_basis_at_quadrature_points(psi_tables, gdim, element_data,
                                         form_prefix, num_vertices, num_cells):
    "Generate code for calling evaluate basis (derivatives) at quadrature points"

    # Prefetch formats to speed up code generation
    f_comment = format["comment"]
    f_declaration = format["declaration"]
    f_static_array = format["static array"]
    f_loop = format["generate loop"]
    f_eval_basis_decl = format["eval_basis_decl"]
    f_eval_basis_init = format["eval_basis_init"]
    f_eval_basis = format["eval_basis"]
    f_eval_basis_copy = format["eval_basis_copy"]
    f_eval_derivs_decl = format["eval_derivs_decl"]
    f_eval_derivs_init = format["eval_derivs_init"]
    f_eval_derivs = format["eval_derivs"]
    f_eval_derivs_copy = format["eval_derivs_copy"]

    code = []

    # Extract prefixes for tables
    prefixes = sorted(set(table.split("_")[0] for table in psi_tables))

    # Use lower case prefix for form name
    form_prefix = form_prefix.lower()

    # The psi_tables used by the quadrature code are for scalar
    # components of specific derivatives, while tabulate_basis_all and
    # tabulate_basis_derivatives_all return data including all
    # possible components and derivatives. We therefore need to
    # iterate over prefixes (= elements), call tabulate_basis_all or
    # tabulate_basis_derivatives all, and then extract the relevant
    # data and fill in the psi_tables. We therefore need to extract
    # for each prefix, which tables need to be filled in.

    # For each unique prefix, check which derivatives and components
    # are used
    used_derivatives_and_components = {}
    for prefix in prefixes:
        used_derivatives_and_components[prefix] = {}
        for table in psi_tables:
            if prefix not in table:
                continue

            # Check for derivative
            if "_D" in table:
                d = table.split("_D")[1].split("_")[0]
                n = sum([int(_d) for _d in d])  # FIXME: Assume at most 9 derivatives...
            else:
                n = 0

            # Check for component
            if "_C" in table:
                c = table.split("_C")[1].split("_")[0]
            else:
                c = None

            # Note that derivative has been used
            if n not in used_derivatives_and_components[prefix]:
                used_derivatives_and_components[prefix][n] = set()
            used_derivatives_and_components[prefix][n].add(c)

    # Generate code for setting quadrature weights
    code += [f_comment("Set quadrature weights")]
    code += [f_declaration("const double*", "W", "quadrature_weights")]
    code += [""]

    # Generate code for calling evaluate_basis_[derivatives_]all
    for prefix in prefixes:

        # Get element data for current element
        counter = int(prefix.split("FE")[1])
        space_dim = element_data[counter]["num_element_dofs"]
        value_size = element_data[counter]["physical_value_size"]
        element_classname = element_data[counter]["classname"]

        # Iterate over derivative orders
        for n, components in sorted_by_key(used_derivatives_and_components[prefix]):
            # components are a set and need to be sorted
            components = sorted(components)

            # Code for evaluate_basis_all (n = 0 which means it's not
            # a derivative)
            if n == 0:

                code += [f_comment("--- Evaluation of basis functions ---")]
                code += [""]

                # Compute variables for code generation
                eval_stride = value_size
                eval_size = space_dim * eval_stride
                table_size = num_cells * space_dim

                # Iterate over components and initialize tables
                for c in components:

                    # Set name of table
                    if c is None:
                        table_name = prefix
                    else:
                        table_name = prefix + "_C%s" % c

                    # Generate code for declaration of table
                    code += [f_comment("Create table %s for basis function values on all cells" % table_name)]
                    code += [f_eval_basis_decl % {"table_name": table_name}]
                    code += [f_eval_basis_init % {"table_name": table_name,
                                                  "table_size": table_size}]
                    code += [""]

                # Iterate over cells in macro element and evaluate basis
                for cell_number in range(num_cells):

                    # Compute variables for code generation
                    eval_name = "%s_values_%d" % (prefix, cell_number)
                    table_offset = cell_number * space_dim
                    vertex_offset = cell_number * num_vertices * gdim

                    # Generate block of code for loop
                    block = []

                    # Generate code for calling evaluate_basis_all
                    block += [f_eval_basis % {"classname": element_classname,
                                              "eval_name": eval_name,
                                              "gdim": gdim,
                                              "vertex_offset": vertex_offset}]

                    # Iterate over components and extract values
                    for c in components:

                        # Set name of table and component offset
                        if c is None:
                            table_name = prefix
                            eval_offset = 0
                        else:
                            table_name = prefix + "_C%s" % c
                            eval_offset = int(c)

                        # Generate code for copying values
                        block += [""]
                        block += [f_eval_basis_copy % {"table_name": table_name,
                                                       "eval_name": eval_name,
                                                       "eval_stride": eval_stride,
                                                       "eval_offset": eval_offset,
                                                       "space_dim": space_dim,
                                                       "table_offset": table_offset}]

                    # Generate code
                    code += [f_comment("Evaluate basis functions on cell %d" % cell_number)]
                    code += [f_static_array("double", eval_name, eval_size)]
                    code += f_loop(block, [("ip", 0, "num_quadrature_points")])
                    code += [""]

            # Code for evaluate_basis_derivatives_all (derivative of degree n > 0)
            else:

                code += [f_comment("--- Evaluation of basis function derivatives of order %d ---" % n)]
                code += [""]

                # FIXME: We extract values for all possible
                # derivatives, even
                # FIXME: if not all are used. (For components, we
                # extract only
                # FIXME: components that are actually used.) This may
                # be optimized
                # FIXME: but the extra cost is likely small.

                # Get derivative tuples
                __, deriv_tuples = compute_derivative_tuples(n, gdim)

                # Generate names for derivatives
                derivs = ["".join(str(_d) for _d in d) for d in deriv_tuples]

                # Compute variables for code generation
                eval_stride = value_size * len(derivs)
                eval_size = space_dim * eval_stride
                table_size = num_cells * space_dim

                # Iterate over derivatives and initialize tables
                seen_derivs = set()
                for d in derivs:

                    # Skip derivative if seen before (d^2/dxdy = d^2/dydx)
                    if d in seen_derivs:
                        continue
                    seen_derivs.add(d)

                    # Iterate over components
                    for c in components:

                        # Set name of table
                        if c is None:
                            table_name = prefix + "_D%s" % d
                        else:
                            table_name = prefix + "_C%s_D%s" % (c, d)

                        # Generate code for declaration of table
                        code += [f_comment("Create table %s for basis function derivatives on all cells" % table_name)]
                        code += [(f_eval_derivs_decl % {"table_name": table_name})]
                        code += [(f_eval_derivs_init % {"table_name": table_name,
                                                        "table_size": table_size})]
                        code += [""]

                # Iterate over cells (in macro element)
                for cell_number in range(num_cells):

                    # Compute variables for code generation
                    eval_name = "%s_dvalues_%d_%d" % (prefix, n, cell_number)
                    table_offset = cell_number * space_dim
                    vertex_offset = cell_number * num_vertices * gdim

                    # Generate block of code for loop
                    block = []

                    # Generate code for calling evaluate_basis_derivatives_all
                    block += [f_eval_derivs % {"classname": element_classname,
                                               "eval_name": eval_name,
                                               "gdim": gdim,
                                               "vertex_offset": vertex_offset,
                                               "n": n}]

                    # Iterate over derivatives and extract values
                    seen_derivs = set()
                    for i, d in enumerate(derivs):

                        # Skip derivative if seen before (d^2/dxdy = d^2/dydx)
                        if d in seen_derivs:
                            continue
                        seen_derivs.add(d)

                        # Iterate over components
                        for c in components:

                            # Set name of table and component offset
                            if c is None:
                                table_name = prefix + "_D%s" % d
                                eval_offset = i
                            else:
                                table_name = prefix + "_C%s_D%s" % (c, d)
                                eval_offset = len(derivs) * int(c) + i

                            # Generate code for copying values
                            block += [""]
                            block += [(f_eval_derivs_copy % {"table_name": table_name,
                                                             "eval_name": eval_name,
                                                             "eval_stride": eval_stride,
                                                             "eval_offset": eval_offset,
                                                             "space_dim": space_dim,
                                                             "table_offset": table_offset})]

                    # Generate code
                    code += [f_comment("Evaluate basis function derivatives on cell %d" % cell_number)]
                    code += [f_static_array("double", eval_name, eval_size)]
                    code += f_loop(block, [("ip", 0, "num_quadrature_points")])
                    code += [""]

                # Add newline
                code += [""]

    return code
