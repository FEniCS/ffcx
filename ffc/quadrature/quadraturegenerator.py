"Code generator for quadrature representation."

# Copyright (C) 2009-2013 Kristian B. Oelgaard
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
# Modified by Anders Logg 2013
# Modified by Martin Alnaes, 2013
#
# First added:  2009-01-07
# Last changed: 2013-02-20

# Python modules.
import functools
import numpy

# UFL modules.
from ufl.algorithms.printing import tree_format

## FFC modules.
from ffc.log import info, debug, ffc_assert
from ffc.cpp import format, remove_unused

from ffc.representationutils import initialize_integral_code

# Utility and optimisation functions for quadraturegenerator.
from symbolics import generate_aux_constants

def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."
    code = initialize_integral_code(ir, prefix, parameters)
    code["tabulate_tensor"] = _tabulate_tensor(ir, parameters)
    code["additional_includes_set"] = ir["additional_includes_set"]
    return code

def _tabulate_tensor(ir, parameters):
    "Generate code for a single integral (tabulate_tensor())."

    f_comment       = format["comment"]
    f_G             = format["geometry constant"]
    f_const_double  = format["assign"]
    f_switch        = format["switch"]
    f_float         = format["float"]
    f_assign        = format["assign"]
    f_A             = format["element tensor"]
    f_r             = format["free indices"][0]
    f_loop          = format["generate loop"]
    f_int           = format["int"]
    f_facet         = format["facet"]

    # Get data.
    opt_par     = ir["optimise_parameters"]
    domain_type = ir["domain_type"]
    gdim        = ir["geometric_dimension"]
    tdim        = ir["topological_dimension"]
    num_facets  = ir["num_facets"]
    num_vertices= ir["num_vertices"]
    prim_idims  = ir["prim_idims"]
    integrals   = ir["trans_integrals"]
    geo_consts  = ir["geo_consts"]
    oriented    = ir["needs_oriented"]

    # Create sets of used variables.
    used_weights    = set()
    used_psi_tables = set()
    used_nzcs       = set()
    trans_set       = set()
    sets = [used_weights, used_psi_tables, used_nzcs, trans_set]

    affine_tables = {} # TODO: This is not populated anywhere, remove?
    quadrature_weights = ir["quadrature_weights"]

    operations = []
    if domain_type == "cell":
        # Update treansformer with facets and generate code + set of used geometry terms.
        tensor_code, mem_code, num_ops = _generate_element_tensor(integrals, sets, \
                                         opt_par)
        tensor_code = "\n".join(tensor_code)

        # Set operations equal to num_ops (for printing info on operations).
        operations.append([num_ops])

        # Generate code for basic geometric quantities
        jacobi_code  = ""
        jacobi_code += format["compute_jacobian"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += format["compute_jacobian_inverse"](tdim, gdim)
        if oriented:
            jacobi_code += format["orientation"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += format["scale factor snippet"]

    elif domain_type == "exterior_facet":
        cases = [None for i in range(num_facets)]
        for i in range(num_facets):
            # Update treansformer with facets and generate case code + set of used geometry terms.
            c, mem_code, ops = _generate_element_tensor(integrals[i], sets, opt_par)
            case = [f_comment("Total number of operations to compute element tensor (from this point): %d" % ops)]
            case += c
            cases[i] = "\n".join(case)

            # Save number of operations (for printing info on operations).
            operations.append([i, ops])

        # Generate tensor code for all cases using a switch.
        tensor_code = f_switch(f_facet(None), cases)

        # Generate code for basic geometric quantities
        jacobi_code  = ""
        jacobi_code += format["compute_jacobian"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += format["compute_jacobian_inverse"](tdim, gdim)
        if oriented:
            jacobi_code += format["orientation"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += "\n\n" + format["facet determinant"](tdim, gdim)
        jacobi_code += "\n\n" + format["generate normal"](tdim, gdim, domain_type)
        jacobi_code += "\n\n" + format["generate facet area"](tdim, gdim)
        if tdim == 3:
            jacobi_code += "\n\n" + format["generate min facet edge length"](tdim, gdim)
            jacobi_code += "\n\n" + format["generate max facet edge length"](tdim, gdim)

    elif domain_type == "interior_facet":
        # Modify the dimensions of the primary indices because we have a macro element
        prim_idims = [d*2 for d in prim_idims]

        cases = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                # Update transformer with facets and generate case code + set of used geometry terms.
                c, mem_code, ops = _generate_element_tensor(integrals[i][j], sets, \
                                                            opt_par)
                case = [f_comment("Total number of operations to compute element tensor (from this point): %d" % ops)]
                case += c
                cases[i][j] = "\n".join(case)

                # Save number of operations (for printing info on operations).
                operations.append([i, j, ops])

        # Generate tensor code for all cases using a switch.
        tensor_code = f_switch(f_facet("+"), [f_switch(f_facet("-"), cases[i]) for i in range(len(cases))])

        # Generate code for basic geometric quantities
        jacobi_code  = ""
        for _r in ["+", "-"]:
            jacobi_code += format["compute_jacobian"](tdim, gdim, r=_r)
            jacobi_code += "\n"
            jacobi_code += format["compute_jacobian_inverse"](tdim, gdim, r=_r)
            if oriented:
                jacobi_code += format["orientation"](tdim, gdim)
            jacobi_code += "\n"
        jacobi_code += "\n\n" + format["facet determinant"](tdim, gdim, r="+")
        jacobi_code += "\n\n" + format["generate normal"](tdim, gdim, domain_type)
        jacobi_code += "\n\n" + format["generate facet area"](tdim, gdim)
        if tdim == 3:
            jacobi_code += "\n\n" + format["generate min facet edge length"](tdim, gdim, r="+")
            jacobi_code += "\n\n" + format["generate max facet edge length"](tdim, gdim, r="+")

    elif domain_type == "point":
        cases = [None for i in range(num_vertices)]
        for i in range(num_vertices):
            # Update treansformer with vertices and generate case code +
            # set of used geometry terms.
            c, mem_code, ops = _generate_element_tensor(integrals[i],
                                                        sets, opt_par)
            case = [f_comment("Total number of operations to compute element tensor (from this point): %d" % ops)]
            case += c
            cases[i] = "\n".join(case)

            # Save number of operations (for printing info on operations).
            operations.append([i, ops])

        # Generate tensor code for all cases using a switch.
        tensor_code = f_switch(format["vertex"], cases)

        # Generate code for basic geometric quantities
        jacobi_code  = ""
        jacobi_code += format["compute_jacobian"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += format["compute_jacobian_inverse"](tdim, gdim)
        if oriented:
            jacobi_code += format["orientation"](tdim, gdim)
        jacobi_code += "\n"
        jacobi_code += "\n\n" + format["facet determinant"](tdim, gdim) # FIXME: This is not defined in a point???

    else:
        error("Unhandled integral type: " + str(integral_type))

    # Add common (for cell, exterior and interior) geo code.
    if domain_type != "point":
        jacobi_code += "\n\n" + format["generate cell volume"](tdim, gdim, domain_type)
        jacobi_code += "\n\n" + format["generate circumradius"](tdim, gdim, domain_type)

    # After we have generated the element code for all facets we can remove
    # the unused transformations and tabulate the used psi tables and weights.
    common = [remove_unused(jacobi_code, trans_set)]
    common += _tabulate_weights([quadrature_weights[p] for p in used_weights])
    name_map = ir["name_map"]
    tables = ir["unique_tables"]
    tables.update(affine_tables) # TODO: This is not populated anywhere, remove?
    common += _tabulate_psis(tables, used_psi_tables, name_map, used_nzcs, opt_par)

    # Reset the element tensor (array 'A' given as argument to tabulate_tensor() by assembler)
    # Handle functionals.
    common += [f_comment("Reset values in the element tensor.")]
    value = f_float(0)
    if prim_idims == []:
        common += [f_assign(f_A(f_int(0)), f_float(0))]
    else:
        dim = functools.reduce(lambda v,u: v*u, prim_idims)
        common += f_loop([f_assign(f_A(f_r), f_float(0))], [(f_r, 0, dim)])

    # Create the constant geometry declarations (only generated if simplify expressions are enabled).
    geo_ops, geo_code = generate_aux_constants(geo_consts, f_G, f_const_double)
    if geo_code:
        common += [f_comment("Number of operations to compute geometry constants: %d." % geo_ops)]
        common += [format["declaration"](format["float declaration"], f_G(len(geo_consts)))]
        common += geo_code

    # Add comments.
    common += ["", f_comment("Compute element tensor using UFL quadrature representation")]
    common += [f_comment("Optimisations: %s" % ", ".join([str((k, opt_par[k]))\
                for k in sorted(opt_par.keys())]))]

    # Print info on operation count.
    message = {"cell":           "Cell, number of operations to compute tensor: %d",
               "exterior_facet": "Exterior facet %d, number of operations to compute tensor: %d",
               "interior_facet": "Interior facets (%d, %d), number of operations to compute tensor: %d",
               "point": "Point %d, number of operations to compute tensor: %d"}
    for ops in operations:
        # Add geo ops count to integral ops count for writing info.
        ops[-1] += geo_ops
        info(message[domain_type] % tuple(ops))
    return "\n".join(common) + "\n" + tensor_code

def _generate_element_tensor(integrals, sets, optimise_parameters):
    "Construct quadrature code for element tensors."

    # Prefetch formats to speed up code generation.
    f_comment    = format["comment"]
    f_ip         = format["integration points"]
    f_I          = format["ip constant"]
    f_loop       = format["generate loop"]
    f_ip_coords  = format["generate ip coordinates"]
    f_coords     = format["vertex_coordinates"]
    f_double     = format["float declaration"]
    f_decl       = format["declaration"]
    f_X          = format["ip coordinates"]
    f_C          = format["conditional"]


    # Initialise return values.
    element_code     = []
    tensor_ops_count = 0

    # TODO: KBO: The members_code was used when I generated the load_table.h
    # file which could load tables of basisfunction. This feature has not
    # been reimplemented. However, with the new design where we only
    # tabulate unique tables (and only non-zero entries) it doesn't seem to
    # be necessary. Should it be deleted?
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
            ops, coord_code = f_ip_coords(gdim, points, name, ip, r)
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
#        ip_const_ops, ip_const_code = generate_aux_constants(ip_consts, f_I,\
#                                        format["const float declaration"], True)
        ip_const_ops, ip_const_code = generate_aux_constants(ip_consts, f_I,\
                                        format["assign"], True)
        num_ops += ip_const_ops
        if ip_const_code:
            ip_code += ["", f_comment("Number of operations to compute ip constants: %d" %ip_const_ops)]
            ip_code += [format["declaration"](format["float declaration"], f_I(len(ip_consts)))]
            ip_code += ip_const_code

        # Generate code to evaluate the element tensor.
        integral_code, ops = _generate_integral_code(points, terms, sets, optimise_parameters)
        num_ops += ops
        tensor_ops_count += num_ops*points
        ip_code += integral_code

        element_code.append(f_comment\
            ("Number of operations to compute element tensor for following IP loop = %d" %(num_ops*points)) )

        # Loop code over all IPs.
        if points > 1:
            element_code += f_loop(ip_code, [(f_ip, 0, points)])
        else:
            element_code.append(f_comment("Only 1 integration point, omitting IP loop."))
            element_code += ip_code

    return (element_code, members_code, tensor_ops_count)

def _generate_functions(functions, sets):
    "Generate declarations for functions and code to compute values."

    f_comment      = format["comment"]
    f_double       = format["float declaration"]
    f_F            = format["function value"]
    f_float        = format["floating point"]
    f_decl         = format["declaration"]
    f_r            = format["free indices"][0]
    f_iadd         = format["iadd"]
    f_loop         = format["generate loop"]

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
    for loop_range, list_of_functions in function_list.items():
        function_expr = {}
        function_numbers = []
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

            # Convert function (that might be a symbol) to a simple string and save.
            function = str(function)
            function_expr[number] = function

            # Get number of operations to compute entry and add to function operations count.
            func_ops += (ops + 1)*range_i

        # Add function operations to total count
        total_ops += func_ops
        code += ["", f_comment("Total number of operations to compute function values = %d" % func_ops)]

        # Sort the functions according to name and create loop to compute the function values.
        lines = [f_iadd(f_F(n), function_expr[n]) for n in sorted(function_expr.keys())]
        code += f_loop(lines, [(f_r, 0, loop_range)]) # TODO: If loop_range == 1, this loop may be unneccessary. Not sure if it's safe to just skip it.

    return code, total_ops

def _generate_integral_code(points, terms, sets, optimise_parameters):
    "Generate code to evaluate the element tensor."

    # Prefetch formats to speed up code generation.
    f_comment       = format["comment"]
    f_mul           = format["mul"]
    f_scale_factor  = format["scale factor"]
    f_iadd          = format["iadd"]
    f_add           = format["add"]
    f_A             = format["element tensor"]
    f_loop          = format["generate loop"]
    f_B             = format["basis constant"]

    # Initialise return values.
    code = []
    num_ops = 0
    loops = {}

    # Extract sets.
    used_weights, used_psi_tables, used_nzcs, trans_set = sets

    # Loop terms and create code.
    for loop, (data, entry_vals) in terms.items():
        # If we don't have any entry values, there's no need to generate the
        # loop.
        if not entry_vals:
            continue

        # Get data.
        t_set, u_weights, u_psi_tables, u_nzcs, basis_consts = data

        # If we have a value, then we also need to update the sets of used variables.
        trans_set.update(t_set)
        used_weights.update(u_weights)
        used_psi_tables.update(u_psi_tables)
        used_nzcs.update(u_nzcs)

        # Generate code for basis constant declarations.
#        basis_const_ops, basis_const_code = generate_aux_constants(basis_consts, f_B,\
#                                        format["const float declaration"], True)
        basis_const_ops, basis_const_code = generate_aux_constants(basis_consts, f_B,\
                                        format["assign"], True)
        decl_code = []
        if basis_consts:
            decl_code = [format["declaration"](format["float declaration"], f_B(len(basis_consts)))]
        loops[loop] = [basis_const_ops, decl_code + basis_const_code]

        for entry, value, ops in entry_vals:
            # Compute number of operations to compute entry
            # (add 1 because of += in assignment).
            entry_ops = ops + 1

            # Create comment for number of operations
            entry_ops_comment = f_comment("Number of operations to compute entry: %d" % entry_ops)

            entry_code = f_iadd(f_A(entry), value)
            loops[loop][0] += entry_ops
            loops[loop][1] += [entry_ops_comment, entry_code]

    # Write all the loops of basis functions.
    for loop, ops_lines in loops.items():
        ops, lines = ops_lines
        prim_ops = functools.reduce(lambda i, j: i*j, [ops] + [l[2] for l in loop])
        # Add number of operations for current loop to total count.
        num_ops += prim_ops
        code += ["", f_comment("Number of operations for primary indices: %d" % prim_ops)]
        code += f_loop(lines, loop)

    return code, num_ops

def _tabulate_weights(quadrature_weights):
    "Generate table of quadrature weights."

    # Prefetch formats to speed up code generation.
    f_float     = format["floating point"]
    f_table     = format["static const float declaration"]
    f_sep       = format["list separator"]
    f_weight    = format["weight"]
    f_component = format["component"]
    f_group     = format["grouping"]
    f_decl      = format["declaration"]
    f_tensor    = format["tabulate tensor"]
    f_comment   = format["comment"]
    f_int       = format["int"]

    code = ["", f_comment("Array of quadrature weights.")]

    # Loop tables of weights and create code.
    for weights, points in quadrature_weights:
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

        # Tabulate the quadrature points (uncomment for different parameters).
        # 1) Tabulate the points as: p0, p1, p2, with p0 = (x0, y0, z0) etc.
        # Use f_float to format the value (enable variable precision).
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

#            # All points have the same number of coordinates.
#            num_coord = len(points[0])

#            # All points have x-coordinates.
#            xs = [f_float(p[0]) for p in points]
#            comment = "X: " + f_sep.join(xs)
#            code += [format["comment"](comment)]

#            ys = []
#            zs = []
#            # Tabulate y-coordinate if we have 2 or more coordinates.
#            if num_coord >= 2:
#                ys = [f_float(p[1]) for p in points]
#                comment = "Y: " + f_sep.join(ys)
#                code += [format["comment"](comment)]
#            # Only tabulate z-coordinate if we have 3 coordinates.
#            if num_coord == 3:
#                zs = [f_float(p[2]) for p in points]
#                comment = "Z: " + f_sep.join(zs)
#                code += [format["comment"](comment)]

        code += [""]

    return code

def _tabulate_psis(tables, used_psi_tables, inv_name_map, used_nzcs, optimise_parameters):
    "Tabulate values of basis functions and their derivatives at quadrature points."

    # Prefetch formats to speed up code generation.
    f_comment     = format["comment"]
    f_table       = format["static const float declaration"]
    f_component   = format["component"]
    f_const_uint  = format["static const uint declaration"]
    f_nzcolumns   = format["nonzero columns"]
    f_list        = format["list"]
    f_decl        = format["declaration"]
    f_tensor      = format["tabulate tensor"]
    f_new_line    = format["new line"]
    f_int         = format["int"]

    # FIXME: Check if we can simplify the tabulation
    code = []
    code += [f_comment("Value of basis functions at quadrature points.")]

    # Get list of non zero columns, if we ignore ones, ignore columns with one component.
    if optimise_parameters["ignore ones"]:
        nzcs = []
        for key, val in inv_name_map.items():
            # Check if we have a table of ones or if number of non-zero columns
            # is larger than one.
            if val[1] and len(val[1][1]) > 1 or not val[3]:
                nzcs.append(val[1])
    else:
        nzcs = [val[1] for key, val in inv_name_map.items()\
                                        if val[1]]

    # TODO: Do we get arrays that are not unique?
    new_nzcs = []
    for nz in nzcs:
        # Only get unique arrays.
        if not nz in new_nzcs:
            new_nzcs.append(nz)

    # Construct name map.
    name_map = {}
    if inv_name_map:
        for name in inv_name_map:
            if inv_name_map[name][0] in name_map:
                name_map[inv_name_map[name][0]].append(name)
            else:
                name_map[inv_name_map[name][0]] = [name]

    # Loop items in table and tabulate.
    for name in sorted(list(used_psi_tables)):
        # Only proceed if values are still used (if they're not remapped).
        vals = tables[name]
        if not vals is None:
            # Add declaration to name.
            ip, dofs = numpy.shape(vals)
            decl_name = f_component(name, [ip, dofs])

            # Generate array of values.
            value = f_tensor(vals)
            code += [f_decl(f_table, decl_name, f_new_line + value), ""]

        # Tabulate non-zero indices.
        if optimise_parameters["eliminate zeros"]:
            if name in name_map:
                for n in name_map[name]:
                    if inv_name_map[n][1] and inv_name_map[n][1] in new_nzcs:
                        i, cols = inv_name_map[n][1]
                        if not i in used_nzcs:
                            continue
                        code += [f_comment("Array of non-zero columns")]
                        value = f_list([f_int(c) for c in list(cols)])
                        name_col = f_component(f_nzcolumns(i), len(cols))
                        code += [f_decl(f_const_uint, name_col, value), ""]

                        # Remove from list of columns.
                        new_nzcs.remove(inv_name_map[n][1])
    return code
