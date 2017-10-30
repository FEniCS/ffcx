

# CURRENTLY UNUSED CODE, CONVERT TO CREATE NEW ufc::finite_element::apply_mappings functions


def _generate_apply_mapping_to_computed_values(L, dof_data):
    mapping = dof_data["mapping"]
    num_components = dof_data["num_components"]
    reference_offset = dof_data["reference_offset"]
    physical_offset = dof_data["physical_offset"]

    #physical_values[num_points][num_dofs][physical_value_size]
    #reference_values[num_points][num_dofs][reference_value_size]

    code = []

    # FIXME: Define values numbering
    if mapping == "affine":
        # Just copy values
        if num_components == 1:
            code += [
                L.ForRange(ip, 0, num_points, body=
                    L.Assign(physical_values[ip][physical_offset],
                             reference_values[ip][reference_offset])),
                ]
        else:
            code += [
                L.ForRange(ip, 0, num_points, body=
                    L.ForRange(k, 0, num_components, body=
                        L.Assign(physical_values[ip][physical_offset + k],
                                 reference_values[ip][reference_offset + k]))),
                ]
        return code

    else:
        FIXME

    return code


def __ffc_implementation_of__generate_apply_mapping_to_computed_values(L):
    # Apply transformation if applicable.
    mapping = dof_data["mapping"]
    num_components = dof_data["num_components"]
    reference_offset = dof_data["reference_offset"]
    physical_offset = dof_data["physical_offset"]

    if mapping == "affine":
        pass

    elif mapping == "contravariant piola":
        code += ["", f_comment("Using contravariant Piola transform to map values back to the physical element")]

        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + offset))
                 for i in range(num_components)]

        # Create names for inner product.
        basis_col = [f_tmp_ref(j) for j in range(tdim)]
        for i in range(gdim):
            # Create Jacobian.
            jacobian_row = [f_trans("J", i, j, gdim, tdim, None) for j in range(tdim)]
            # Create inner product and multiply by inverse of Jacobian.
            inner = f_group(f_inner(jacobian_row, basis_col)) # sum_j J[i,j], values[offset + j])
            value = f_mul([f_inv(f_detJ(None)), inner])
            name = f_component(f_values, i + offset)
            # FIXME: This is writing values[offset+:] = M[:,:] * values[offset+:],
            #        i.e. offset must be physical (unless there's a bug).
            #        We want to use reference offset for evaluate_reference_basis,
            #        and to make this mapping read from one array using reference offset
            #        and write to another array using physical offset.
            code += [f_assign(name, value)]

    elif mapping == "covariant piola":
        code += ["", f_comment("Using covariant Piola transform to map values back to the physical element")]
        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + offset))
                 for i in range(num_components)]
        # Create names for inner product.
        tdim = data["topological_dimension"]
        gdim = data["geometric_dimension"]
        basis_col = [f_tmp_ref(j) for j in range(tdim)]
        for i in range(gdim):
            # Create inverse of Jacobian.
            inv_jacobian_column = [f_trans("JINV", j, i, tdim, gdim, None) for j in range(tdim)]

            # Create inner product of basis values and inverse of Jacobian.
            value = f_group(f_inner(inv_jacobian_column, basis_col))
            name = f_component(f_values, i + offset)
            code += [f_assign(name, value)]

    elif mapping == "double covariant piola":
        code += ["", f_comment("Using double covariant Piola transform to map values back to the physical element")]
        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + offset))
                 for i in range(num_components)]
        # Create names for inner product.
        tdim = data["topological_dimension"]
        gdim = data["geometric_dimension"]
        basis_col = [f_tmp_ref(j) for j in range(num_components)]
        for p in range(num_components):
            # unflatten the indices
            i = p // tdim
            l = p % tdim  # noqa: E741
            # g_il = K_ji G_jk K_kl
            value = f_group(f_inner(
                [f_inner([f_trans("JINV", j, i, tdim, gdim, None)
                          for j in range(tdim)],
                         [basis_col[j * tdim + k] for j in range(tdim)])
                 for k in range(tdim)],
                [f_trans("JINV", k, l, tdim, gdim, None)
                 for k in range(tdim)]))
            name = f_component(f_values, p + offset)
            code += [f_assign(name, value)]

    elif mapping == "double contravariant piola":
        code += ["", f_comment("Pullback of a matrix-valued funciton as contravariant 2-tensor mapping values back to the physical element")]
        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + offset))
                 for i in range(num_components)]
        # Create names for inner product.
        tdim = data["topological_dimension"]
        gdim = data["geometric_dimension"]
        basis_col = [f_tmp_ref(j) for j in range(num_components)]
        for p in range(num_components):
            # unflatten the indices
            i = p // tdim
            l = p % tdim  # noqa: E741
            # g_il = (detJ)^(-2) J_ij G_jk J_lk
            value = f_group(f_inner(
                [f_inner([f_trans("J", i, j, tdim, gdim, None)
                          for j in range(tdim)],
                         [basis_col[j * tdim + k] for j in range(tdim)])
                 for k in range(tdim)],
                [f_trans("J", l, k, tdim, gdim, None)
                 for k in range(tdim)]))
            value = f_mul([f_inv(f_detJ(None)), f_inv(f_detJ(None)), value])
            name = f_component(f_values, p + offset)
            code += [f_assign(name, value)]

    else:
        error("Unknown mapping: %s" % mapping)
