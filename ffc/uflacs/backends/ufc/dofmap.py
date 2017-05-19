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


from ffc.uflacs.backends.ufc.generator import ufc_generator
from ffc.uflacs.backends.ufc.utils import generate_return_new_switch, generate_return_sizet_switch, generate_return_bool_switch


class ufc_dofmap(ufc_generator):
    "Each function maps to a keyword in the template. See documentation of ufc_generator."
    def __init__(self):
        ufc_generator.__init__(self, "dofmap")


    def topological_dimension(self, L, topological_dimension):
        return L.Return(topological_dimension)


    def needs_mesh_entities(self, L, needs_mesh_entities):
        "needs_mesh_entities is a list of num dofs per entity."
        return generate_return_bool_switch(L, "d", needs_mesh_entities, False)


    def global_dimension(self, L, ir):
        # A list of num dofs per entity
        num_entity_dofs, num_global_dofs = ir["global_dimension"]

        # Input array
        entities = L.Symbol("num_global_entities")

        # Accumulate number of dofs per entity times given number of entities in mesh
        entity_dofs = []
        for dim, num in enumerate(num_entity_dofs):
            if num == 1:
                entity_dofs.append(entities[dim])
            elif num > 1:
                entity_dofs.append(num * entities[dim])

        if entity_dofs:
            dimension = L.Sum(entity_dofs)
        else:
            dimension = 0

        # Add global dofs if any
        if num_global_dofs:
            dimension = dimension + num_global_dofs

        return L.Return(dimension)


    def num_global_support_dofs(self, L, num_global_support_dofs):
        return L.Return(num_global_support_dofs)


    def num_element_support_dofs(self, L, num_element_support_dofs):
        return L.Return(num_element_support_dofs)


    def num_element_dofs(self, L, num_element_dofs):
        return L.Return(num_element_dofs)


    def num_facet_dofs(self, L, num_facet_dofs):
        return L.Return(num_facet_dofs)


    def num_entity_dofs(self, L, num_entity_dofs):
        return generate_return_sizet_switch(L, "d", num_entity_dofs, 0)


    def num_entity_closure_dofs(self, L, num_entity_closure_dofs):
        return generate_return_sizet_switch(L, "d", num_entity_closure_dofs, 0)


    def tabulate_dofs(self, L, ir):
        # Input arguments
        entity_indices = L.Symbol("entity_indices")
        num_mesh_entities = L.Symbol("num_global_entities")

        # Output arguments
        dofs_variable = L.Symbol("dofs")

        ir = ir["tabulate_dofs"]
        if ir is None:
            # Special case for SpaceOfReals, ir returns None
            # FIXME: This is how the old code did in this case,
            # is it correct? Is there only 1 dof?
            # I guess VectorElement(Real) handled elsewhere?
            code = [L.Assign(dofs_variable[0], 0)]
            return L.StatementList(code)

        # Extract representation
        #(dofs_per_element, num_dofs_per_element, need_offset, fakes) = ir # names from original ffc code
        (subelement_dofs, num_dofs_per_subelement, need_offset, is_subelement_real) = ir

        #len(is_subelement_real) == len(subelement_dofs) == len(num_dofs_per_subelement) == number of combined (flattened) subelements?
        #is_subelement_real[subelement_number] is bool, is subelement Real?
        #num_dofs_per_subelement[subelement_number] is int, number of dofs for each (flattened) subelement
        #subelement_dofs[subelement_number] is list of entity_dofs for each (flattened) subelement
        #len(entity_dofs) == tdim+1
        #len(entity_dofs[cell_entity_dim]) == "num_cell_entities[dim]" (i.e. 3 vertices, 3 edges, 1 face for triangle)
        #len(entity_dofs[cell_entity_dim][entity_index[dim][k]]) == number of dofs for this entity
        #num = entity_dofs[dim]

        #entity_dofs =? entity_dofs[d][i][:] # dofs on entity (d,i)

        # Collect code pieces in list
        code = []

        # Declare offset if needed
        if need_offset:
            offset = L.Symbol("offset")
            code.append(L.VariableDecl("std::size_t", offset, value=0))
        else:
            offset = 0

        # Generate code for each element
        subelement_offset = 0
        for (subelement_index, entity_dofs) in enumerate(subelement_dofs):

            # Handle is_subelement_real (Space of reals)
            if is_subelement_real[subelement_index]:
                assert num_dofs_per_subelement[subelement_index] == 1
                code.append(L.Assign(dofs_variable[subelement_offset], offset))
                if need_offset:
                    code.append(L.AssignAdd(offset, 1))
                subelement_offset += 1
                continue

            # Generate code for each degree of freedom for each dimension
            for (cell_entity_dim, dofs_on_cell_entity) in enumerate(entity_dofs):
                num_dofs_per_mesh_entity = len(dofs_on_cell_entity[0])
                assert all(num_dofs_per_mesh_entity == len(dofs)
                           for dofs in dofs_on_cell_entity)

                # Ignore if no dofs for this dimension
                if num_dofs_per_mesh_entity == 0:
                    continue

                # For each cell entity of this dimension
                for (cell_entity_index, dofs) in enumerate(dofs_on_cell_entity):
                    # dofs is a list of the local dofs that live on this cell entity

                    # find offset for this particular mesh entity
                    entity_offset = len(dofs) * entity_indices[cell_entity_dim, cell_entity_index]

                    for (j, dof) in enumerate(dofs):
                        # dof is the local dof index on the subelement
                        # j is the local index of dof among the dofs
                        # on this particular cell/mesh entity
                        local_dof_index = subelement_offset + dof
                        global_dof_index = offset + entity_offset + j
                        code.append(L.Assign(dofs_variable[local_dof_index], global_dof_index))

                # Update offset corresponding to mesh entity:
                if need_offset:
                    value = num_dofs_per_mesh_entity * num_mesh_entities[cell_entity_dim]
                    code.append(L.AssignAdd(offset, value))

            subelement_offset += num_dofs_per_subelement[subelement_index]

        return L.StatementList(code)


    def tabulate_facet_dofs(self, L, ir):
        all_facet_dofs = ir["tabulate_facet_dofs"]

        # Input arguments
        facet = L.Symbol("facet")
        dofs = L.Symbol("dofs")

        # For each facet, copy all_facet_dofs[facet][:] into output argument array dofs[:]
        cases = []
        for f, single_facet_dofs in enumerate(all_facet_dofs):
            assignments = [L.Assign(dofs[i], dof)
                           for (i, dof) in enumerate(single_facet_dofs)]
            if assignments:
                cases.append((f, L.StatementList(assignments)))
        if cases:
            return L.Switch(facet, cases, autoscope=False)
        else:
            return L.NoOp()


    def tabulate_entity_dofs(self, L, ir):
        entity_dofs, num_dofs_per_entity = ir["tabulate_entity_dofs"]

        # Output argument array
        dofs = L.Symbol("dofs")

        # Input arguments
        d = L.Symbol("d")
        i = L.Symbol("i")

        # TODO: Removed check for (d <= tdim + 1)
        tdim = len(num_dofs_per_entity) - 1

        # Generate cases for each dimension:
        all_cases = []
        for dim in range(tdim + 1):

            # Ignore if no entities for this dimension
            if num_dofs_per_entity[dim] == 0:
                continue

            # Generate cases for each mesh entity
            cases = []
            for entity in range(len(entity_dofs[dim])):
                casebody = []
                for (j, dof) in enumerate(entity_dofs[dim][entity]):
                    casebody += [L.Assign(dofs[j], dof)]
                cases.append((entity, L.StatementList(casebody)))

            # Generate inner switch
            # TODO: Removed check for (i <= num_entities-1)
            inner_switch = L.Switch(i, cases, autoscope=False)
            all_cases.append((dim, inner_switch))

        if all_cases:
            return L.Switch(d, all_cases, autoscope=False)
        else:
            return L.NoOp()


    def tabulate_entity_closure_dofs(self, L, ir):
        # Extract variables from ir
        entity_closure_dofs, entity_dofs, num_dofs_per_entity = \
            ir["tabulate_entity_closure_dofs"]

        # Output argument array
        dofs = L.Symbol("dofs")

        # Input arguments
        d = L.Symbol("d")
        i = L.Symbol("i")

        # TODO: Removed check for (d <= tdim + 1)
        tdim = len(num_dofs_per_entity) - 1

        # Generate cases for each dimension:
        all_cases = []
        for dim in range(tdim + 1):
            num_entities = len(entity_dofs[dim])

            # Generate cases for each mesh entity
            cases = []
            for entity in range(num_entities):
                casebody = []
                for (j, dof) in enumerate(entity_closure_dofs[(dim, entity)]):
                    casebody += [L.Assign(dofs[j], dof)]
                cases.append((entity, L.StatementList(casebody)))

            # Generate inner switch unless empty
            if cases:
                # TODO: Removed check for (i <= num_entities-1)
                inner_switch = L.Switch(i, cases, autoscope=False)
                all_cases.append((dim, inner_switch))

        if all_cases:
            return L.Switch(d, all_cases, autoscope=False)
        else:
            return L.NoOp()


    def num_sub_dofmaps(self, L, num_sub_dofmaps):
        return L.Return(num_sub_dofmaps)


    def create_sub_dofmap(self, L, ir):
        classnames = ir["create_sub_dofmap"]
        return generate_return_new_switch(L, "i", classnames, factory=ir["jit"])
