# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
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

"""Main algorithm for building the uflacs intermediate representation."""

import numpy

from ufl.classes import CellCoordinate, FacetCoordinate

from ffc.uflacs.analysis.modified_terminals import analyse_modified_terminal

from ffc.uflacs.representation.compute_expr_ir import compute_expr_ir


def build_uflacs_ir(ir, integrands, coefficient_numbering, table_provider):
    uflacs_ir = {}

    # { ufl coefficient: count }
    uflacs_ir["coefficient_numbering"] = coefficient_numbering

    # { num_points: expr_ir for one integrand }
    uflacs_ir["expr_irs"] = {}

    # Build the core uflacs expression ir for each num_points/integrand
    # TODO: Better to compute joint IR for all integrands
    #       and deal with num_points later?
    #       I.e. common_expr_ir = compute_common_expr_ir(integrands)
    #       If we want to adjoint quadrature rules for subterms
    #       automatically anyway, num_points should be advisory.
    #       For now, expecting multiple num_points to be rare.
    for num_points in sorted(integrands.keys()):
        expr_ir = compute_expr_ir(integrands[num_points])

        uflacs_ir["expr_irs"][num_points] = expr_ir

        # Build set of modified terminal ufl expressions
        V = expr_ir["V"]
        modified_terminals = [analyse_modified_terminal(V[i])
                              for i in expr_ir["modified_terminal_indices"]]
        terminal_data = modified_terminals + expr_ir["modified_arguments"]

        # FIXME: For custom integrals, skip table building but set up
        # the necessary table names and classname mappings instead

        # FIXME: Want table information earlier, even before scalar
        # rebuilding! Must split compute_expr_ir to achieve this.
        unique_tables, mt_table_ranges, table_types = \
          table_provider.build_optimized_tables(num_points, ir["entitytype"], terminal_data)


        # Figure out if we need to access CellCoordinate to
        # avoid generating quadrature point table otherwise
        if ir["integral_type"] == "cell":
            expr_ir["need_points"] = any(isinstance(mt.terminal, CellCoordinate)
                                         for mt in modified_terminals)
        elif ir["integral_type"] in ("interior_facet", "exterior_facet"):
            expr_ir["need_points"] = any(isinstance(mt.terminal, FacetCoordinate)
                                         for mt in modified_terminals)
        else:
            expr_ir["need_points"] = False


        # Ordered table data
        terminal_table_ranges = [mt_table_ranges.get(mt) for mt in terminal_data]

        # Split into arguments and other terminals before storing in expr_ir
        # TODO: Some tables are associated with num_points, some are not
        #       (i.e. piecewise constant, averaged and x0).
        #       It will be easier to deal with that if we can join
        #       the expr_ir for all num_points as mentioned above.
        n = len(expr_ir["modified_terminal_indices"])
        m = len(expr_ir["modified_arguments"])
        assert len(terminal_data) == n + m
        assert len(terminal_table_ranges) == n + m
        expr_ir["modified_terminal_table_ranges"] = terminal_table_ranges[:n]
        expr_ir["modified_argument_table_ranges"] = terminal_table_ranges[n:]

        # Store table data in V indexing, this is used in integralgenerator
        expr_ir["table_ranges"] = numpy.empty(len(V), dtype=object)
        expr_ir["table_ranges"][expr_ir["modified_terminal_indices"]] = \
            expr_ir["modified_terminal_table_ranges"]

        # FIXME: Drop tables for Real and DG0 elements (all 1.0 for each dof)

        # FIXME: Replace coefficients with empty dofrange with zero (which are these?)
        # FIXME: Propagate constants

        # Drop factorization terms where table dof range is
        # empty for any of the modified arguments
        AF = expr_ir["argument_factorization"]
        MATR = expr_ir["modified_argument_table_ranges"]
        for mas in list(AF.keys()):
            for j in mas:
                dofrange = MATR[j][1:3]
                if dofrange[0] == dofrange[1]:
                    del AF[mas]
                    break
        # FIXME: Propagate dependencies back to remove expressions
        # not used anymore after dropping factorization terms

        # Drop tables not referenced from modified terminals
        # and and tables of zeros and ones
        used_table_names = set()
        for tabledata in terminal_table_ranges:
            if tabledata is not None:
                name, begin, end = tabledata
                if table_types[name] not in ("zeros", "ones"):
                    used_table_names.add(name)
        if None in used_table_names:
            used_table_names.remove(None)
        unique_tables = { name: unique_tables[name] for name in used_table_names }

        # Store the tables and ranges
        expr_ir["table_types"] = table_types
        expr_ir["unique_tables"] = unique_tables

    return uflacs_ir
