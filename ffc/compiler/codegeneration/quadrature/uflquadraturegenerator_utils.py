"Utility functions for UFL quadrature code generation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-02-09 -- 2009-02-09"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Python modules
from numpy import shape, transpose

# UFL modules
from ufl.classes import AlgebraOperator, FormArgument

# Utility and optimisation functions for quadraturegenerator
from quadraturegenerator_utils import unique_psi_tables


def create_psi_tables(tables, optimise_level, format):
    "Create names and maps for tables and non-zero entries if appropriate."

    print "\nQG-utils, psi_tables:\n", tables

    flat_tables = flatten_psi_tables(tables)
#    print "\nQG-utils, psi_tables, flat_tables:\n", flat_tables

    name_map, unique_tables = unique_psi_tables(flat_tables, optimise_level, format)
#    name_map, new_tables = unique_psi_tables(flat_tables, 1, format)

#    print "\nQG-utils, psi_tables, unique_tables:\n", unique_tables
#    print "\nQG-utils, psi_tables, name_map:\n", name_map

    return (name_map, unique_tables)

def flatten_psi_tables(tables):
    "Create a 'flat' dictionary of tables with unique names."

#    print "\nQG-utils, flatten_psi_tables:\n", tables

    flat_tables = {}
    counter = 0
    # Loop quadrature points and get element dictionary {elem: {tables}}
    for point, elem_dict in tables.items():
#        print "\nQG-utils, flatten_tables, points:\n", point
#        print "\nQG-utils, flatten_tables, elem_dict:\n", elem_dict

        # Loop all elements and get all their tables
        for elem, elem_tables in elem_dict.items():
#            print "\nQG-utils, flatten_tables, elem:\n", elem
#            print "\nQG-utils, flatten_tables, elem[0].value_rank():\n", elem[0].value_rank()
#            print "\nQG-utils, flatten_tables, elem_tables:\n", elem_tables
            # If the element value rank != 0, we must loop the components
            # before the derivatives
            if elem[0].value_rank() != 0:
                for num_comp, comp in enumerate(elem_tables):
                    for num_deriv in comp:
                        for derivs, psi_table in num_deriv.items():
#                            print "\nQG-utils, flatten_tables, derivs:\n", derivs
#                            print "\nQG-utils, flatten_tables, psi_table:\n", psi_table
                            # Verify shape of basis (can be omitted for speed
                            # if needed I think)
                            if shape(psi_table) != 2 and shape(psi_table)[1] != point:
                                raise RuntimeError(psi_table, "Something is wrong with this table")

                            name = generate_psi_name(counter, elem[1], num_comp, derivs)
#                            print "Name: ", name
                            if name in flat_tables:
                                raise RuntimeError(name, "Name is not unique, something is wrong")
                            flat_tables[name] = transpose(psi_table)
            else:
                for num_deriv in elem_tables:
                    for derivs, psi_table in num_deriv.items():
#                        print "\nQG-utils, flatten_tables, derivs:\n", derivs
#                        print "\nQG-utils, flatten_tables, psi_table:\n", psi_table
                        # Verify shape of basis (can be omitted for speed
                        # if needed I think)
                        if shape(psi_table) != 2 and shape(psi_table)[1] != point:
                            raise RuntimeError(psi_table, "Something is wrong with this table")
                        name = generate_psi_name(counter, elem[1], None, derivs)
#                        print "Name: ", name
                        if name in flat_tables:
                            raise RuntimeError(name, "Name is not unique, something is wrong")
                        flat_tables[name] = transpose(psi_table)
            counter += 1

    return flat_tables

def generate_psi_name(counter, restriction, component, derivatives):
    """Generate a name for the psi table of the form:
    FE#_R#_C#_D###, where '#' will be an integer value.

    FE  - is a simple counter to distinguish the various basis, it will be
          assigned in an arbitrary fashion.

    R   - denotes restrictions if applicable, 0 or 1.

    C   - is the component number if any (this does not yet take into account
          tensor valued functions)

    D   - is the number of derivatives in each spatial direction if any. If the
          element is defined in 3D, then D012 means d^3(*)/dydz^2."""

    name = "FE%d" %counter
    if restriction:
        name += "_R%d" % restriction
    if component != None:
        name += "_C%d" % component
    if any(derivatives):
        name += "_D" + "".join([str(d) for d in derivatives])

    return name

def generate_code(raw_terms, geo_terms, optimise_level, Indent, format):
    """Generate code from a UFL integral type.
    This function implements the different optimisation strategies."""

    return ("", 0)






