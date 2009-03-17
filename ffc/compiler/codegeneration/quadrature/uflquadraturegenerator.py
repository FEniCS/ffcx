"Code generator for quadrature representation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2008-09-08"
__copyright__ = "Copyright (C) 2007-2008 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg 2007

# Python modules
from numpy import shape

# FFC common modules
#from ffc.common.constants import *
#from ffc.common.utils import *
from ffc.common.debug import *

# FFC fem modules
from ffc.fem.finiteelement import FiniteElement as FIATFiniteElement
from ffc.fem.vectorelement import VectorElement as FIATVectorElement
from ffc.fem.mixedelement import MixedElement as FIATMixedElement

# FFC language modules
from ffc.compiler.language.integral import Integral as FFCIntegral
#from ffc.compiler.language.index import *
#from ffc.compiler.language.restriction import *

# FFC code generation modules
#from ffc.compiler.codegeneration.common.codegenerator import *
from ffc.compiler.codegeneration.common.utils import *
from ffc.compiler.codegeneration.common.evaluatebasis import IndentControl

from ffc.fem.createelement import create_element

# FFC tensor representation modules
#from ffc.compiler.representation.tensor.multiindex import *

# Utility and optimisation functions for quadraturegenerator
from quadraturegenerator_utils import generate_loop
from uflquadraturegenerator_utils import generate_code, QuadratureTransformer
#from quadraturegenerator_optimisation import *
#import reduce_operations

# FFC format modules
from ffc.compiler.format.removeunused import remove_unused

# UFL modules
from ufl.classes import FiniteElement, MixedElement, VectorElement, FiniteElementBase, Measure
#from ufl.algorithms import *
from ufl.algorithms.analysis import extract_basis_functions, extract_elements, extract_unique_elements
#from ufl.algorithms.transformations import *
from ufl.algorithms.printing import tree_format

#class QuadratureGenerator(CodeGenerator):
class QuadratureGenerator:
    "Code generator for quadrature representation"

    def __init__(self):
        "Constructor"

        # TODO: Set this through OPTIONS
        self.optimise_options = {"non zero columns": True,
                                 "ignore ones": False,
                                 "remove zero terms": False,
                                 "simplify expressions": False,
                                 "ignore zero tables": False}

        self.reset_code = ""
        self.reset_code_restricted = ""

    def generate_integrals(self, form_representation, format):
        "Generate code for all integrals."

        code = {}

        # Set represenation
        code["representation"] = "quadrature"

        # Generate code for cell integrals
        code.update(self.generate_cell_integrals(form_representation, format))

        # Generate code for exterior facet integrals
        code.update(self.generate_exterior_facet_integrals(form_representation, format))
        
        # Generate code for interior facet integrals
        code.update(self.generate_interior_facet_integrals(form_representation, format))

        return code

    def generate_cell_integrals(self, form_representation, format):
        code = {}
        if not form_representation.cell_integrals:
            return code

        # Create transformer
        transformer = QuadratureTransformer(form_representation, Measure.CELL,\
                                            self.optimise_options, format)

        # Generate code for cell integral
        debug("Generating code for cell integrals using quadrature representation...")
        for subdomain, integrals in form_representation.cell_integrals.items():
            transformer.reset()
            code[("cell_integral", subdomain)] =\
                 self.generate_cell_integral(form_representation, transformer, integrals, format)
        debug("done")
        return code

    def generate_exterior_facet_integrals(self, form_representation, format):
        code = {}
        if not form_representation.exterior_facet_integrals:
            return code

        # Create transformer
        transformer = QuadratureTransformer(form_representation, Measure.EXTERIOR_FACET,\
                                            self.optimise_options, format)

        # Generate code for cell integral
        debug("Generating code for exterior facet integrals using quadrature representation...")
        for subdomain, integrals in form_representation.exterior_facet_integrals.items():
            transformer.reset()
            code[("exterior_facet_integral", subdomain)] =\
                 self.generate_exterior_facet_integral(form_representation, transformer, integrals, format)

        debug("done")
        return code

    def generate_interior_facet_integrals(self, form_representation, format):
        code = {}
        if not form_representation.interior_facet_integrals:
            return code

        # Create transformer
        transformer = QuadratureTransformer(form_representation, Measure.INTERIOR_FACET,\
                                            self.optimise_options, format)

        # Generate code for cell integral
        debug("Generating code for interior facet integrals using quadrature representation...")
        for subdomain, integrals in form_representation.interior_facet_integrals.items():
            transformer.reset()
            code[("interior_facet_integral", subdomain)] =\
                 self.generate_interior_facet_integral(form_representation, transformer, integrals, format)

        debug("done")
        return code

    def generate_cell_integral(self, form_representation, transformer, integrals, format):
        """Generate dictionary of code for cell integrals on a given subdomain
        from the given form representation according to the given format."""

        # Object to control the code indentation
        Indent = IndentControl()

        debug("")
#        print "\nQG, cell_integral, integrals:\n", integrals

        # FIXME: Get one of the elements, they should all be defined on the same Cell?
        # TODO: Is it faster/better to just generate it on the fly?
#        fiat_element = form_representation.fiat_elements_map[list(extract_unique_elements(integrals[0]))[0]]

        # Update treansformer with facets
        transformer.update_facets(None, None)

        # Generate element code + set of used geometry terms
        element_code, members_code, num_ops =\
          self.__generate_element_tensor(form_representation, transformer,\
                                         integrals, Indent, format)

        # Get Jacobian snippet
        # FIXME: This will most likely have to change if we support e.g., 2D elements in 3D space
        jacobi_code = [format["generate jacobian"](transformer.geo_dim, FFCIntegral.CELL)]

        # Remove unused declarations
        code = self.__remove_unused(jacobi_code, transformer.trans_set, format)

        # After we have generated the element code we know which psi tables and
        # weights will be used so we can tabulate them.

        # Tabulate weights at quadrature points
        code += self.__tabulate_weights(transformer, Indent, format)

        # Tabulate values of basis functions and their derivatives.
        code += self.__tabulate_psis(transformer, Indent, format)

        # If needed, compute the code to reset the element tensor
        # FIXME: It should be OK to pick first?
        if not self.reset_code:
            self.reset_code = self.__reset_element_tensor(integrals.items()[0][1], transformer, Indent, format, False)

        # Add element code
        code += ["", format["comment"]("Compute element tensor using UFL quadrature representation"),\
                 format["comment"]("Optimisations: %s" % ", ".join([str(i) for i in self.optimise_options.items()])),\
                 format["comment"]("Total number of operations to compute element tensor (from this point): %d" %num_ops)]
        code += element_code
        debug("Number of operations to compute tensor: %d" % num_ops)

        return {"tabulate_tensor": code, "members": members_code}


    def generate_exterior_facet_integral(self, form_representation, transformer, integrals, format):
        """Generate dictionary of code for exterior facet integral from the given
        form representation according to the given format"""

        # Object to control the code indentation
        Indent = IndentControl()

        debug("")
#        print "\nQG, exterior_facet_integral, integral:\n", integrals

        # Prefetch formats to speed up code generation
        format_comment      = format["comment"]
        format_block_begin  = format["block begin"]
        format_block_end    = format["block end"]

        # FIXME: Get one of the elements, they should all be defined on the same Cell?
        ffc_element = create_element(list(extract_unique_elements(integrals.items()[0][1]))[0])
        num_facets = ffc_element.num_facets()
        cases = [None for i in range(num_facets)]
        for i in range(num_facets):

            # Update treansformer with facets
            transformer.update_facets(i, None)
            
            case = [format_block_begin]
            c, members_code, num_ops =\
                self.__generate_element_tensor(form_representation, transformer,\
                                               integrals, Indent, format)

            case += [format_comment("Total number of operations to compute element tensor (from this point): %d" %num_ops)] + c
            case += [format_block_end]
            cases[i] = case
            debug("Number of operations to compute tensor for facet %d: %d" % (i, num_ops))

        # Get Jacobian snippet
        jacobi_code = [format["generate jacobian"](transformer.geo_dim, FFCIntegral.EXTERIOR_FACET)]

        # Remove unused declarations
        common = self.__remove_unused(jacobi_code, transformer.trans_set, format)

        # After we have generated the element code we know which psi tables and
        # weights will be used so we can tabulate them.

        # Tabulate weights at quadrature points
        common += self.__tabulate_weights(transformer, Indent, format)

        # Tabulate values of basis functions and their derivatives.
        common += self.__tabulate_psis(transformer, Indent, format)

        # If needed, compute the code to reset the element tensor
        # FIXME: It should be OK to pick first?
        if not self.reset_code:
            self.reset_code = self.__reset_element_tensor(integrals.items()[0][1], transformer, Indent, format, False)

        # Add element code
        common += ["", format["comment"]("Compute element tensor using UFL quadrature representation"),\
                 format["comment"]("Optimisations: %s" % ", ".join([str(i) for i in self.optimise_options.items()]))]

        return {"tabulate_tensor": (common, cases), "members": members_code}
    
    def generate_interior_facet_integral(self, form_representation, transformer, integrals, format):
        """Generate dictionary of code for interior facet integral from the given
        form representation according to the given format"""

        # Object to control the code indentation
        Indent = IndentControl()

        debug("")
#        print "\nQG, exterior_facet_integral, integral:\n", integrals

        # Prefetch formats to speed up code generation
        format_comment      = format["comment"]
        format_block_begin  = format["block begin"]
        format_block_end    = format["block end"]

        # FIXME: Get one of the elements, they should all be defined on the same Cell?
        ffc_element = create_element(list(extract_unique_elements(integrals.items()[0][1]))[0])
        num_facets = ffc_element.num_facets()
        cases = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                # Update treansformer with facets
                transformer.update_facets(i, j)

                case = [format_block_begin]
                c, members_code, num_ops =\
                    self.__generate_element_tensor(form_representation, transformer,\
                                                   integrals, Indent, format)
                case += [format_comment("Total number of operations to compute element tensor (from this point): %d" %num_ops)] + c
                case += [format_block_end]
                cases[i][j] = case
                debug("Number of operations to compute tensor for facets (%d, %d): %d" % (i, j, num_ops))

        # Get Jacobian snippet
        jacobi_code = [format["generate jacobian"](transformer.geo_dim, FFCIntegral.INTERIOR_FACET)]

        # Remove unused declarations
        common = self.__remove_unused(jacobi_code, transformer.trans_set, format)

        # Tabulate weights at quadrature points
        common += self.__tabulate_weights(transformer, Indent, format)

        # Tabulate values of basis functions and their derivatives.
        common += self.__tabulate_psis(transformer, Indent, format)

        # If needed, compute the code to reset the element tensor
        # FIXME: It should be OK to pick first?
        if not self.reset_code_restricted:
            self.reset_code_restricted = self.__reset_element_tensor(integrals.items()[0][1], transformer, Indent, format, True)

        # Add element code
        common += ["", format["comment"]("Compute element tensor using UFL quadrature representation"),\
                 format["comment"]("Optimisations: %s" % ", ".join([str(i) for i in self.optimise_options.items()]))]

        return {"tabulate_tensor": (common, cases), "constructor":"// Do nothing", "members":members_code}

    def __generate_element_tensor(self, form_representation, transformer, integrals, Indent, format):
        "Construct quadrature code for element tensors"

        # Prefetch formats to speed up code generation
        format_comment      = format["comment"]
        format_ip           = format["integration points"]
        format_G            = format["geometry tensor"]
        format_const_float  = format["const float declaration"]
        format_weight       = format["weight"]
        format_scale_factor = format["scale factor"]

        # Initialise return values.
        # FIXME: The members_code was used when I generated the load_table.h
        # file which could load tables of basisfunction. This feature has not
        # been reimplemented. However, with the new design where we only
        # tabulate unique tables (and only non-zero entries) it doesn't seem to
        # be necessary.
        members_code     = ""
        element_code     = []
        tensor_ops_count = 0

        # Since the form_representation holds common tables for all integrals,
        # I need to keep track of which tables are actually used for the current
        # subdomain and then only tabulate those.
        # The same holds true for the quadrature weights.
        # I therefore need to generate the actual code to compute the element
        # tensors first, and then create the auxiliary code.

        transformer.disp()

        # We receive a dictionary {num_points: integral,}
        # Loop points and integrals
        for points, integral in integrals.items():

#            print "\nIntegral: ", integral
#            print "\nIntegral: ", str(integral)
#            print "\nIntegral tree_format: ", tree_format(integral)

            ip_code = ["", Indent.indent(format_comment\
                ("Loop quadrature points for integral: %s" % repr(integral)))]
#                ("Loop quadrature points for integral: %s" % integral.__repr__()))]

            # Update transformer to the current number of quadrature points
            transformer.update_points(points)

            # Generate code for all terms according to optimisation level
            integral_code, num_ops =\
                generate_code(integral.integrand(), transformer, Indent, format)

            # Get number of operations to compute entries for all terms when
            # looping over all IPs and update tensor count
            num_operations = num_ops*points
            tensor_ops_count += num_operations

            ip_code.append(format_comment\
                ("Number of operations to compute element tensor for following IP loop = %d" %(num_operations)) )

            # Loop code over all IPs
            if points > 1:
                ip_code += generate_loop(integral_code, [(format_ip, 0, points)], Indent, format)
            else:
                ip_code.append(format_comment("Only 1 integration point, omitting IP loop."))
                ip_code += integral_code

            # Add integration point code element code
            element_code += ip_code

#        # Tabulate geometry code, sort according to number
#        geo_code = []
#        items = geo_terms.items()
#        items = [(int(v.replace(format_G, "")), k) for (k, v) in items]
#        items.sort()
#        items = [(k, format_G + str(v)) for (v, k) in items]
#        geo_ops = 0
#        for key, val in items:
#            declaration = red_ops(exp_ops(key, format), format)
#            # Get number of operations needed to compute geometry declaration
#            geo_ops += count_ops(declaration, format)
#            geo_code += [(format_const_float + val, declaration)]
#        geo_code.append("")
#        if geo_ops:
#            geo_code = ["", format_comment("Number of operations to compute geometry constants = %d" % geo_ops)] + geo_code
#        else:
#            geo_code = [""] + geo_code

#        # Add operation count
#        tensor_ops_count += geo_ops

#        element_code = geo_code + element_code

#        return (tabulate_code + element_code, members_code, tensor_ops_count)
        return (element_code, members_code, tensor_ops_count)


    def __tabulate_weights(self, transformer, Indent, format):
        "Generate table of quadrature weights"

        # Prefetch formats to speed up code generation
        format_float    = format["floating point"]
        format_table    = format["table declaration"]
        format_block    = format["block"]
        format_sep      = format["separator"]
        format_weight   = format["weight"]
        format_array    = format["array access"]

        code = ["", Indent.indent(format["comment"]("Array of quadrature weights"))]

        # Loop tables of weights and create code
        for points in transformer.used_weights:
            weights = transformer.quadrature_weights[points]

            # FIXME: For now, raise error if we don't have weights.
            # We might want to change this later
            if not weights.any():
                raise RuntimeError(weights, "No weights")

            # Create name and value
            name = format_table + format_weight(points)
            value = format_float(weights[0])
            if len(weights) > 1:
                name += format_array(str(points))
                value = format_block(format_sep.join([format_float(w)\
                                                      for w in weights]))
            code += [(Indent.indent(name), value), ""]

        return code

    def __tabulate_psis(self, transformer, Indent, format):
        "Tabulate values of basis functions and their derivatives at quadrature points"

        # Prefetch formats to speed up code generation
        format_comment    = format["comment"]
        format_float      = format["floating point"]
        format_block      = format["block"]
        format_table      = format["table declaration"]
        format_matrix     = format["matrix access"]
        format_array      = format["array access"]
        format_const_uint = format["static const uint declaration"]
        format_nzcolumns  = format["nonzero columns"]
        format_sep        = format["separator"]

        code = []
        # FIXME: Check if we can simplify the tabulation

        inv_name_map = transformer.name_map
        tables = transformer.unique_tables

        # Get list of non zero columns, if we ignore ones ignore columns with
        # one component
        if self.optimise_options["ignore ones"]:
            nzcs = [val[1] for key, val in inv_name_map.items()\
                                           if val[1] and len(val[1][1]) > 1]
        else:
            nzcs = [val[1] for key, val in inv_name_map.items()\
                                           if val[1]]

        # TODO: Do we get arrays that are not unique?
        new_nzcs = []
        for nz in nzcs:
            # Only get unique arrays
            if not nz in new_nzcs:
                new_nzcs.append(nz)

        # Construct name map
        name_map = {}
        if inv_name_map:
            for name in inv_name_map:
                if inv_name_map[name][0] in name_map:
                    name_map[inv_name_map[name][0]].append(name)
                else:
                    name_map[inv_name_map[name][0]] = [name]

        # Loop items in table and tabulate 
#        for name, vals in tables.items():
        for name in transformer.used_psi_tables:
            # Only proceed if values are still used (if they're not remapped)
            vals = tables[name]
            if not vals == None:
                # Add declaration to name
                ip, dofs = shape(vals)
                decl_name = format_table + name + format_matrix(ip, dofs)

                # Generate array of values
                value = tabulate_matrix(vals, format)
                code += ["", (Indent.indent(decl_name), Indent.indent(value))]

            # Tabulate non-zero indices
            if self.optimise_options["non zero columns"]:
                if name in name_map:
                    for n in name_map[name]:
                        if inv_name_map[n][1] and inv_name_map[n][1] in new_nzcs:
                            code += [Indent.indent(format_comment("Array of non-zero columns") )]
                            i, cols = inv_name_map[n][1]
                            value = format_block(format_sep.join(["%d" %c for c in list(cols)]))
                            name_col = format_const_uint + format_nzcolumns(i) + format_array(len(cols))
                            code += [(Indent.indent(name_col), value)]
                            # Remove from list of columns
                            new_nzcs.remove(inv_name_map[n][1])

        # Tabulate remaining non-zero columns for tables that might have been deleted
        new_nzcs = [nz for nz in new_nzcs if nz and len(nz[1]) > 1]
        if self.optimise_options["non zero columns"] and new_nzcs:
            code += [Indent.indent(format_comment("Array of non-zero columns for arrays that might have been deleted (on purpose)") )]
            for i, cols in new_nzcs:
                value = format_block(format_sep.join(["%d" %c for c in list(cols)]))
                name_col = format_const_uint + format_nzcolumns(i) + format_array(len(cols))
                code += [(Indent.indent(name_col), value)]

        return code

    def __reset_element_tensor(self, integral, transformer, Indent, format, interior_integrals):
        "Reset the entries of the local element tensor"

        code = [""]

        # Comment
        code.append(Indent.indent(format["comment"]\
                ("Reset values of the element tensor block")))

        # Get basisfunctions
        basis = extract_basis_functions(integral)

        # Create FIAT elements for each basisfunction. There should be one and
        # only one element per basisfunction so it is OK to pick first.
        elements = [create_element(extract_elements(b)[0]) for b in basis]
        #print "\nQG, reset_element_tensor, Elements:\n", elements

        # Create the index range for resetting the element tensor by
        # multiplying the element space dimensions
        # FIXME: I don't think restricted basisfunctions on e.g., interior
        # facet integrals are handled correctly yet. Should multiply by 2
        # somewhere.
        index_range = 1
        restrict_range = 1
        if interior_integrals:
            restrict_range = 2
        for element in elements:
            index_range *= restrict_range*element.space_dimension()

        # Create loop
        # FIXME: It is general to create a loop, however, for a functional it
        # is not strictly needed. On the other hand, it is not expected to have
        # any influence on the runtime performance.
        value = format["floating point"](0.0)
        name =  format["element tensor quad"] + format["array access"](format["first free index"])
        lines = [(name, value)]
        code += generate_loop(lines, [(format["first free index"], 0, index_range)], Indent, format)
        code.append("")

        return code

    def __remove_unused(self, code, trans_set, format):
        "Remove unused variables so that the compiler will not complain"

        # Normally, the removal of unused variables should happen at the
        # formatting stage, but since the code for the tabulate_tensor()
        # function may grow to considerable size, we make an exception and
        # remove unused variables here when we know the names of the used
        # variables. No searching necessary and much, much, much faster.

        if code:
            # Generate body of code, using the format
            lines = format["generate body"](code)

            # Generate auxiliary code line that uses all members of the set
            # (to trick remove_unused)
            line_set = format["add equal"]("A", format["multiply"](trans_set))
            lines += "\n" + line_set

            # Remove unused Jacobi declarations
            code = remove_unused(lines)

            # Delete auxiliary line
            code = code.replace("\n" + line_set, "")

            return [code]
        else:
            return code
