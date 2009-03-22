"Common base class for code generators"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-03-06 -- 2007-04-27"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Code generation modules
from dofmap import *
from finiteelement import *
from form import *

class CodeGenerator:

    def __init__(self):
        "Constructor"

        # Common methods
        self.generate_dof_map = generate_dof_map
        self.generate_finite_element = generate_finite_element
        self.generate_form = generate_form

    def generate_form_code(self, form_data, form_representation, format):
        "Generator form code according to given format"

        code = {}

        # Generate code for finite elements
        code["finite_elements"] = self.__generate_finite_elements(form_data.elements, format)

        # Generate code for dof maps
        code["dof_maps"] = self.__generate_dof_maps(form_data.dof_maps, format)

        # Generate code for cell integral
        debug("Generating code for cell integrals...")
        if form_data.num_cell_integrals:
            code["cell_integrals"] = []
        for i in range(form_data.num_cell_integrals):
            code["cell_integrals"].append(("%d" %i, self.generate_cell_integral(form_representation, i, format)))
        debug("done")

        # Generate code for exterior facet integral
        debug("Generating code for exterior facet integrals...")
        if form_data.num_exterior_facet_integrals:
            code["exterior_facet_integrals"] = []

        print "Check:", form_data.num_exterior_facet_integrals
        
        for i in range(form_data.num_exterior_facet_integrals):
            code["exterior_facet_integrals"].append(("%d"%i, self.generate_exterior_facet_integral(\
                                                       form_representation, i, format)))
        debug("done")

        # Generate code for interior facet integral
        debug("Generating code for interior facet integrals...")
        if form_data.num_interior_facet_integrals:
            code["interior_facet_integrals"] = []
        for i in range(form_data.num_interior_facet_integrals):
            code["interior_facet_integrals"].append(("%d"%i, self.generate_interior_facet_integral(\
                                                       form_representation, i, format)))
        debug("done")

        # Generate code for form
        debug("Generating code for form...")
        code["form"] = generate_form(form_data, format)
        debug("done")

        return code

    def generate_element_code(self, element_data, format):
        "Generate element code according to given format"

        code = {}

        # Generate code for finite elements
        code["finite_elements"] = self.__generate_finite_elements(element_data.elements, format)

        # Generate code for dof maps
        code["dof_maps"] = self.__generate_dof_maps(element_data.dof_maps, format)

        return code

    def __generate_finite_elements(self, elements, format):
        "Generate code for finite elements, including recursively nested sub elements"

        debug("Generating code for finite elements...")
        code = []

        # Iterate over form elements
        for i in range(len(elements)):

            # Extract sub elements
            sub_elements = self.__extract_sub_elements(elements[i], (i,))

            # Generate code for each element
            for (label, sub_element) in sub_elements:
                code += [(label, self.generate_finite_element(sub_element, format))]
                
        debug("done")
        return code

    def __generate_dof_maps(self, dof_maps, format):
        "Generate code for dof maps, including recursively nested dof maps"

        debug("Generating code for finite dof maps...")
        code = []

        # Iterate over form dof maps
        for i in range(len(dof_maps)):

            # Extract sub dof maps
            sub_dof_maps = self.__extract_sub_dof_maps(dof_maps[i], (i,))

            # Generate code for each dof map
            for (label, sub_dof_map) in sub_dof_maps:
                code += [(label, self.generate_dof_map(sub_dof_map, format))]
                
        debug("done")
        return code

    def __extract_sub_elements(self, element, parent):
        """Recursively extract sub elements as a list of tuples where
        each tuple consists of a tuple labeling the sub element and
        the sub element itself"""
        if element.num_sub_elements() == 1:
            return [(parent, element)]
        sub_elements = []
        for i in range(element.num_sub_elements()):
            sub_elements += self.__extract_sub_elements(element.sub_element(i), parent + (i,))
        return sub_elements + [(parent, element)]

    def __extract_sub_dof_maps(self, dof_map, parent):
        """Recursively extract sub dof maps as a list of tuples where
        each tuple consists of a tuple labeling the sub dof map and
        the sub dof map itself"""
        if dof_map.num_sub_dof_maps() == 1:
            return [(parent, dof_map)]
        sub_dof_maps = []
        for i in range(dof_map.num_sub_dof_maps()):
            sub_dof_maps += self.__extract_sub_dof_maps(dof_map.sub_dof_map(i), parent + (i,))
        return sub_dof_maps + [(parent, dof_map)]
