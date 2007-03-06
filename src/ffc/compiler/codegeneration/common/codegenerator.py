"Common base class for code generators"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-03-06 -- 2007-02-06"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Code generation modules
from dofmap import *
from finiteelement import *
from form import *

class CodeGenerator:

    def __init__(self, form_data, form_representation, format):
        "Constructor"

        # Common methods
        self.generate_dof_map = generate_dof_map
        self.generate_finite_element = generate_finite_element
        self.generate_form = generate_form

        # Save data
        self.form_data = form_data
        self.form_representation = form_representation
        self.format = format
        
    def generate_code(self):
        "Generator code according to given format"

        code = {}

        form_data = self.form_data
        form_representation = self.form_representation
        format = self.format

        # Generate code for finite elements
        debug("Generating code for finite elements...")
        for i in range(len(form_data.elements)):
            code[("finite_element", i)] = self.generate_finite_element(form_data.elements[i], format)
        debug("done")

        # Generate code for dof maps
        debug("Generating code for dof maps...")
        for i in range(len(form_data.dof_maps)):
            code[("dof_map", i)] = self.generate_dof_map(form_data.dof_maps[i], format)
        debug("done")

        # Generate code for cell integral
        debug("Generating code for cell integrals...")
        for i in range(form_data.num_cell_integrals):
            code[("cell_integral", i)] = self.generate_cell_integral(form_representation, i, format)
        debug("done")

        # Generate code for exterior facet integral
        debug("Generating code for exterior facet integrals...")
        for i in range(form_data.num_exterior_facet_integrals):
            code[("exterior_facet_integral", i)] = self.generate_exterior_facet_integral(form_representation, i, format)
        debug("done")

        # Generate code for interior facet integral
        debug("Generating code for interior facet integrals...")
        for i in range(form_data.num_interior_facet_integrals):
            code[("interior_facet_integral", i)] = self.generate_interior_facet_integral(form_representation, i, format)
        debug("done")

        # Generate code for form
        debug("Generating code for form...")
        code["form"] = generate_form(form_data, format)
        debug("done")

        return code
