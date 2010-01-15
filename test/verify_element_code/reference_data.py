"""
References for known values for generated element and dof map code
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU GPL version 3 or any later version"

reference = {}

# 1D
element = "FiniteElement('CG', 'interval', 1)"
ref = {
    "space_dimension": 2,
    "value_dimension": 1
    }
reference[element] = ref

element = "FiniteElement('CG', 'interval', 4)"
ref = {
    "space_dimension": 5,
    "value_dimension": 1
    }
reference[element] = ref

element = "FiniteElement('DG', 'interval', 3)"
ref = {
    "space_dimension": 4,
    "value_dimension": 1
    }
reference[element] = ref

# 2D
element = "FiniteElement('CG', 'triangle', 1)"
ref = {
    "space_dimension": 3,
    "value_dimension": 1
    }
reference[element] = ref


element = "FiniteElement('CG', 'triangle', 2)"
ref = {
    "space_dimension": 6,
    "value_dimension": 1
    }
reference[element] = ref


element = "VectorElement('CG', 'triangle', 2)"
ref = {
    "space_dimension": 12,
    "value_dimension": 2
    }
reference[element] = ref


element = "FiniteElement('RT', 'triangle', 2)"
ref = {
    "space_dimension": 8,
    "value_dimension": 2
    }
reference[element] = ref

element = "FiniteElement('Nedelec 1st kind H(curl)', 'triangle', 3)"
ref = {
    "space_dimension": 15,
    "value_dimension": 2
    }
reference[element] = ref


element = "FiniteElement('BDM', 'triangle', 1)"
ref = {
    "space_dimension": 6,
    "value_dimension": 2
    }
reference[element] = ref
