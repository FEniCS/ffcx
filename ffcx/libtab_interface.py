
import ufl
import libtab


def _create_element(element):
    el = libtab.create_element(element.family(),
                               element.cell().cellname(),
                               element.degree())
    return el


def create_libtab_elements(ufl_element):

    elements = []

    print('FiniteElement ', isinstance(ufl_element, ufl.FiniteElement))
    print('VectorElement ', isinstance(ufl_element, ufl.VectorElement))
    print('MixedElement ', isinstance(ufl_element, ufl.MixedElement))

    def rextract(els):
        for e in els:
            if isinstance(e, ufl.MixedElement) \
               and not isinstance(e, ufl.VectorElement) \
               and not isinstance(e, ufl.TensorElement):
                rextract(e.sub_elements())
            else:
                elements.append(e)

    if isinstance(ufl_element, ufl.FiniteElement):
        elements.append(_create_element(ufl_element))
    elif isinstance(ufl_element, ufl.MixedElement):
        rextract(ufl_element.sub_elements())
        elements = list(map(_create_element, elements))

    return elements
