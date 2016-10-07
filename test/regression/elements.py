# -*- coding: utf-8 -*-
interval_2D = "Cell('interval', geometric_dimension=2)"
interval_3D = "Cell('interval', geometric_dimension=3)"
triangle_3D = "Cell('triangle', geometric_dimension=3)"

elements = ["FiniteElement('N1curl', triangle, 2)",
            "MixedElement([FiniteElement('Lagrange', triangle, 3), \
                           VectorElement('Lagrange', triangle, 3)['facet']])",
            "VectorElement('R', triangle, 0, 3)",
            "VectorElement('DG', %s, 1)" % interval_2D,
            "VectorElement('DG', %s, 1)" % interval_3D,
            "VectorElement('DG', %s, 1)" % triangle_3D,
            "MixedElement([VectorElement('CG', %s, 2), \
                           FiniteElement('CG', %s, 1)])" % (interval_2D,
                                                            interval_2D),
            "MixedElement([VectorElement('CG', %s, 2), \
                           FiniteElement('CG', %s, 1)])" % (interval_3D,
                                                            interval_3D),
            "MixedElement([VectorElement('CG', %s, 2), \
                           FiniteElement('CG', %s, 1)])" % (triangle_3D,
                                                            triangle_3D),
            "MixedElement([FiniteElement('RT', %s, 2), \
                           FiniteElement('BDM', %s, 1), \
                           FiniteElement('N1curl', %s, 1), \
                           FiniteElement('DG', %s, 1)])" % (triangle_3D,
                                                            triangle_3D,
                                                            triangle_3D,
                                                            triangle_3D),
            "FiniteElement('Regge', triangle, 2)",
            "MixedElement([FiniteElement('HHJ', triangle, 2), \
                           FiniteElement('CG', triangle, 3)])",
            ]
