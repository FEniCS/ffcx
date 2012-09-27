from uflacs.geometry.cellcodegen import CellGeometryCG

class TetrahedronGeometryCG(CellGeometryCG):
    def __init__(self, restriction=''):
        CellGeometryCG.__init__(self, 'tetrahedron', 3, 3, restriction)
