from uflacs.geometry.cellcodegen import CellGeometryCG

class HexahedronGeometryCG(CellGeometryCG):
    def __init__(self, restriction=''):
        CellGeometryCG.__init__(self, 'hexahedron', 3, 3, restriction)
