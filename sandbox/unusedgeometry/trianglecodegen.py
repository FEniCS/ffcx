from uflacs.geometry.cellcodegen import CellGeometryCG

class TriangleGeometryCG(CellGeometryCG):
    def __init__(self, restriction=''):
        CellGeometryCG.__init__(self, 'triangle', 2, 2, restriction)
