from uflacs.geometry.cellcodegen import CellGeometryCG

class QuadrilateralGeometryCG(CellGeometryCG):
    def __init__(self, restriction=''):
        CellGeometryCG.__init__(self, 'quadrilateral', 2, 2, restriction)
