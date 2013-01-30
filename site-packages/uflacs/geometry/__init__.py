
from uflacs.geometry.cellcodegen          import CellGeometryNames
from uflacs.geometry.intervalcodegen      import IntervalGeometryCG
from uflacs.geometry.trianglecodegen      import TriangleGeometryCG
from uflacs.geometry.tetrahedroncodegen   import TetrahedronGeometryCG
from uflacs.geometry.quadrilateralcodegen import QuadrilateralGeometryCG
from uflacs.geometry.hexahedroncodegen    import HexahedronGeometryCG

def make_cellcg(cell): # FIXME: This is broken by design, need full cell with t and g dimensions
    if not isinstance(cell, str):
        cellname = cell.cellname()
    cgclass = {
        'interval':       IntervalGeometryCG,
        'triangle':       TriangleGeometryCG,
        'tetrahedron':    TetrahedronGeometryCG,
        'quadrilateral':  QuadrilateralGeometryCG,
        'hexahedron':     HexahedronGeometryCG,
        }[cellname]
    return cgclass()
