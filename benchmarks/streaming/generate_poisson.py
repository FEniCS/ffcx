"""Generate prescribed-quadrature Poisson kernels and independent tensors."""

# Generated C fragments retain their original source layout.

import hashlib
import json
import shutil
from pathlib import Path

import basix
import basix.ufl
import numpy as np
import ufl

import ffcx
from ffcx.compiler import compile_ufl_objects
from ffcx.options import get_options

root = Path("/tmp/poisson-vs-reference-current")
root.mkdir(exist_ok=True)
shutil.copy(Path(ffcx.__file__).parent / "codegeneration/ufcx.h", root / "ufcx.h")
manifest = {
    "ffcx": ffcx.__version__,
    "basix": basix.__version__,
    "ufl": ufl.__version__,
    "cases": {},
}
for cell in ["triangle", "quadrilateral", "tetrahedron", "hexahedron"]:
    dim = 2 if cell in ["triangle", "quadrilateral"] else 3
    coordinate = basix.ufl.element("Lagrange", cell, 1, shape=(dim,))
    mesh = ufl.Mesh(coordinate)
    for degree in [1, 2, 3]:
        e = basix.ufl.element(
            "Lagrange", cell, degree, lagrange_variant=basix.LagrangeVariant.equispaced
        )
        V = ufl.FunctionSpace(mesh, e)
        u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
        qdegree = 2 * (degree - 1) if cell in ["triangle", "tetrahedron"] else 2 * degree
        points, weights = basix.make_quadrature(getattr(basix.CellType, cell), qdegree)
        form = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx(
            metadata={
                "quadrature_rule": "custom",
                "quadrature_points": points,
                "quadrature_weights": weights,
            }
        )
        sources, suffixes = compile_ufl_objects(
            [form], options=get_options({"scalar_type": "float64"})
        )
        case = f"{cell}_p{degree}"
        for s, ext in zip(sources, suffixes):
            (root / f"ffcx_{case}{ext}").write_text(s)
        matrices = []
        coords = []
        expected = []
        derivatives = e.basix_element.tabulate(1, points)[1 : dim + 1, :, :, 0]
        for k in range(32):
            B = (np.diag([1.1, 0.92, 1.2][:dim]) + 0.04 * np.ones((dim, dim))) * (1 + 0.001 * k)
            matrices.append(B)
            x = np.zeros((len(coordinate.basix_element.points), 3))
            x[:, :dim] = coordinate.basix_element.points @ B.T + 0.01 * k
            coords.append(x.ravel())
            gradient = np.einsum("dqi,dc->qic", derivatives, np.linalg.inv(B))
            expected.append(
                np.einsum("q,qic,qjc->ij", weights * np.linalg.det(B), gradient, gradient)
            )
        np.savez(
            root / f"{case}.npz",
            points=points,
            weights=weights,
            nodes=e.basix_element.points,
            coordinate_nodes=coordinate.basix_element.points,
            matrices=matrices,
            ffcx_coordinates=coords,
            expected=expected,
        )
        manifest["cases"][case] = {
            "cell": cell,
            "dimension": dim,
            "degree": degree,
            "dofs": e.dim,
            "quadrature_degree": qdegree,
            "quadrature_points": len(points),
            "source_sha256": hashlib.sha256(sources[1].encode()).hexdigest(),
        }
        print(case, e.dim, len(points), flush=True)
(root / "manifest-ffcx.json").write_text(json.dumps(manifest, indent=2) + "\n")
