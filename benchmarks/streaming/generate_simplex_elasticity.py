"""Generate simplex elasticity kernels and independent quadrature tensors."""

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

root = Path("/tmp/simplex-elasticity-reference")
shutil.copy(Path(ffcx.__file__).parent / "codegeneration/ufcx.h", root / "ufcx.h")
manifest = {
    "ffcx": ffcx.__version__,
    "basix": basix.__version__,
    "ufl": ufl.__version__,
    "mu": 1.0,
    "lambda": 1.25,
    "cases": {},
}
for cell, dim in [("triangle", 2), ("tetrahedron", 3)]:
    ce = basix.ufl.element("Lagrange", cell, 1, shape=(dim,))
    mesh = ufl.Mesh(ce)
    for degree in [1, 2, 3]:
        scalar = basix.ufl.element(
            "Lagrange", cell, degree, lagrange_variant=basix.LagrangeVariant.equispaced
        )
        e = basix.ufl.blocked_element(scalar, shape=(dim,))
        V = ufl.FunctionSpace(mesh, e)
        u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
        qdegree = 2 * (degree - 1)
        points, weights = basix.make_quadrature(getattr(basix.CellType, cell), qdegree)
        form = (
            2 * ufl.inner(ufl.sym(ufl.grad(u)), ufl.sym(ufl.grad(v)))
            + 1.25 * ufl.div(u) * ufl.div(v)
        ) * ufl.dx(
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
        derivatives = scalar.basix_element.tabulate(1, points)[1 : dim + 1, :, :, 0]
        matrices = []
        coords = []
        expected = []
        for k in range(32):
            B = (np.diag([1.1, 0.92, 1.2][:dim]) + 0.04 * np.ones((dim, dim))) * (1 + 0.001 * k)
            matrices.append(B)
            x = np.zeros((len(ce.basix_element.points), 3))
            x[:, :dim] = ce.basix_element.points @ B.T + 0.01 * k
            coords.append(x.ravel())
            g = np.einsum("dqi,dc->qic", derivatives, np.linalg.inv(B))
            w = weights * abs(np.linalg.det(B))
            tensor = (
                np.einsum("q,qic,qjc,ab->iajb", w, g, g, np.eye(dim))
                + np.einsum("q,qib,qja->iajb", w, g, g)
                + 1.25 * np.einsum("q,qia,qjb->iajb", w, g, g)
            )
            expected.append(tensor.reshape(e.dim, e.dim))
        np.savez(
            root / f"{case}.npz",
            points=points,
            weights=weights,
            nodes=scalar.basix_element.points,
            coordinate_nodes=ce.basix_element.points,
            matrices=matrices,
            ffcx_coordinates=coords,
            expected=expected,
        )
        manifest["cases"][case] = {
            "cell": cell,
            "dimension": dim,
            "degree": degree,
            "dofs": e.dim,
            "scalar_dofs": scalar.dim,
            "quadrature_degree": qdegree,
            "quadrature_points": len(points),
            "source_sha256": hashlib.sha256(sources[1].encode()).hexdigest(),
        }
        print(case, e.dim, len(points), flush=True)
(root / "manifest-ffcx.json").write_text(json.dumps(manifest, indent=2) + "\n")
