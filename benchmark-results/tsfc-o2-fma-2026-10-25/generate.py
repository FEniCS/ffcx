"""Generate equivalent standalone cell kernels and nodal-order metadata."""

import argparse
import itertools
import json
import re
from pathlib import Path

import numpy as np
import ufl

FORMS = ("mass", "laplace", "vector_laplace", "elasticity", "mass_action", "laplace_action")


def physical_points(points):
    """Map reference nodes to the shared non-affine physical cell."""
    x, y = points[:, 0], points[:, 1]
    if points.shape[1] == 2:
        return np.column_stack((x + 0.2 * x * y + 0.1 * y, y + 0.1 * x + 0.4 * x * y))
    z = points[:, 2]
    return np.column_stack((x + 0.1 * y * z, y + 0.05 * x, z + 0.1 * x * y))


def coefficient_values(points):
    """Sample a nonconstant coefficient in each compiler's nodal ordering."""
    return 0.3 + points @ np.arange(1, points.shape[1] + 1) * 0.17 + 0.2 * points[:, 0] ** 2


def basix_element(cellname, degree):
    """Return a tensor-product Basix element and its basis-ordered nodes."""
    import basix
    import basix.ufl

    raw = basix.create_tp_element(
        basix.ElementFamily.P,
        getattr(basix.CellType, cellname),
        degree,
        basix.LagrangeVariant.gll_warped,
    )
    values = raw.tabulate(0, raw.points)[0, :, :, 0]
    indices = np.argmax(np.abs(values), axis=0)
    np.testing.assert_allclose(values[indices, np.arange(raw.dim)], 1, atol=1e-12)
    return basix.ufl.wrap_element(raw), raw.points[indices]


def finat_element(cellname, degree):
    """Return a development FInAT element, requiring all coordinate factors."""
    import finat.ufl as fu
    from finat.element_factory import create_element

    element = fu.FiniteElement("Q", getattr(ufl, cellname), degree, variant="spectral")
    converted = create_element(element)
    dim = getattr(ufl, cellname).topological_dimension
    assert converted.index_shape == (degree + 1,) * dim
    assert len(converted.product.factors) == dim
    nodes = [
        [tuple(next(iter(node.get_point_dict()))) for node in factor.fiat_equivalent.dual.nodes]
        for factor in converted.product.factors
    ]
    points = np.array([sum(p, ()) for p in itertools.product(*nodes)])
    return element, points


def generate(compiler, out, degrees):
    """Emit six equivalent forms on each cell at every requested degree."""
    out.mkdir(parents=True, exist_ok=True)
    for cellname, dim in (("quadrilateral", 2), ("hexahedron", 3)):
        for degree in degrees:
            make = basix_element if compiler == "ffcx" else finat_element
            scalar, points = make(cellname, degree)
            coordinate_scalar, coordinate_points = make(cellname, 1)
            if compiler == "ffcx":
                import basix.ufl

                vector = basix.ufl.blocked_element(scalar, shape=(dim,))
                coordinate = basix.ufl.blocked_element(coordinate_scalar, shape=(dim,))
            else:
                import finat.ufl as fu

                vector = fu.VectorElement(scalar, dim=dim)
                coordinate = fu.VectorElement(coordinate_scalar, dim=dim)
            mesh = ufl.Mesh(coordinate)
            V, W = ufl.FunctionSpace(mesh, scalar), ufl.FunctionSpace(mesh, vector)
            u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
            uw, vw = ufl.TrialFunction(W), ufl.TestFunction(W)
            w = ufl.Coefficient(V)
            dx = ufl.dx(metadata={"quadrature_degree": 2 * degree})
            forms = {
                "mass": u * v * dx,
                "laplace": ufl.inner(ufl.grad(u), ufl.grad(v)) * dx,
                "vector_laplace": ufl.inner(ufl.grad(uw), ufl.grad(vw)) * dx,
                "elasticity": (
                    ufl.inner(ufl.sym(ufl.grad(uw)), ufl.sym(ufl.grad(vw)))
                    + ufl.div(uw) * ufl.div(vw)
                )
                * dx,
                "mass_action": w * v * dx,
                "laplace_action": ufl.inner(ufl.grad(w), ufl.grad(v)) * dx,
            }
            coords = physical_points(coordinate_points)
            if compiler == "ffcx":
                coords = np.pad(coords, ((0, 0), (0, 3 - dim)))
            for name, form in forms.items():
                tag = f"{cellname}_p{degree}_{name}"
                if compiler == "ffcx":
                    import ffcx.compiler
                    import ffcx.options

                    code, _ = ffcx.compiler.compile_ufl_objects(
                        [form], ffcx.options.get_options({"sum_factorization": True})
                    )
                    text = "\n".join(code)
                    function = re.search(r"\.tabulate_tensor_float64 = (\w+)", text)[1]
                    # Actions must evaluate the coefficient by directional contractions.
                    if name.endswith("action"):
                        assert "_c0[" in text
                else:
                    import loopy
                    import tsfc

                    (kernel,) = tsfc.compile_form(form, prefix="k", parameters={"mode": "spectral"})
                    text = loopy.generate_code_v2(kernel.ast).device_code()
                    function = kernel.name
                (out / f"{tag}.c").write_text(text)
                components = dim if name in ("vector_laplace", "elasticity") else 1
                ndofs = len(points) * components
                metadata = {
                    "compiler": compiler,
                    "cell": cellname,
                    "dim": dim,
                    "degree": degree,
                    "form": name,
                    "function": function,
                    "points": points.tolist(),
                    "coords": coords.ravel().tolist(),
                    "coefficients": coefficient_values(points).tolist(),
                    "components": components,
                    "ndofs": ndofs,
                    "rank": len(form.arguments()),
                    "tensor_size": ndofs ** len(form.arguments()),
                    "tensor_index_shape": [degree + 1] * dim,
                }
                (out / f"{tag}.json").write_text(json.dumps(metadata, indent=2))
                print(f"{compiler}: {tag}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("compiler", choices=("ffcx", "tsfc"))
    parser.add_argument("out", type=Path)
    parser.add_argument("--degrees", type=int, nargs="+", default=list(range(1, 9)))
    args = parser.parse_args()
    generate(args.compiler, args.out, args.degrees)
