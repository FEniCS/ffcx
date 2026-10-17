"""Benchmark linear and bilinear cell forms across FFCx revisions.

Generate each revision into its own directory, then compare, for example::

    python benchmarks/forms.py generate /tmp/forms/before --revision before
    python benchmarks/forms.py generate /tmp/forms/after --revision after
    python benchmarks/forms.py compare --variant before=/tmp/forms/before \
        --variant after=/tmp/forms/after --output /tmp/forms/results.json

The suite covers scalar loads, mass, diffusion, transport, vector mass, mixed
divergence, and H(div)/H(curl) operators on simplex and tensor-product cells.
Kernels and timing drivers are separate translation units. All output entries
are checked on 32 cells before collecting interleaved samples, including output
zeroing. Each revision uses its default quadrature rules. Cross-release results
therefore include changes to quadrature selection in the dependency stack.
These are local tensor assembly timings, not global assembly timings.
"""

import argparse
import ctypes
import hashlib
import json
import os
import platform
import re
import shlex
import shutil
import statistics
import subprocess
from pathlib import Path

import numpy as np


def cases(cells, degrees):
    """Yield forms together with deterministic coefficient and coordinate data."""
    import basix
    import basix.ufl
    import ufl

    for cell in cells:
        dim = len(basix.topology(getattr(basix.CellType, cell))) - 1
        coordinate = basix.ufl.element("Lagrange", cell, 1, shape=(dim,))
        mesh = ufl.Mesh(coordinate)
        vertices = coordinate.basix_element.points
        B = np.eye(dim) + 0.08 * np.ones((dim, dim))
        coordinates = []
        for k in range(32):
            physical = vertices @ (B * (1 + 0.001 * k)).T + 0.01 * k
            padded = np.zeros((len(vertices), 3))
            padded[:, :dim] = physical
            coordinates.append(padded.ravel())
        coordinates = np.array(coordinates)
        for degree in degrees:
            if cell == "hexahedron" and degree > 2:
                continue
            scalar = basix.ufl.element(
                "Lagrange", cell, degree, lagrange_variant=basix.LagrangeVariant.equispaced
            )
            V = ufl.FunctionSpace(mesh, scalar)
            u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
            f = ufl.Coefficient(V)
            nodes = scalar.basix_element.points
            coefficients = np.array(
                [
                    1 + 0.2 * np.sin((nodes @ (B * (1 + 0.001 * k)).T).sum(axis=1) + 0.1 * k)
                    for k in range(32)
                ]
            )
            direction = ufl.as_vector([1 + 0.3 * i for i in range(dim)])
            tensor = ufl.as_matrix(
                [[2.0 if i == j else 0.2 for j in range(dim)] for i in range(dim)]
            )
            forms = {
                "load": v * ufl.dx,
                "weighted_load": f * v * ufl.dx,
                "gradient_load": ufl.inner(direction, ufl.grad(v)) * ufl.dx,
                "mass": ufl.inner(u, v) * ufl.dx,
                "weighted_mass": f * ufl.inner(u, v) * ufl.dx,
                "diffusion": ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx,
                "anisotropic_diffusion": ufl.inner(tensor * ufl.grad(u), ufl.grad(v)) * ufl.dx,
                "weighted_diffusion": f * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx,
                "transport": f * ufl.dot(direction, ufl.grad(u)) * v * ufl.dx,
                "reaction_diffusion": (f * u * v + ufl.inner(ufl.grad(u), ufl.grad(v))) * ufl.dx,
            }
            vector = basix.ufl.blocked_element(scalar, shape=(dim,))
            W = ufl.FunctionSpace(mesh, vector)
            a, b = ufl.TrialFunction(W), ufl.TestFunction(W)
            forms["vector_mass"] = ufl.inner(a, b) * ufl.dx
            Q = ufl.FunctionSpace(mesh, basix.ufl.element("DG", cell, degree - 1))
            q = ufl.TestFunction(Q)
            forms["mixed_divergence"] = ufl.div(a) * q * ufl.dx
            if cell in ("triangle", "tetrahedron") and degree <= 2:
                for family, operator, label in [
                    ("RT", ufl.div, "hdiv"),
                    ("N1curl", ufl.curl, "hcurl"),
                ]:
                    E = basix.ufl.element(family, cell, degree)
                    X = ufl.FunctionSpace(mesh, E)
                    p, r = ufl.TrialFunction(X), ufl.TestFunction(X)
                    forms[label] = (ufl.inner(p, r) + ufl.inner(operator(p), operator(r))) * ufl.dx
            for name, form in forms.items():
                shape = [arg.ufl_function_space().ufl_element().dim for arg in form.arguments()]
                yield f"{cell}_p{degree}_{name}", form, shape, coefficients, coordinates


def generate(args):
    """Write sources and input arrays using the active FFCx revision."""
    import basix
    import ufl

    import ffcx
    from ffcx.compiler import compile_ufl_objects
    from ffcx.options import get_options

    root = args.directory
    root.mkdir(parents=True, exist_ok=True)
    shutil.copy(Path(ffcx.__file__).parent / "codegeneration/ufcx.h", root / "ufcx.h")
    manifest = {
        "revision": args.revision,
        "ffcx": ffcx.__version__,
        "basix": basix.__version__,
        "ufl": ufl.__version__,
        "numpy": np.__version__,
        "python": platform.python_version(),
        "cases": {},
    }
    for name, form, shape, coefficients, coordinates in cases(args.cells, args.degrees):
        sources, suffixes = compile_ufl_objects(
            [form], options=get_options({"scalar_type": "float64"})
        )
        for source, suffix in zip(sources, suffixes):
            (root / f"{name}{suffix}").write_text(source)
        np.savez(root / f"{name}.npz", coefficients=coefficients, coordinates=coordinates)
        manifest["cases"][name] = {
            "shape": shape,
            "source_sha256": hashlib.sha256(sources[1].encode()).hexdigest(),
            "inputs_sha256": hashlib.sha256((root / f"{name}.npz").read_bytes()).hexdigest(),
        }
        print(name, shape, flush=True)
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def load_kernel(root, name, shape, data, cc, flags):
    """Compile a float64 kernel and a separate timing driver."""
    source = root / f"{name}.c"
    function = re.search(r"void (tabulate_tensor_\w+)\(", source.read_text())[1]
    build = root / ".build" / Path(cc).name
    build.mkdir(parents=True, exist_ok=True)
    size = int(np.prod(shape))
    nw, nx = data["coefficients"].shape[1], data["coordinates"].shape[1]
    driver = build / f"{name}_driver.c"
    driver.write_text(f"""#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <time.h>
void {function}(double*,const double*,const double*,const double*,const int*,const uint8_t*,void*);
double checksum;
void evaluate(double *A,const double *w,const double *x) {{
  {function}(A,w,NULL,x,NULL,NULL,NULL);
}}
double measure(long count,const double *ws,const double *xs) {{
  double A[{size}],check=0;
  struct timespec t0,t1;
  clock_gettime(CLOCK_MONOTONIC,&t0);
  for(long k=0;k<count;++k) {{
    const double *w=ws+(k%32)*{nw}, *x=xs+(k%32)*{nx};
    memset(A,0,sizeof(A));
    {function}(A,w,NULL,x,NULL,NULL,NULL);
    check+=A[k%{size}];
  }}
  clock_gettime(CLOCK_MONOTONIC,&t1);
  checksum=check;
  return (t1.tv_sec-t0.tv_sec)+1e-9*(t1.tv_nsec-t0.tv_nsec);
}}
""")
    library = build / f"{name}.so"
    process = subprocess.run(
        [
            cc,
            *shlex.split(flags),
            "-shared",
            "-fPIC",
            f"-I{root.resolve()}",
            str(source),
            str(driver),
            "-o",
            str(library),
            "-lm",
        ],
        capture_output=True,
        text=True,
    )
    if process.returncode:
        raise RuntimeError(process.stderr)
    module = ctypes.CDLL(str(library.resolve()))
    pointer = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    module.evaluate.argtypes = [pointer, pointer, pointer]
    module.evaluate.restype = None
    module.measure.argtypes = [ctypes.c_long, pointer, pointer]
    module.measure.restype = ctypes.c_double
    return module


def compare(args):
    """Compare complete tensors and time variants in rotating order."""
    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    roots = {label: Path(path) for label, path in (v.split("=", 1) for v in args.variant)}
    labels = list(roots)
    manifests = {
        label: json.loads((root / "manifest.json").read_text()) for label, root in roots.items()
    }
    names = list(manifests[labels[0]]["cases"])
    for manifest in manifests.values():
        assert list(manifest["cases"]) == names
    report = {
        "flags": args.flags,
        "cpu": args.cpu,
        "platform": platform.platform(),
        "repetitions": args.repetitions,
        "duration_s": args.duration,
        "cell_inputs": 32,
        "zero_output_per_call": True,
        "scalar_type": "float64",
        "cpu_model": next(
            (
                line.split(":", 1)[1].strip()
                for line in Path("/proc/cpuinfo").read_text().splitlines()
                if line.startswith("model name")
            ),
            platform.processor(),
        ),
        "variants": manifests,
        "compilers": {},
        "results": {},
    }
    for cc in args.cc:
        report["compilers"][cc] = subprocess.check_output(
            [cc, "--version"], text=True
        ).splitlines()[0]
        results = report["results"][cc] = {}
        for name in names:
            if args.filter and not re.search(args.filter, name):
                continue
            modules, data, counts, outputs = {}, {}, {}, {}
            shape = manifests[labels[0]]["cases"][name]["shape"]
            for label, root in roots.items():
                assert manifests[label]["cases"][name]["shape"] == shape
                data[label] = np.load(root / f"{name}.npz")
                ws, xs = data[label]["coefficients"], data[label]["coordinates"]
                assert len(ws) == len(xs) == 32
                np.testing.assert_array_equal(ws, data[labels[0]]["coefficients"])
                np.testing.assert_array_equal(xs, data[labels[0]]["coordinates"])
                m = modules[label] = load_kernel(root, name, shape, data[label], cc, args.flags)
                tensors = []
                for w, x in zip(ws, xs):
                    tensor = np.zeros(int(np.prod(shape)))
                    m.evaluate(tensor, w, x)
                    assert np.isfinite(tensor).all()
                    tensors.append(tensor)
                outputs[label] = np.array(tensors)
                np.testing.assert_allclose(
                    outputs[label], outputs[labels[0]], rtol=2e-10, atol=2e-11
                )
                count = 1000
                while m.measure(count, ws, xs) < args.duration:
                    count *= 2
                counts[label] = count
            samples = {label: [] for label in labels}
            for repetition in range(args.repetitions):
                order = labels[repetition % len(labels) :] + labels[: repetition % len(labels)]
                for label in order:
                    ws, xs = data[label]["coefficients"], data[label]["coordinates"]
                    samples[label].append(
                        modules[label].measure(counts[label], ws, xs) / counts[label] * 1e9
                    )
                    assert np.isfinite(ctypes.c_double.in_dll(modules[label], "checksum").value)
            results[name] = {
                label: {
                    "best_ns": min(samples[label]),
                    "median_ns": statistics.median(samples[label]),
                    "samples_ns": samples[label],
                    "calls": counts[label],
                    "max_absolute_tensor_error": float(
                        np.max(np.abs(outputs[label] - outputs[labels[0]]))
                    ),
                }
                for label in labels
            }
            print(
                cc,
                name,
                {label: round(results[name][label]["best_ns"], 2) for label in labels},
                flush=True,
            )
            args.output.write_text(json.dumps(report, indent=2) + "\n")


def main():
    """Generate the suite or compare generated kernels."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    gen = commands.add_parser("generate")
    gen.add_argument("directory", type=Path)
    gen.add_argument("--revision", required=True)
    gen.add_argument(
        "--cells", nargs="+", default=["triangle", "tetrahedron", "quadrilateral", "hexahedron"]
    )
    gen.add_argument("--degrees", type=int, nargs="+", default=[1, 2, 3])
    timing = commands.add_parser("compare")
    timing.add_argument("--variant", action="append", required=True)
    timing.add_argument("--output", type=Path, required=True)
    timing.add_argument("--cc", nargs="+", default=["gcc-16", "clang-22"])
    timing.add_argument("--flags", default="-O3 -march=native -ffp-contract=fast")
    timing.add_argument("--cpu", type=int, default=2)
    timing.add_argument("--duration", type=float, default=0.05)
    timing.add_argument("--repetitions", type=int, default=7)
    timing.add_argument("--filter")
    args = parser.parse_args()
    generate(args) if args.command == "generate" else compare(args)


if __name__ == "__main__":
    main()
