"""Generate and compare tetrahedral elasticity kernels from FFCx revisions.

Use each revision's Python environment for generation, sharing the quadrature
rules from the first directory. Then compare the generated sources::

    python benchmarks/elasticity.py generate /tmp/elasticity/change --revision change
    main-python benchmarks/elasticity.py generate /tmp/elasticity/main \
        --revision main --quadrature-from /tmp/elasticity/change
    release-python benchmarks/elasticity.py generate /tmp/elasticity/release \
        --revision v0.11.0 --quadrature-from /tmp/elasticity/change
    python benchmarks/elasticity.py compare --variant change=/tmp/elasticity/change \
        --variant main=/tmp/elasticity/main --variant v0.11.0=/tmp/elasticity/release \
        --cc gcc-16 clang-22 --cpu 2 --output results.json

Separate translation units prevent inlining into the timing driver. Every call
zeros the actual output size and cycles through 32 distorted cells. Full tensors
are checked across revisions after matching interpolation nodes. Timing samples
rotate the variant order to reduce ordering bias.
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


def forms(element, coordinate_element, degree, rules):
    """Build linear elasticity and neo-Hookean residual and tangent forms."""
    import ufl

    mesh = ufl.Mesh(coordinate_element)
    space = ufl.FunctionSpace(mesh, element)
    trial, test = ufl.TrialFunction(space), ufl.TestFunction(space)
    state = ufl.Coefficient(space)

    def measure(name, quadrature_degree):
        points, weights = rules[name]
        return ufl.dx(
            metadata={
                "quadrature_degree": quadrature_degree,
                "quadrature_rule": "custom",
                "quadrature_points": points,
                "quadrature_weights": weights,
            }
        )

    linear = (
        2 * ufl.inner(ufl.sym(ufl.grad(trial)), ufl.sym(ufl.grad(test)))
        + 1.25 * ufl.div(trial) * ufl.div(test)
    ) * measure("linear", 2 * (degree - 1))
    F = ufl.Identity(3) + ufl.grad(state)
    J = ufl.det(F)
    energy = ((ufl.inner(F, F) - 3) / 2 - ufl.ln(J) + 0.625 * ufl.ln(J) ** 2) * measure(
        "residual", 2 * degree
    )
    residual = ufl.derivative(energy, state, test)
    return {
        "linear": linear,
        "residual": residual,
        "tangent": ufl.derivative(residual, state, trial),
    }


def generate(args):
    """Save source, header, interpolation nodes, quadrature and stack versions."""
    import basix.ufl
    import ufl

    import ffcx
    from ffcx.compiler import compile_ufl_objects
    from ffcx.options import get_options

    root = args.directory
    root.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(Path(ffcx.__file__).parent / "codegeneration" / "ufcx.h", root / "ufcx.h")
    for degree in args.degrees:
        element = basix.ufl.element(
            "Lagrange",
            "tetrahedron",
            degree,
            shape=(3,),
            lagrange_variant=basix.LagrangeVariant.equispaced,
        )
        coordinate_element = basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,))
        np.save(root / f"nodes{degree}.npy", element.basix_element.points)
        rules = {}
        for name in ["linear", "residual"]:
            quadrature_degree = 2 * (degree - 1) if name == "linear" else 2 * degree
            if args.quadrature_from is not None:
                data = np.load(args.quadrature_from / f"quadrature_{name}{degree}.npz")
                points, weights = data["points"], data["weights"]
            else:
                points, weights = basix.make_quadrature(
                    basix.CellType.tetrahedron, quadrature_degree
                )
            np.savez(root / f"quadrature_{name}{degree}.npz", points=points, weights=weights)
            rules[name] = points, weights
        for name, form in forms(element, coordinate_element, degree, rules).items():
            sources, suffixes = compile_ufl_objects([form], options=get_options())
            for source, suffix in zip(sources, suffixes):
                (root / f"{name}{degree}{suffix}").write_text(source)
    manifest = {
        "revision": args.revision,
        "ffcx": ffcx.__version__,
        "basix": basix.__version__,
        "ufl": ufl.__version__,
        "numpy": np.__version__,
        "python": platform.python_version(),
        "source_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.glob("*.c"))
        },
        "header_sha256": hashlib.sha256((root / "ufcx.h").read_bytes()).hexdigest(),
    }
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


def load_kernel(root, case, size, ndofs, compiler, flags):
    """Compile a kernel with its own header and a separate C timing wrapper."""
    source = root / f"{case}.c"
    function = re.search(r"void (tabulate_tensor_\w+)\(", source.read_text())[1]
    build = root / ".build" / Path(compiler).name
    build.mkdir(parents=True, exist_ok=True)
    driver = build / f"{case}_driver.c"
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
    const double *w=ws+(k%32)*{ndofs}, *x=xs+(k%32)*12;
    memset(A,0,sizeof(A));
    {function}(A,w,NULL,x,NULL,NULL,NULL);
    check+=A[k%{size}];
  }}
  clock_gettime(CLOCK_MONOTONIC,&t1);
  checksum=check;
  return (t1.tv_sec-t0.tv_sec)+1e-9*(t1.tv_nsec-t0.tv_nsec);
}}
""")
    library = build / f"{case}.so"
    command = [
        compiler,
        *shlex.split(flags),
        "-shared",
        "-fPIC",
        f"-I{root.resolve()}",
        str(source),
        str(driver),
        "-o",
        str(library),
        "-lm",
    ]
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode:
        raise RuntimeError(f"Compilation failed: {command}\n{process.stderr}")
    module = ctypes.CDLL(str(library.resolve()))
    pointer = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    module.evaluate.argtypes = [pointer, pointer, pointer]
    module.evaluate.restype = None
    module.measure.argtypes = [ctypes.c_long, pointer, pointer]
    module.measure.restype = ctypes.c_double
    return module


def inputs(nodes):
    """Create finite, orientation-preserving cells and displacement coefficients."""
    coefficients, coordinates = [], []
    for k in range(32):
        B = np.array([[1.1, 0.1, 0.04], [0.03, 0.92, 0.08], [0.02, 0.05, 1.2]]) * (1 + 0.001 * k)
        coordinates.append(np.vstack([np.zeros(3), B.T]).ravel())
        physical = nodes @ B.T
        coefficients.append(
            (0.03 * np.sin(physical + np.arange(3) * 0.4) + 0.005 * physical**2).ravel()
        )
    return np.array(coefficients), np.array(coordinates)


def compare(args):
    """Verify complete tensors and collect interleaved timings across revisions."""
    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    variants = dict(item.split("=", 1) for item in args.variant)
    roots = {label: Path(path) for label, path in variants.items()}
    labels = list(roots)
    manifest = {
        label: json.loads((root / "manifest.json").read_text()) for label, root in roots.items()
    }
    model = next(
        (
            line.split(":", 1)[1].strip()
            for line in Path("/proc/cpuinfo").read_text().splitlines()
            if line.startswith("model name")
        ),
        platform.processor(),
    )
    report = {
        "environment": {
            "platform": platform.platform(),
            "cpu_model": model,
            "cpu_affinity": args.cpu,
            "flags": args.flags,
            "duration_s": args.duration,
            "repetitions": args.repetitions,
        },
        "scope": {
            "dimension": 3,
            "cell": "tetrahedron",
            "geometry_degree": 1,
            "degrees": args.degrees,
            "scalar_type": "float64",
            "cell_inputs": 32,
            "mu": 1.0,
            "lambda": 1.25,
            "zero_output_per_call": True,
            "hyperelastic_model": "compressible neo-Hookean",
        },
        "variants": manifest,
        "compilers": {},
        "results": {},
    }
    for compiler in args.cc:
        report["compilers"][compiler] = subprocess.check_output(
            [compiler, "--version"], text=True
        ).splitlines()[0]
        results = report["results"][compiler] = {}
        for degree in args.degrees:
            target_nodes = np.load(roots[labels[0]] / f"nodes{degree}.npy")
            ndofs = 3 * len(target_nodes)
            for name in ["linear", "residual", "tangent"]:
                case = f"{name}{degree}"
                size = ndofs if name == "residual" else ndofs**2
                modules, data, outputs, counts = {}, {}, {}, {}
                for label, root in roots.items():
                    nodes = np.load(root / f"nodes{degree}.npy")
                    distances = np.linalg.norm(target_nodes[:, None] - nodes[None, :], axis=2)
                    permutation = np.argmin(distances, axis=1)
                    assert np.max(np.min(distances, axis=1)) < 1e-12
                    assert len(np.unique(permutation)) == len(nodes)
                    permutation = (3 * permutation[:, None] + np.arange(3)).ravel()
                    module = modules[label] = load_kernel(
                        root, case, size, ndofs, compiler, args.flags
                    )
                    ws, xs = data[label] = inputs(nodes)
                    tensors = []
                    for w, x in zip(ws, xs):
                        tensor = np.zeros(size)
                        module.evaluate(tensor, w, x)
                        assert np.isfinite(tensor).all()
                        tensor = (
                            tensor[permutation]
                            if name == "residual"
                            else tensor.reshape(ndofs, ndofs)[np.ix_(permutation, permutation)]
                        )
                        tensors.append(tensor.ravel())
                    outputs[label] = np.array(tensors)
                    for rule in ["linear", "residual"]:
                        common = np.load(roots[labels[0]] / f"quadrature_{rule}{degree}.npz")
                        current = np.load(root / f"quadrature_{rule}{degree}.npz")
                        for key in ["points", "weights"]:
                            np.testing.assert_array_equal(common[key], current[key])
                    counts[label] = 1000
                    while module.measure(counts[label], ws, xs) < args.duration:
                        counts[label] *= 2
                errors = {}
                for label in labels:
                    np.testing.assert_allclose(
                        outputs[label], outputs[labels[0]], atol=2e-11, rtol=2e-10
                    )
                    errors[label] = float(np.max(np.abs(outputs[label] - outputs[labels[0]])))
                samples = {label: [] for label in labels}
                for repetition in range(args.repetitions):
                    order = labels[repetition % len(labels) :] + labels[: repetition % len(labels)]
                    for label in order:
                        ws, xs = data[label]
                        samples[label].append(
                            modules[label].measure(counts[label], ws, xs) / counts[label] * 1e9
                        )
                        assert np.isfinite(ctypes.c_double.in_dll(modules[label], "checksum").value)
                results[case] = {
                    label: {
                        "best_ns": min(samples[label]),
                        "median_ns": statistics.median(samples[label]),
                        "samples_ns": samples[label],
                        "calls": counts[label],
                        "max_absolute_tensor_error": errors[label],
                    }
                    for label in labels
                }
                print(
                    compiler,
                    case,
                    {label: round(results[case][label]["best_ns"], 2) for label in labels},
                    flush=True,
                )
                args.output.write_text(json.dumps(report, indent=2) + "\n")


def main():
    """Parse generation or comparison options."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    gen = commands.add_parser("generate")
    gen.add_argument("directory", type=Path)
    gen.add_argument("--revision", required=True)
    gen.add_argument("--quadrature-from", type=Path)
    timing = commands.add_parser("compare")
    timing.add_argument("--variant", action="append", required=True, metavar="LABEL=DIRECTORY")
    timing.add_argument("--output", type=Path, required=True)
    timing.add_argument("--cc", nargs="+", default=["gcc-16", "clang-22"])
    timing.add_argument("--flags", default="-O3 -march=native -ffp-contract=fast")
    timing.add_argument("--cpu", type=int)
    timing.add_argument("--duration", type=float, default=0.1)
    timing.add_argument("--repetitions", type=int, default=9)
    for command in [gen, timing]:
        command.add_argument("--degrees", type=int, nargs="+", default=[1, 2, 3])
    args = parser.parse_args()
    generate(args) if args.command == "generate" else compare(args)


if __name__ == "__main__":
    main()
