"""Test alternate float64 basis-table layouts without changing FFCx defaults.

Generate sources with ``forms.py`` or ``elasticity.py``, transform them here,
then compare both directories using the corresponding benchmark driver::

    python benchmarks/transpose_tables.py /tmp/before /tmp/transposed
    python benchmarks/forms.py compare --variant before=/tmp/before \
        --variant transposed=/tmp/transposed --output /tmp/layout-results.json

The ``small-linear-math`` policy restricts transposition to linear forms with
math calls and at most ten basis dofs. Both policies require at least eight
quadrature points. They are experimental policies, not compiler defaults.
"""

import argparse
import ast
import hashlib
import json
import re
import shutil
from pathlib import Path

import numpy as np

_TABLE = re.compile(
    r"static const double (FE\w+)\[(\d+)\]\[(\d+)\]\[(\d+)\]\[(\d+)\] = (.*?);",
    re.DOTALL,
)
_MATH = re.compile(
    r"\b(?:acos|acosh|asin|asinh|atan|atanh|cos|cosh|erf|exp|log|sin|sinh|tan|tanh)\("
)


def transpose(source: str, max_dofs: int | None = None) -> tuple[str, list[str]]:
    """Swap point/dof axes in declarations and every corresponding table read."""
    names = []
    for match in list(_TABLE.finditer(source)):
        name, permutations, entities, points, dofs, initializer = match.groups()
        npoints, ndofs = int(points), int(dofs)
        if int(permutations) != 1 or int(entities) != 1:
            continue
        if npoints < 8 or ndofs < 2 or (max_dofs is not None and ndofs > max_dofs):
            continue
        values = np.array(ast.literal_eval(initializer.replace("{", "[").replace("}", "]")))
        assert values.shape == (int(permutations), int(entities), npoints, ndofs)
        initializer = repr(values.swapaxes(2, 3).tolist()).replace("[", "{").replace("]", "}")
        declaration = (
            f"static const double {name}[{permutations}][{entities}][{dofs}][{points}]"
            f" = {initializer};"
        )
        placeholder = f"@TABLE_{name}@"
        source = source.replace(match.group(0), placeholder)
        access = re.compile(
            re.escape(name) + r"\[([^\[\]]+)\]\[([^\[\]]+)\]\[([^\[\]]+)\]\[([^\[\]]+)\]"
        )
        source = access.sub(lambda m: f"{name}[{m[1]}][{m[2]}][{m[4]}][{m[3]}]", source)
        source = source.replace(placeholder, declaration)
        names.append(name)
    return source, names


def main():
    """Copy generated benchmark inputs and change only their C table layout."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--policy", choices=["all", "small-linear-math"], default="all")
    args = parser.parse_args()
    if args.source.resolve() == args.destination.resolve():
        parser.error("Source and destination directories must differ")
    manifest = json.loads((args.source / "manifest.json").read_text())
    args.destination.mkdir(parents=True, exist_ok=True)
    changed = {}
    for path in args.source.iterdir():
        if not path.is_file() or path.suffix not in {".c", ".h", ".npz", ".npy"}:
            continue
        destination = args.destination / path.name
        if path.suffix != ".c":
            shutil.copyfile(path, destination)
            continue
        source = path.read_text()
        if "cases" in manifest:
            linear = len(manifest["cases"][path.stem]["shape"]) == 1
        else:
            linear = path.stem.startswith("residual")
        eligible = args.policy == "all" or (linear and _MATH.search(source))
        result, names = (
            transpose(source, 10 if args.policy == "small-linear-math" else None)
            if eligible
            else (source, [])
        )
        destination.write_text(result)
        changed[path.stem] = names
        digest = hashlib.sha256(result.encode()).hexdigest()
        if "cases" in manifest:
            manifest["cases"][path.stem]["source_sha256"] = digest
        else:
            manifest["source_sha256"][path.name] = digest
    manifest["table_layout_experiment"] = {"policy": args.policy, "transposed_tables": changed}
    (args.destination / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
