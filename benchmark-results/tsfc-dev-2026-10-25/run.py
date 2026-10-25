"""Build, validate and time the generated kernels with one shared C driver."""

import argparse
import concurrent.futures
import csv
import json
import random
import resource
import statistics
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FLAGS = {
    "O2": ["-O2", "-march=native"],
    "O3fast": ["-O3", "-march=native", "-ffast-math"],
}


def command(args, **kwargs):
    """Run a checked subprocess and capture its output."""
    return subprocess.run(args, check=True, capture_output=True, text=True, **kwargs)


def build_one(root, flagset, compiler, tag, cc, reuse):
    """Compile one kernel with the common driver and save its full tensor."""
    source = root / compiler / f"{tag}.c"
    meta = json.loads(source.with_suffix(".json").read_text())
    out = root / "build" / flagset / compiler / tag
    out.mkdir(parents=True, exist_ok=True)
    declarations = []
    for name, values in (
        ("coords", meta["coords"]),
        ("coefficients", meta["coefficients"]),
        ("constants", [0]),
    ):
        declarations.append(
            f"static const double {name}[] = {{"
            + ",".join(format(v, ".17g") for v in values)
            + "};"
        )
    (out / "inputs.h").write_text("\n".join(declarations))
    args = [
        cc,
        *FLAGS[flagset],
        "-std=c11",
        "-D_POSIX_C_SOURCE=200809L",
        f"-DKERNEL={meta['function']}",
        f"-DTSFC={int(compiler == 'tsfc')}",
        f"-DACTION={int(meta['rank'] == 1)}",
        f"-DTENSOR_SIZE={meta['tensor_size']}",
        "-I",
        str(out),
        "-I",
        str(REPO / "ffcx/codegeneration"),
        str(source),
        str(HERE / "driver.c"),
        "-o",
        str(out / "bench"),
        "-lm",
    ]
    saved = out / "command.json"
    if not (
        reuse
        and (out / "bench").exists()
        and saved.exists()
        and json.loads(saved.read_text()) == args
    ):
        proc = subprocess.run(args, capture_output=True, text=True)
        (out / "compile.log").write_text(proc.stdout + proc.stderr)
        saved.write_text(json.dumps(args))
        if proc.returncode:
            raise RuntimeError(f"Build failed: {flagset}/{compiler}/{tag}: {proc.stderr}")
    command([str(out / "bench"), str(out / "tensor.bin")])
    return flagset, compiler, tag


def validate(root, tags):
    """Compare complete tensors after matching scalar and vector nodal orders."""
    rows = []
    for flagset in FLAGS:
        for tag in tags:
            fmeta, tmeta = [
                json.loads((root / compiler / f"{tag}.json").read_text())
                for compiler in ("ffcx", "tsfc")
            ]
            fp, tp = np.array(fmeta["points"]), np.array(tmeta["points"])
            distances = np.linalg.norm(fp[:, None, :] - tp[None, :, :], axis=2)
            permutation = np.argmin(distances, axis=1)
            assert len(set(permutation)) == len(permutation)
            assert np.max(distances[np.arange(len(fp)), permutation]) < 1e-12
            components = fmeta["components"]
            permutation = (permutation[:, None] * components + np.arange(components)).ravel()
            f, t = [
                np.fromfile(root / "build" / flagset / compiler / tag / "tensor.bin")
                for compiler in ("ffcx", "tsfc")
            ]
            if fmeta["rank"] == 2:
                t = t.reshape(fmeta["ndofs"], fmeta["ndofs"])[
                    np.ix_(permutation, permutation)
                ].ravel()
            else:
                t = t[permutation]
            scale = np.max(np.abs(f))
            error = np.max(np.abs(f - t)) / scale
            relative_l2 = np.linalg.norm(f - t) / np.linalg.norm(f)
            np.testing.assert_allclose(
                f, t, rtol=2e-10, atol=1e-11 * scale, err_msg=f"{flagset}/{tag}"
            )
            rows.append(
                {
                    "flagset": flagset,
                    "case": tag,
                    "relative_max_error": float(error),
                    "relative_l2_error": float(relative_l2),
                }
            )
    return rows


def summary(records, root, output):
    """Write medians of the best samples from independent timing rounds."""
    groups = {}
    for record in records:
        key = record["flagset"], record["case"], record["compiler"]
        groups.setdefault(key, []).append(min(record["ns"]))
    rows = []
    for flagset, tag in sorted({(flagset, tag) for flagset, tag, _ in groups}):
        if any((flagset, tag, compiler) not in groups for compiler in ("ffcx", "tsfc")):
            continue
        meta = json.loads((root / "ffcx" / f"{tag}.json").read_text())
        f, t = [statistics.median(groups[flagset, tag, compiler]) for compiler in ("ffcx", "tsfc")]
        rows.append(
            {
                "flags": flagset,
                "cell": meta["cell"],
                "degree": meta["degree"],
                "form": meta["form"],
                "ffcx_us": f / 1000,
                "tsfc_us": t / 1000,
                "ffcx_over_tsfc": f / t,
                "ffcx_best_ns_by_round": json.dumps(groups[flagset, tag, "ffcx"]),
                "tsfc_best_ns_by_round": json.dumps(groups[flagset, tag, "tsfc"]),
            }
        )
    if rows:
        with (output / "results.csv").open("w") as file:
            writer = csv.DictWriter(file, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def main():
    """Build and validate all cases before timing them in shuffled order."""
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--cc", default="gcc-14")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--cpu", type=int, default=2)
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--reuse-build", action="store_true")
    args = parser.parse_args()
    # TSFC's high-degree elasticity kernels need more than the default 8 MiB.
    # Apply the same process stack limit to both compilers' executables.
    _, hard_limit = resource.getrlimit(resource.RLIMIT_STACK)
    resource.setrlimit(resource.RLIMIT_STACK, (64 * 1024 * 1024, hard_limit))
    root = args.root.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    tags = sorted(p.stem for p in (root / "ffcx").glob("*.json"))
    assert tags == sorted(p.stem for p in (root / "tsfc").glob("*.json"))
    jobs = [
        (flags, compiler, tag) for flags in FLAGS for tag in tags for compiler in ("ffcx", "tsfc")
    ]
    print(f"Building {len(jobs)} executables", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(build_one, root, *job, args.cc, args.reuse_build) for job in jobs]
        for index, future in enumerate(concurrent.futures.as_completed(futures), 1):
            result = future.result()
            if index % 16 == 0:
                print(f"Built {index}/{len(jobs)}: {result}", flush=True)
    validation = validate(root, tags)
    (args.output / "validation.json").write_text(json.dumps(validation, indent=2))
    print(
        f"Validated {len(validation)} tensor comparisons; max relative error "
        f"{max(r['relative_max_error'] for r in validation):.3g}",
        flush=True,
    )
    if args.validate_only:
        return
    records = []
    for sweep in range(args.rounds):
        random.Random(100 + sweep).shuffle(jobs)
        for index, (flagset, compiler, tag) in enumerate(jobs, 1):
            executable = root / "build" / flagset / compiler / tag / "bench"
            result = command(["taskset", "-c", str(args.cpu), str(executable)], timeout=120)
            record = json.loads(result.stdout)
            record.update(
                {"round": sweep + 1, "flagset": flagset, "compiler": compiler, "case": tag}
            )
            records.append(record)
            (args.output / "samples.json").write_text(json.dumps(records, indent=2))
            if index % 16 == 0:
                print(f"Timing round {sweep + 1}/{args.rounds}: {index}/{len(jobs)}", flush=True)
        summary(records, root, args.output)
    summary(records, root, args.output)


if __name__ == "__main__":
    main()
