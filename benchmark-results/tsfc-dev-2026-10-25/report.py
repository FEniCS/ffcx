# ruff: noqa: E501
"""Create tables and a figure from the recorded benchmark samples."""

import csv
import json
import math
import statistics
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FORMS = ("mass", "laplace", "vector_laplace", "elasticity", "mass_action", "laplace_action")
LABELS = ("Mass", "Laplace", "Vector Laplace", "Elasticity", "Mass action", "Laplace action")


def main():
    """Save a standalone comparison figure, full tables and aggregate statistics."""
    with (HERE / "results.csv").open() as file:
        rows = list(csv.DictReader(file))
    for row in rows:
        row["degree"] = int(row["degree"])
        for key in ("ffcx_us", "tsfc_us", "ffcx_over_tsfc"):
            row[key] = float(row[key])
    lookup = {(r["flags"], r["cell"], r["form"], r["degree"]): r for r in rows}
    assert len(lookup) == 192
    manifest = json.loads((HERE / "environment.json").read_text())
    validation = json.loads((HERE / "validation.json").read_text())
    audit = json.loads((HERE / "factorisation-audit.json").read_text())
    spreads = []
    for row in rows:
        for compiler in ("ffcx", "tsfc"):
            minima = json.loads(row[f"{compiler}_best_ns_by_round"])
            assert len(minima) == 3
            spreads.append(max(minima) / min(minima) - 1)
    fallback = {r["case"] for r in audit if not r["ffcx_directional_contraction"]}
    aggregates = []
    for cell in ("quadrilateral", "hexahedron"):
        for flags in ("O2", "O3fast"):
            subset = [r for r in rows if r["cell"] == cell and r["flags"] == flags]
            directional = [
                r for r in subset if f"{cell}_p{r['degree']}_{r['form']}" not in fallback
            ]
            ratios = [r["ffcx_over_tsfc"] for r in subset]
            aggregates.append(
                {
                    "cell": cell,
                    "flags": flags,
                    "cases": len(subset),
                    "ffcx_faster_cases": sum(r < 1 for r in ratios),
                    "geometric_mean_ffcx_over_tsfc": math.exp(
                        statistics.mean(map(math.log, ratios))
                    ),
                    "geometric_mean_directional_only": math.exp(
                        statistics.mean(math.log(r["ffcx_over_tsfc"]) for r in directional)
                    ),
                    "min_ffcx_over_tsfc": min(ratios),
                    "max_ffcx_over_tsfc": max(ratios),
                }
            )
    (HERE / "aggregates.json").write_text(json.dumps(aggregates, indent=2))
    fig, axes = plt.subplots(1, 2, figsize=(15, 8), constrained_layout=True)
    norm = colors.TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
    for ax, flags in zip(axes, ("O2", "O3fast")):
        values = np.array(
            [
                [lookup[flags, cell, form, p]["ffcx_over_tsfc"] for p in range(1, 9)]
                for cell in ("quadrilateral", "hexahedron")
                for form in FORMS
            ]
        )
        image = ax.imshow(np.log2(values), norm=norm, cmap="RdYlGn_r", aspect="auto")
        for i, (cell, form) in enumerate(
            (cell, form) for cell in ("quadrilateral", "hexahedron") for form in FORMS
        ):
            for p in range(1, 9):
                ratio = values[i, p - 1]
                text = f"{ratio:.2f}" if ratio < 10 else f"{ratio:.0f}"
                if f"{cell}_p{p}_{form}" in fallback:
                    text += "*"
                ax.text(
                    p - 1,
                    i,
                    text,
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white" if abs(math.log2(ratio)) > 2.3 else "black",
                )
        ax.set_xticks(range(8), range(1, 9))
        ax.set_yticks(
            range(12), [f"{cell}: {label}" for cell in ("Quad", "Hex") for label in LABELS]
        )
        ax.axhline(5.5, color="black", linewidth=1.5)
        ax.set_xlabel("Polynomial degree p")
        ax.set_title("-O2 -march=native" if flags == "O2" else "-O3 -march=native -ffast-math")
    bar = fig.colorbar(image, ax=axes, shrink=0.7, ticks=range(-3, 4))
    bar.ax.set_yticklabels(["≤0.125", "0.25", "0.5", "1", "2", "4", "≥8"])
    bar.set_label("FFCx time / TSFC time (lower is better for FFCx)")
    fig.suptitle(
        "FFCx sum-factorisation-o2 vs development TSFC spectral\n"
        "* FFCx stack-budget fallback; median of 3 rounds, each best of 7 samples",
        fontsize=13,
    )
    fig.savefig(HERE / "comparison.png", dpi=180)
    fig.savefig(HERE / "comparison.svg")
    plt.close(fig)

    report = [
        "# FFCx vs development TSFC: quad and hex cell kernels",
        "",
        "Date: 25 October 2026. These are standalone, hot element-kernel timings; "
        "they do not include Python, mesh traversal, assembly, or solver costs.",
        "",
        f"FFCx: `dokken/sum-factorisation-o2`, commit `{manifest['ffcx_commit']}`. "
        f"TSFC: Firedrake development HEAD `{manifest['tsfc_commit']}`. "
        f"FIAT/FInAT/GEM development HEAD `{manifest['fiat_finat_gem_commit']}`. "
        "Both development revisions were checked against upstream HEAD before generation.",
        "",
        "## Method",
        "",
        "- Intel Xeon w7-3445; GCC 14.3.0; one process pinned to CPU 2.",
        "- Scalar float64 GLL Lagrange Q_p bases; p=1-8; bilinear Q1 quads and trilinear Q1 hexes with non-affine geometry.",
        "- Degree 2p Gauss-Legendre quadrature, giving p+1 points in every direction. This is identical for both compilers; no GLL collocation shortcut.",
        "- FFCx uses `sum_factorization=True`. TSFC uses development `spectral` mode on native quad/hex cells. FInAT's converted hex basis has three independent index extents `(p+1, p+1, p+1)`, asserted for every degree. No extruded-cell replacement or release TSFC is used.",
        "- Mass and Laplace bilinear forms, vector Laplace, elasticity `inner(sym(grad(u)), sym(grad(v))) + div(u)*div(v)`, and mass/Laplace coefficient actions.",
        "- Actions use the same nonconstant scalar field, sampled at matching reference nodes, rather than a constant input.",
        "- Same C driver, matching geometry and coefficient values, separate kernel/driver translation units, no LTO. Kernel calls accumulate into an already allocated tensor; allocation and zeroing are excluded from timing.",
        "- Each executable calibrates to at least 10 ms per sample and takes 7 samples. Three complete sweeps shuffle compiler/case/flag ordering. Reported time is the median of the three per-sweep minima. Compilation and correctness validation finish before timing starts.",
        "- CPU frequency is not fixed. Across the three sweeps, the median per-executable range "
        f"(maximum round minimum / minimum round minimum - 1) is {statistics.median(spreads):.1%}; "
        f"the 95th percentile is {np.percentile(spreads, 95):.1%}. Near-parity differences should be treated accordingly.",
        "- Standalone TSFC imports use no-op citation registration and an in-process memoisation shim for `pyop3.cache.serial_cache`. TSFC, FInAT and GEM code generation and optimization passes are unmodified; no Firedrake/PETSc runtime is timed.",
        "- Process stack limit is 64 MiB for both compilers. Development TSFC's p8 hex elasticity kernel segfaults with the default 8 MiB stack; its O2 frame reserves about 10.15 MiB. FFCx retains its internal 256 KiB buffer guard and falls back to fused quadrature for hex elasticity at p7 and p8.",
        "",
        "## Correctness",
        "",
        f"All {len(validation)} complete element tensor comparisons passed after mapping scalar basis and interleaved vector component orderings. "
        f"Maximum absolute difference divided by the largest reference tensor entry: {max(v['relative_max_error'] for v in validation):.3g}. "
        f"Maximum relative L2 difference: {max(v['relative_l2_error'] for v in validation):.3g}.",
        "",
        "## Aggregate comparison",
        "",
        "Ratios are FFCx time / TSFC time: below 1 favors FFCx. Geometric means give equal weight to each of the 48 sampled form/degree cases per cell and flag set; they are not application speedups.",
        "",
        "| Cell | Flags | FFCx faster | Geometric mean ratio | Directional cases only | Ratio range |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for a in aggregates:
        report.append(
            f"| {a['cell']} | {a['flags']} | {a['ffcx_faster_cases']}/{a['cases']} | {a['geometric_mean_ffcx_over_tsfc']:.3f} | {a['geometric_mean_directional_only']:.3f} | {a['min_ffcx_over_tsfc']:.3f}-{a['max_ffcx_over_tsfc']:.3f} |"
        )
    report += [
        "",
        "![Timing ratio heatmap](comparison.png)",
        "",
        "## Per-case times",
        "",
        "Each entry is **FFCx / TSFC** in microseconds per kernel call. `*` marks FFCx's fused-quadrature fallback.",
    ]
    for cell in ("quadrilateral", "hexahedron"):
        for flags in ("O2", "O3fast"):
            report += [
                "",
                f"### {cell}, {flags}",
                "",
                "| p | " + " | ".join(LABELS) + " |",
                "|---:" + "|---:" * 6 + "|",
            ]
            for p in range(1, 9):
                entries = []
                for form in FORMS:
                    r = lookup[flags, cell, form, p]
                    marker = "*" if f"{cell}_p{p}_{form}" in fallback else ""
                    entries.append(f"{r['ffcx_us']:.3g} / {r['tsfc_us']:.3g}{marker}")
                report.append(f"| {p} | " + " | ".join(entries) + " |")
    report += [
        "",
        "## Reproduction and artifacts",
        "",
        "`environment.json` records revisions, package versions, hardware, compiler flags, interpreter paths, cache shims and work directory. "
        "`generate.py` generates each compiler's kernels; `run.py` compiles, checks and benchmarks them; `driver.c` is the shared timing driver; "
        "`report.py` regenerates this report and figures. "
        "`generated-kernels.tar.gz` contains the exact generated sources and nodal metadata. "
        "`samples.json` preserves every timing sample and calibrated iteration count. "
        "`results.csv`, `validation.json`, and `factorisation-audit.json` preserve the summaries, numerical errors and source hashes.",
        "",
        "```bash",
        "# Re-time the exact captured kernels without regenerating them.",
        "WORK=/tmp/ffcx-tsfc-reproduce",
        'mkdir -p "$WORK"',
        'tar -xzf generated-kernels.tar.gz -C "$WORK"',
        'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$FFCX_PYTHON" run.py "$WORK" --cpu 2 --rounds 3',
        "",
        "# Use the interpreter paths and TSFC_PYTHONPATH recorded in environment.json.",
        "WORK=/tmp/ffcx-tsfc-dev-benchmark-20261025/full",
        '"$FFCX_PYTHON" generate.py ffcx "$WORK/ffcx"',
        'XDG_CACHE_HOME=/tmp/tsfc-benchmark-cache PYTHONPATH="$TSFC_PYTHONPATH" "$TSFC_PYTHON" generate.py tsfc "$WORK/tsfc"',
        'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$FFCX_PYTHON" run.py "$WORK" --cpu 2 --rounds 3',
        'MPLCONFIGDIR=/tmp/tsfc-benchmark-matplotlib "$PLOT_PYTHON" report.py',
        "```",
        "",
    ]
    (HERE / "report.md").write_text("\n".join(report))
    print(json.dumps(aggregates, indent=2))


if __name__ == "__main__":
    main()
