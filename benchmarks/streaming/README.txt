Experimental streaming investigation checkpoint

These scripts do not change default code generation. Symmetry prototypes assume
that the supplied benchmark form is symmetric; they are not general optimisers.
The scalar-scaling and alias passes are experimental monkeypatches. None has
passed a broad no-regression acceptance check yet.

Run from the repository root using a Python environment with FFCx dependencies:

  export PYTHONPATH="$PWD"
  export FFCX_STREAM_WORKDIR=/tmp/ffcx-streaming-investigation
  python benchmarks/streaming/generate.py baseline
  python benchmarks/streaming/generate.py scaling
  python benchmarks/streaming/generate.py fullscaling
  python benchmarks/streaming/generate.py aliases
  python benchmarks/streaming/build.py
  python benchmarks/streaming/add_symmetry.py
  python benchmarks/streaming/add_blocksym.py

Screen one family of variants, pinned to CPU 2 by bench.py:

  STREAM_MODES=baseline,aligned,balanced STREAM_OUTPUT=screen.json \
    python benchmarks/streaming/bench.py

Use STREAM_MIB=1024 and STREAM_REPS=9 for larger-array repeats. STREAM_FILTER is
a regular expression selecting case names. Run only one timing process at once.
Checks include all stored output entries, arbitrary nonzero output accumulation,
and reflected geometry. The earlier Poisson/elasticity results used 1 GiB per
variant and nine samples; the new exploratory screens used 256 MiB and three.

Findings so far:
- Width-2 tiling and staging have substantial regressions on some forms.
- Factoring Jacobian scaling, aligning tables and balancing sums have mixed,
  mostly small gains. Several possible gains still need independent repeats.
- Resolving scalar aliases reduces some contraction terms, but has not shown
  a consistent streaming gain.
- Direct triangular symmetry reuse often hurts vectorisation; it helps some
  cases substantially but regresses others. Blocked variants are the next test.

The results file records completed experiments and pending validation. This is
an investigation checkpoint, not evidence of a regression-free optimisation.

Memory-traffic bound (roofline.py)

roofline.py times generated simplex kernels for Poisson, linear elasticity and
hyperelasticity over DRAM-resident cell records, together with a null kernel
that reads every input and updates every output entry without arithmetic. The
null time is the per-core DRAM bound for the record layout. Each kernel is also
timed on one record held in L1 ("hot"), which isolates its arithmetic cost.
All output entries are checked, with zero and nonzero initial output, against
independent quadrature (Poisson, elasticity) or the first variant (hyper):

  python benchmarks/streaming/generate.py head        # in each revision
  python benchmarks/streaming/roofline.py --variant head=$FFCX_STREAM_WORKDIR/head \
      --variant new=$FFCX_STREAM_WORKDIR/new --mib 1024 --reps 7 --output roofline.json

The driver zeroes each output record before calling the kernel. Arithmetic
performed before a kernel first touches A largely overlaps with that traffic,
while arithmetic interleaved with or after updates of A adds to it.

Refreshed performance tables

The 20 October 2026 refresh is pinned to kernel source 9e736b9f. Its streaming
baseline is 30757e95, which already includes local block accumulation and
point-major evaluation. The comparison therefore measures the later row-blocked
updates and symmetric-diagonal contraction, not the complete PR versus main.
Both streaming variants use matching dependencies, quadrature and input arrays.
The results retain individual streaming/hot samples and numerical errors.
benchmarks/refresh-current-results.json maps all four published tables to their
current measurements; refresh-results.json retains the previous refresh.
