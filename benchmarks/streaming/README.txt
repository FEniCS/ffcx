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
