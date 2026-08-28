# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Run the generated-code performance benchmark suite and print a results table.

Usage: python -m bench.run [name-substring ...]

With no arguments, runs the full suite (bench/kernels.py). Otherwise runs
only cases whose name contains one of the given substrings. The timing loop
currently requires a POSIX C toolchain.
"""

import sys
import tempfile
from pathlib import Path

from bench.harness import _build_timer_library, run_case
from bench.kernels import suite


def main(argv: list[str]) -> None:
    """Run the benchmark suite, filtered by optional name substrings."""
    cases = suite()
    if argv:
        cases = [c for c in cases if any(pat in c.name for pat in argv)]
        if not cases:
            print(f"No benchmark case name matches any of {argv!r}")
            return

    timer = _build_timer_library()

    # FFCx's JIT cache key describes the form and compilation options, not
    # the FFCx source revision. A fresh cache per benchmark invocation avoids
    # accidentally timing a compiled module from a branch checked out earlier
    # in the same worktree, while still sharing modules across this suite run.
    with tempfile.TemporaryDirectory(prefix="ffcx-bench-") as tmp:
        cache_dir = Path(tmp)
        name_width = max(len(c.name) for c in cases)
        print(f"{'case':<{name_width}}  {'ns/cell':>12}  {'stack (B)':>10}  {'tensor size':>11}")
        for case in cases:
            result = run_case(case, timer, cache_dir)
            stack = "n/a" if result.stack_bytes is None else str(result.stack_bytes)
            print(
                f"{result.name:<{name_width}}  {result.ns_per_call:>12.1f}  "
                f"{stack:>10}  {result.tensor_size:>11}"
            )


if __name__ == "__main__":
    main(sys.argv[1:])
