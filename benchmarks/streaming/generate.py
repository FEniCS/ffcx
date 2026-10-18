"""Generate benchmark forms for an experimental compiler variant."""

# Generated C fragments retain their original source layout.

import argparse
import os
import sys
from pathlib import Path

from benchmarks.elasticity import generate

root = Path(os.environ.get("FFCX_STREAM_WORKDIR", "/tmp/ffcx-streaming-investigation"))
mode = sys.argv[1]
out = root / mode
out.mkdir(exist_ok=True)
if mode == "aliases":
    from aliases import install

    install()
if mode == "fullscaling":
    from scaling import install_full

    install_full()
if mode == "scaling":
    from scaling import install

    install()
# Reuse the exact original forms, quadrature and input construction.
for suite, source in [
    ("poisson", str(Path(__file__).with_name("generate_poisson.py"))),
    ("elasticity", str(Path(__file__).with_name("generate_simplex_elasticity.py"))),
]:
    path = out / suite
    path.mkdir(exist_ok=True)
    text = Path(source).read_text()
    old = (
        "/tmp/poisson-vs-reference-current"
        if suite == "poisson"
        else "/tmp/simplex-elasticity-reference"
    )
    exec(compile(text.replace(old, str(path)), source, "exec"), {"__name__": "__main__"})


generate(
    argparse.Namespace(
        directory=out / "hyper",
        revision=mode,
        degrees=[1, 2, 3],
        quadrature_from=None if mode == "baseline" else root / "baseline" / "hyper",
    )
)
