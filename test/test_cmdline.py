# Copyright (C) 2018-2025 Chris N. Richardson and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import importlib
import shutil
import subprocess
from pathlib import Path

import pytest

import ffcx.codegeneration


def test_cmdline_simple():
    dir = Path(__file__).parent
    subprocess.run(["ffcx", "-o", dir, dir / "poisson.py"], check=True)


@pytest.mark.skipif(
    importlib.util.find_spec("pygraphviz") is None, reason="pygraphviz not installed."
)
@pytest.mark.parametrize("args,suffix", [([], "pdf"), (["--visualise-format", "png"], "png")])
def test_visualise(tmp_path, args, suffix):
    """Each integrand's graphs are named after the integral they belong to.

    The graphs are written relative to the working directory, so the
    command is run inside `tmp_path`.
    """
    dir = Path(__file__).parent
    subprocess.run(
        ["ffcx", "-d", str(tmp_path), "--visualise", *args, str(dir / "poisson.py")],
        check=True,
        cwd=tmp_path,
    )
    s_graphs = sorted(tmp_path.glob(f"*_S.{suffix}"))
    f_graphs = sorted(tmp_path.glob(f"*_F.{suffix}"))
    assert len(s_graphs) == len(f_graphs) > 0

    # poisson.py holds a bilinear and a linear form, which must not
    # overwrite each other's graphs
    assert len(s_graphs) > 1
    assert all(p.stat().st_size > 0 for p in s_graphs + f_graphs)

    # The old, fixed names are no longer used
    assert not (tmp_path / f"S.{suffix}").exists()


def test_visualise_format_normalisation():
    """`True` is still accepted, and means the default format."""
    from ffcx.compiler import _visualise_format

    assert _visualise_format(True) == "pdf"
    assert _visualise_format(False) is None
    assert _visualise_format(None) is None
    assert _visualise_format("png") == "png"


def test_generated_code_without_complex(tmp_path):
    """Generated code must compile where C99 complex is unavailable.

    `ufcx.h` declares the complex kernel members only when the compiler
    supports them, so the generated initialisers have to be guarded the
    same way. MSVC is the case that matters.
    """
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if cc is None:
        pytest.skip("no C compiler found.")

    dir = Path(__file__).parent
    subprocess.run(["ffcx", "-d", str(tmp_path), str(dir / "poisson.py")], check=True)

    include = Path(ffcx.codegeneration.__file__).parent
    for extra in ([], ["-D__STDC_NO_COMPLEX__"]):
        subprocess.run(
            [
                cc,
                "-c",
                f"-I{include}",
                f"-I{tmp_path}",
                *extra,
                str(tmp_path / "poisson.c"),
                "-o",
                str(tmp_path / "poisson.o"),
            ],
            check=True,
        )
