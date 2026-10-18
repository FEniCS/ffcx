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
def test_visualise():
    dir = Path(__file__).parent
    subprocess.run(["ffcx", "-o", dir, "--visualise", dir / "poisson.py"], check=True)
    assert (Path.cwd() / "S.pdf").is_file()
    assert (Path.cwd() / "F.pdf").is_file()


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
