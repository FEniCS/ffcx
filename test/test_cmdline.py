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


def test_generated_code_included_in_cpp(tmp_path):
    """Generated code must compile when included in C++ without C99 complex.

    DOLFINx's custom kernel demo includes generated C in a C++ file, which
    MSVC compiles without `_Complex`. `-D__STDC_NO_COMPLEX__` reproduces
    that on other compilers.
    """
    cxx = shutil.which("c++") or shutil.which("g++") or shutil.which("clang++")
    if cxx is None:
        pytest.skip("no C++ compiler found.")

    dir = Path(__file__).parent
    subprocess.run(["ffcx", "-d", str(tmp_path), str(dir / "poisson.py")], check=True)

    # ufcx.h undefines its restrict workaround, so define it after the header
    wrapper = tmp_path / "wrapper.cpp"
    wrapper.write_text(
        '#include "poisson.h"\n#define restrict __restrict__\n#include "poisson.c"\n'
    )

    include = Path(ffcx.codegeneration.__file__).parent
    for extra in ([], ["-D__STDC_NO_COMPLEX__"]):
        subprocess.run(
            [
                cxx,
                "-std=c++20",
                "-c",
                f"-I{include}",
                f"-I{tmp_path}",
                *extra,
                str(wrapper),
                "-o",
                str(tmp_path / "wrapper.o"),
            ],
            check=True,
        )
