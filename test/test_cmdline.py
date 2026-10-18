# Copyright (C) 2018-2025 Chris N. Richardson and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import importlib
import shutil
import subprocess
import sys
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


def _find_compiler(cpp: bool) -> str | None:
    """Find a C or C++ compiler, preferring MSVC on Windows."""
    names = ["cl"] if sys.platform == "win32" else []
    names += ["c++", "g++", "clang++"] if cpp else ["cc", "gcc", "clang"]
    return next((n for n in names if shutil.which(n)), None)


def _compile(compiler: str, source: Path, include_dirs: list[Path], defines: list[str], std: str):
    if Path(compiler).stem == "cl":
        cmd = [compiler, "/nologo", "/c", f"/std:{std}", f"/Fo{source.with_suffix('.obj')}"]
        cmd += [f"/I{d}" for d in include_dirs] + [f"/D{d}" for d in defines]
    else:
        cmd = [compiler, "-c", f"-std={std}", "-o", str(source.with_suffix(".o"))]
        cmd += [f"-I{d}" for d in include_dirs] + [f"-D{d}" for d in defines]
    result = subprocess.run([*cmd, str(source)], capture_output=True, text=True)
    if result.returncode != 0:
        pytest.fail(f"{' '.join(cmd)} failed:\n{result.stdout}{result.stderr}")


def test_generated_code_without_complex(tmp_path):
    """Generated code must compile where C99 complex is unavailable.

    `ufcx.h` declares the complex kernel members only when the compiler
    supports them, so the generated initialisers have to be guarded the
    same way. MSVC is the case that matters.
    """
    cc = _find_compiler(cpp=False)
    if cc is None:
        pytest.skip("no C compiler found.")

    dir = Path(__file__).parent
    subprocess.run(["ffcx", "-d", str(tmp_path), str(dir / "poisson.py")], check=True)

    include = Path(ffcx.codegeneration.__file__).parent
    msvc = Path(cc).stem == "cl"
    for defines in [[]] if msvc else [[], ["__STDC_NO_COMPLEX__"]]:
        _compile(cc, tmp_path / "poisson.c", [include, tmp_path], defines, "c17")


def test_generated_code_included_in_cpp(tmp_path):
    """Generated code must compile when included in C++.

    DOLFINx's custom kernel demo includes generated C in a C++ file. MSVC
    has no `_Complex`, and `__STDC_NO_COMPLEX__` takes the same path on
    other compilers.
    """
    cxx = _find_compiler(cpp=True)
    if cxx is None:
        pytest.skip("no C++ compiler found.")

    dir = Path(__file__).parent
    subprocess.run(["ffcx", "-d", str(tmp_path), str(dir / "poisson.py")], check=True)

    # ufcx.h undefines its restrict workaround, so define it after the header
    wrapper = tmp_path / "wrapper.cpp"
    wrapper.write_text(
        '#include "poisson.h"\n'
        "#ifdef _MSC_VER\n#define restrict __restrict\n"
        "#else\n#define restrict __restrict__\n#endif\n"
        '#include "poisson.c"\n'
    )

    include = Path(ffcx.codegeneration.__file__).parent
    msvc = Path(cxx).stem == "cl"
    for defines in [[]] if msvc else [[], ["__STDC_NO_COMPLEX__"]]:
        _compile(cxx, wrapper, [include, tmp_path], defines, "c++20")
