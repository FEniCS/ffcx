"""Test demos."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

import ffcx.main

demo_dir = Path(__file__).parent

ufl_files = [
    f
    for f in demo_dir.iterdir()
    if f.suffix == ".py" and not f.stem.endswith("_numba") and f != Path(__file__)
]

skip_complex = ["BiharmonicHHJ", "BiharmonicRegge", "StabilisedStokes"]


def skip_unsupported(test):
    """Dcecorate test case to skip unsupported cases."""

    def check_skip(file, scalar_type):
        """Skip scalar_type file combinations not supported."""
        if "complex" in scalar_type and file.stem in skip_complex:
            pytest.skip(reason="Not implemented for complex types")
        elif "Complex" in file.stem and scalar_type in ["float64", "float32"]:
            pytest.skip(reason="Not implemented for real types")

        return test(file, scalar_type)

    return check_skip


@pytest.mark.parametrize("file", ufl_files)
@pytest.mark.parametrize("scalar_type", ["float64", "float32", "complex128", "complex64"])
@skip_unsupported
def test_C(file, scalar_type):
    """Test a demo."""
    if sys.platform.startswith("win32") and "complex" in scalar_type:
        # Skip complex demos on win32
        pytest.skip(reason="_Complex not supported on Windows")

    assert ffcx.main.main(["--scalar_type", scalar_type, str(file)]) == 0

    if sys.platform.startswith("win32"):
        extra_flags = "/std:c17"
        for compiler in ["cl.exe", "clang-cl.exe"]:
            subprocess.run(
                [
                    compiler,
                    "/I",
                    f"{demo_dir.parent / 'ffcx/codegeneration'}",
                    *extra_flags.split(" "),
                    "/c",
                    file.with_suffix(".c"),
                ],
                cwd=demo_dir,
                check=True,
            )
    else:
        cc = os.environ.get("CC", "cc")
        extra_flags = (
            "-std=c17 -Wunused-variable -Werror -fPIC -Wno-error=implicit-function-declaration"
        )
        subprocess.run(
            [
                cc,
                f"-I{demo_dir.parent / 'ffcx/codegeneration'}",
                *extra_flags.split(" "),
                "-c",
                file.with_suffix(".c"),
            ],
            cwd=demo_dir,
            check=True,
        )


@pytest.mark.parametrize("file", ufl_files)
@pytest.mark.parametrize("scalar_type", ["float64", "float32", "complex128", "complex64"])
@skip_unsupported
def test_numba(file, scalar_type):
    """Test numba generation."""
    opts = f"--language numba --scalar_type {scalar_type}"
    assert ffcx.main.main([*opts.split(" "), str(file)]) == 0
    subprocess.run(["python", file], cwd=demo_dir, check=True)
