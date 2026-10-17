# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Build and run the minimal CMake + C example."""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

repo_dir = Path(__file__).parent.parent
demo_dir = repo_dir / "demo" / "c"


@pytest.mark.skipif(shutil.which("cmake") is None, reason="CMake not found.")
@pytest.mark.skipif(
    sys.platform.startswith("win32"),
    reason="Generated code initialises the complex kernel members "
    "unconditionally, which MSVC rejects. Remove this skip once #917 lands.",
)
def test_c_demo(tmp_path):
    """The example builds against an installed UFCx and gives the right matrix."""
    prefix = tmp_path / "ufcx"
    subprocess.run(
        [
            "cmake",
            "-S",
            str(repo_dir / "cmake"),
            "-B",
            str(tmp_path / "build-ufcx"),
            f"-DCMAKE_INSTALL_PREFIX={prefix}",
        ],
        check=True,
    )
    subprocess.run(["cmake", "--install", str(tmp_path / "build-ufcx")], check=True)

    build = tmp_path / "build"
    subprocess.run(
        ["cmake", "-S", str(demo_dir), "-B", str(build), f"-DCMAKE_PREFIX_PATH={prefix}"],
        check=True,
    )
    subprocess.run(["cmake", "--build", str(build)], check=True)

    # The example returns non-zero if the element matrix is wrong
    subprocess.run([str(build / "poisson")], check=True)
