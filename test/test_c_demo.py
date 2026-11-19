# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Build and run the minimal CMake + C example."""

import shutil
import subprocess
from pathlib import Path

import pytest

repo_dir = Path(__file__).parent.parent
demo_dir = repo_dir / "demo" / "c"


@pytest.mark.skipif(shutil.which("cmake") is None, reason="CMake not found.")
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
    subprocess.run(["cmake", "--build", str(build), "--config", "Release"], check=True)

    # Run it through CTest rather than by path: multi-config generators put the
    # binary in a per-config subdirectory, and Windows adds a suffix. The
    # example returns non-zero if the element matrix is wrong.
    subprocess.run(
        [
            "ctest",
            "--test-dir",
            str(build),
            "--build-config",
            "Release",
            "--output-on-failure",
        ],
        check=True,
    )
