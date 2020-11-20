import os
import sys
import pytest

demo_dir = os.path.dirname(os.path.realpath(__file__))

ufl_files = []
for file in os.listdir(demo_dir):
    if file.endswith(".ufl") and not file.startswith("."):
        ufl_files.append(file[:-4])


@pytest.mark.parametrize("file", ufl_files)
def test_demo(file):
    assert os.system(f"cd {demo_dir} && ffcx {file}.ufl") == 0
    assert os.system(f"cd {demo_dir} && "
                     "CPATH=../ffcx/codegeneration/ "
                     f"gcc -I/usr/include/python{sys.version_info.major}.{sys.version_info.minor} -fPIC "
                     f"-shared {file}.c -o {file}.so") == 0
