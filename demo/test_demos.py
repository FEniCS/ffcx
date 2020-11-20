import os
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
                     "CPATH=/usr/local/lib/python3.8/dist-packages/ffcx/codegeneration/ "
                     "gcc -I/usr/include/python3.8 -fPIC -shared {file}.c "
                     "-o {file}.so")
