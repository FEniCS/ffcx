import os
import sys
import pytest
from ffcx.main import main

demo_dir = os.path.dirname(os.path.realpath(__file__))

ufl_files = []
for file in os.listdir(demo_dir):
    if file.endswith(".ufl") and not file.startswith("."):
        ufl_files.append(file[:-4])


@pytest.mark.parametrize("file", ufl_files)
def test_demo(file):
    if file in [
        "MixedPoissonDual",  # Discontinuous Raviart-Thomas
        "BiharmonicHHJ",  # Hellan-Herrmann-Johnson
        "NodalMini",  # NodalEnrichedElement
        "Mini",  # EnrichedElement
        "MixedGradient", "TraceElement",  # HDiv Trace
        "MassHdiv_2D_1", "MassHdiv_2D_3", "MixedPoisson", "MassHdiv_2D_2",  # Brezzi-Douglas-Marini
    ]:
        # Skip demos that use elements not yet implemented in basix
        pytest.skip()

    os.chdir(f"{demo_dir}")
    main([f"{file}.ufl"])
#    assert os.system(f"cd {demo_dir} && ffcx {file}.ufl") == 0

    assert os.system(f"cd {demo_dir} && "
                     "CPATH=../ffcx/codegeneration/ "
                     f"gcc -I/usr/include/python{sys.version_info.major}.{sys.version_info.minor} -fPIC "
                     f"-shared {file}.c -o {file}.so") == 0
