from ufl import FiniteElement
from ufl import interval, triangle, tetrahedron
from ffc.compiler import compile_element

# Compile Lagrane elements of differing degress in 1D, 2D and 3D
f = open('lagrange_elements.h', 'w')
#for cell in (interval, triangle, tetrahedron):
for cell in (triangle,):
    for p in range(1, 2):
        element = FiniteElement("Lagrange", cell, p)
        header, implementation = compile_element(element,
                                                 prefix=str(cell) + "_" + str(p))

        f.write(header)

f.close()
