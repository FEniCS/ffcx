#include <dolfin.h>
#include "ProjectionManifold.h"

using namespace dolfin;

int main()
{
  // Create mesh and function space
  UnitCircleMesh mesh(2);
  BoundaryMesh boundary(mesh);

  ProjectionManifold::FunctionSpace V(boundary);

  // Define variational forms
  ProjectionManifold::BilinearForm a(V, V);
  ProjectionManifold::LinearForm L(V);

  Function f(V);
  (*f.vector()) = 1.0;
  L.f = f;

  // Assemble
  Matrix A;
  Vector b;
  assemble(A, a);
  assemble(b, L);

  // Solve
  Vector x;
  solve(A, x, b);
  info(x, true);

  return 0;
}
