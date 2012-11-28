#include <dolfin.h>
#include "ProjectionManifold.h"

using namespace dolfin;

int main()
{
  // Create mesh and function space
  UnitCircleMesh mesh(16);
  BoundaryMesh boundary(mesh);

  ProjectionManifold::FunctionSpace V(boundary);

  // Define variational forms
  ProjectionManifold::BilinearForm a(V, V);
  Matrix A;
  assemble(A, a);

  info(A, true);
  plot(mesh);

  interactive();

  return 0;
}
