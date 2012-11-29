#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

// Source term (right-hand side)
class Source : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double dx = x[0] - 0.5;
    double dy = x[1] - 0.5;
    values[0] = 10*exp(-(dx*dx + dy*dy) / 0.02);
  }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS;
  }
};

int main()
{
  // Create mesh and function space
  UnitCircleMesh mesh(16);
  BoundaryMesh boundary(mesh);

  Poisson::FunctionSpace V(boundary);

  // Define variational forms
  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  Source f;
  L.f = f;

  // Define some "boundary" (not) condition to non-singularize.
  Constant zero(0.0);
  DirichletBoundary hack;
  DirichletBC bc(V, zero, hack);

  // Assemble
  Matrix A;
  Vector b;
  assemble(A, a);
  assemble(b, L);
  bc.apply(A, b);

  // Solve
  Vector x;
  solve(A, x, b);
  info(x, true);

  return 0;
}
