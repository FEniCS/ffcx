// Minimal example of calling an FFCx-generated kernel from C.
//
// Assembles the P1 Poisson stiffness matrix on a single triangle and
// compares it with the exact result.

#include "poisson.h"

#include <math.h>
#include <stdio.h>
#include <ufcx.h>

int main(void)
{
  // Pick the form's cell integral. A form holds its integrals grouped by
  // type, with form_integral_offsets giving the start of each group.
  const ufcx_form* a = form_poisson_a;
  const int offset = a->form_integral_offsets[cell];
  if (a->form_integral_offsets[cell + 1] - offset != 1)
  {
    fprintf(stderr, "expected exactly one cell integral\n");
    return 1;
  }
  ufcx_tabulate_tensor_float64* tabulate_tensor
      = a->form_integrals[offset]->tabulate_tensor_float64;

  // The reference triangle. Coordinates are always three components per
  // vertex, whatever the geometric dimension of the mesh.
  const double coords[9] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0};

  // A is the 3x3 element matrix, stored row-major. This form has no
  // coefficients or constants, so those arguments are unused.
  double A[9] = {0.0};
  tabulate_tensor(A, NULL, NULL, coords, NULL, NULL, NULL);

  const double expected[9] = {1.0, -0.5, -0.5, -0.5, 0.5, 0.0, -0.5, 0.0, 0.5};

  int status = 0;
  for (int i = 0; i < 9; ++i)
  {
    printf("%8.4f%s", A[i], (i % 3 == 2) ? "\n" : " ");
    if (fabs(A[i] - expected[i]) > 1e-12)
      status = 1;
  }

  if (status)
    fprintf(stderr, "element matrix does not match the expected values\n");

  return status;
}
