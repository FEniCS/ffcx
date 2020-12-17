#include "Grad.h"
#include <cmath>
#include <dolfinx.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/fem/Constant.h>
#include <stdio.h>

using namespace dolfinx;

int main(int argc, char* argv[])
{
  common::SubSystemsManager::init_logging(argc, argv);
  common::SubSystemsManager::init_petsc(argc, argv);

  {
#define INPUT_DOFS 3
#define OUTPUT_DOFS 12

      double A[INPUT_DOFS*OUTPUT_DOFS];
      const double cdofs[6] = {0.0, 0.0, 1.0, 0.0, 0.0, 1.0};
      ufc_expression* expr = create_expression();
      expr->tabulate_expression(A, nullptr, nullptr, cdofs);
      free(expr);

      for (int i = 0; i < OUTPUT_DOFS; ++i)
      {
        for (int j = 0; j < INPUT_DOFS; ++j)
        {
            printf("A[%d][%d]=%f \n", i, j, A[INPUT_DOFS*i + j]);
        }
      }

  }

  common::SubSystemsManager::finalize_petsc();
  return 0;
}