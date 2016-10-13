# -*- coding: utf-8 -*-

import re
from ffc.backends.ufc import *

# TODO: Make cell_orientation a double +1.0|-1.0 instead of the int flag in ffc/ufc/dolfin
# TODO: Simplify ufc templates by introducing 'preamble' keyword in place of members, constructor, destructor

domain_background = """
/// This is just here to document the memory layout of the geometry data arrays
struct geometry_data
{
  // Example dimensions
  std::size_t gdim = 3;
  std::size_t tdim = 2;
  std::size_t num_points = 1;

  // Memory layout of geometry data arrays
  double x[num_points * gdim];         // x[i]   -> x[ip*gdim + i]
  double X[num_points * tdim];         // X[j]   -> X[ip*tdim + j]
  double J[num_points * gdim * tdim];  // J[i,j] -> J[ip*gdim*tdim + i*tdim + j]
  double detJ[num_points];             // detJ   -> detJ[ip]
  double K[num_points * tdim * gdim];  // K[j,i] -> K[ip*tdim*gdim + j*gdim + i]
  double n[num_points * gdim];         // n[i]   -> n[ip*gdim + i]

  // In the affine case we have the relation:
  // x[i] = x0[i] + sum_j J[i,j] X[j]
  // X[j] = sum_i K[j,i] (x[i] - x0[i])

};
"""

def extract_keywords(template):
    r = re.compile(r"%\(([a-zA-Z0-9_]*)\)")
    return set(r.findall(template))
