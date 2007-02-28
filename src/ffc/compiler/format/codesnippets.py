"Code snippets for code generation"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-28 -- 2007-02-28"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Code snippet for computing the Jacobian, its inverse and determinant in 2D
jacobian_2D = """\
// Compute Jacobian of affine map from reference cell
const double J%(restriction)s00 = c.coordinates[1][0] - c.coordinates[0][0];
const double J%(restriction)s01 = c.coordinates[2][0] - c.coordinates[0][0];
const double J%(restriction)s10 = c.coordinates[1][1] - c.coordinates[0][1];
const double J%(restriction)s11 = c.coordinates[2][1] - c.coordinates[0][1];
  
// Compute determinant
double det = J%(restriction)s00*J%(restriction)s11 - J%(restriction)s01*J%(restriction)s10;
  
// Compute inverse of Jacobian
const double J%(restriction)sinv00 =  J%(restriction)s11 / det;
const double J%(restriction)sinv01 = -J%(restriction)s01 / det;
const double J%(restriction)sinv10 = -J%(restriction)s10 / det;
const double J%(restriction)sinv11 =  J%(restriction)s00 / det;

// Take absolute value of determinant
det = std::abs(det);
"""
# Code snippet for computing the Jacobian, its inverse and determinant in 3D
jacobian_3D = """\
// Compute Jacobian of affine map from reference cell
const double J%(restriction)s00 = c.coordinates[1][0] - c.coordinates[0][0];
const double J%(restriction)s01 = c.coordinates[2][0] - c.coordinates[0][0];
const double J%(restriction)s02 = c.coordinates[3][0] - c.coordinates[0][0];
const double J%(restriction)s10 = c.coordinates[1][1] - c.coordinates[0][1];
const double J%(restriction)s11 = c.coordinates[2][1] - c.coordinates[0][1];
const double J%(restriction)s12 = c.coordinates[3][1] - c.coordinates[0][1];
const double J%(restriction)s20 = c.coordinates[1][2] - c.coordinates[0][2];
const double J%(restriction)s21 = c.coordinates[2][2] - c.coordinates[0][2];
const double J%(restriction)s22 = c.coordinates[3][2] - c.coordinates[0][2];
  
// Compute sub determinants
const double d00 = J%(restriction)s11*J%(restriction)s22 - J%(restriction)s12*J%(restriction)s21;
const double d01 = J%(restriction)s12*J%(restriction)s20 - J%(restriction)s10*J%(restriction)s22;
const double d02 = J%(restriction)s10*J%(restriction)s21 - J%(restriction)s11*J%(restriction)s20;

const double d10 = J%(restriction)s02*J%(restriction)s21 - J%(restriction)s01*J%(restriction)s22;
const double d11 = J%(restriction)s00*J%(restriction)s22 - J%(restriction)s02*J%(restriction)s20;
const double d12 = J%(restriction)s01*J%(restriction)s20 - J%(restriction)s00*J%(restriction)s21;

const double d20 = J%(restriction)s01*J%(restriction)s12 - J%(restriction)s02*J%(restriction)s11;
const double d21 = J%(restriction)s02*J%(restriction)s10 - J%(restriction)s00*J%(restriction)s12;
const double d22 = J%(restriction)s00*J%(restriction)s11 - J%(restriction)s01*J%(restriction)s10;
  
// Compute determinant
double det = J%(restriction)s00*d00 + J%(restriction)s10*d10 + J%(restriction)s20*d20;
  
// Compute inverse of Jacobian
const double J%(restriction)s00 = d00 / det;
const double J%(restriction)s01 = d10 / det;
const double J%(restriction)s02 = d20 / det;
const double J%(restriction)s10 = d01 / det;
const double J%(restriction)s11 = d11 / det;
const double J%(restriction)s12 = d21 / det;
const double J%(restriction)s20 = d02 / det;
const double J%(restriction)s21 = d12 / det;
const double J%(restriction)s22 = d22 / det;

// Take absolute value of determinant
det = std::abs(det);
"""
