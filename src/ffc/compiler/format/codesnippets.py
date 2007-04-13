"Code snippets for code generation"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-28 -- 2007-04-04"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Code snippet for computing the Jacobian, its inverse and determinant in 2D
jacobian_2D = """\
// Extract vertex coordinates
const double * const * x%(restriction)s = c%(restriction)s.coordinates;

// Compute Jacobian of affine map from reference cell
const double J%(restriction)s_00 = x%(restriction)s[1][0] - x%(restriction)s[0][0];
const double J%(restriction)s_01 = x%(restriction)s[2][0] - x%(restriction)s[0][0];
const double J%(restriction)s_10 = x%(restriction)s[1][1] - x%(restriction)s[0][1];
const double J%(restriction)s_11 = x%(restriction)s[2][1] - x%(restriction)s[0][1];
  
// Compute determinant of Jacobian
double detJ%(restriction)s = J%(restriction)s_00*J%(restriction)s_11 - J%(restriction)s_01*J%(restriction)s_10;
  
// Compute inverse of Jacobian
const double Jinv%(restriction)s_00 =  J%(restriction)s_11 / detJ%(restriction)s;
const double Jinv%(restriction)s_01 = -J%(restriction)s_01 / detJ%(restriction)s;
const double Jinv%(restriction)s_10 = -J%(restriction)s_10 / detJ%(restriction)s;
const double Jinv%(restriction)s_11 =  J%(restriction)s_00 / detJ%(restriction)s;

// Take absolute value of determinant
detJ%(restriction)s = std::abs(detJ%(restriction)s);
"""

# Code snippet for computing the Jacobian, its inverse and determinant in 3D
jacobian_3D = """\
// Extract vertex coordinates
const double * const * x%(restriction)s = c%(restriction)s.coordinates;

// Compute Jacobian of affine map from reference cell
const double J%(restriction)s_00 = x%(restriction)s[1][0] - x%(restriction)s[0][0];
const double J%(restriction)s_01 = x%(restriction)s[2][0] - x%(restriction)s[0][0];
const double J%(restriction)s_02 = x%(restriction)s[3][0] - x%(restriction)s[0][0];
const double J%(restriction)s_10 = x%(restriction)s[1][1] - x%(restriction)s[0][1];
const double J%(restriction)s_11 = x%(restriction)s[2][1] - x%(restriction)s[0][1];
const double J%(restriction)s_12 = x%(restriction)s[3][1] - x%(restriction)s[0][1];
const double J%(restriction)s_20 = x%(restriction)s[1][2] - x%(restriction)s[0][2];
const double J%(restriction)s_21 = x%(restriction)s[2][2] - x%(restriction)s[0][2];
const double J%(restriction)s_22 = x%(restriction)s[3][2] - x%(restriction)s[0][2];
  
// Compute sub determinants
const double d00 = J%(restriction)s_11*J%(restriction)s_22 - J%(restriction)s_12*J%(restriction)s_21;
const double d01 = J%(restriction)s_12*J%(restriction)s_20 - J%(restriction)s_10*J%(restriction)s_22;
const double d02 = J%(restriction)s_10*J%(restriction)s_21 - J%(restriction)s_11*J%(restriction)s_20;

const double d10 = J%(restriction)s_02*J%(restriction)s_21 - J%(restriction)s_01*J%(restriction)s_22;
const double d11 = J%(restriction)s_00*J%(restriction)s_22 - J%(restriction)s_02*J%(restriction)s_20;
const double d12 = J%(restriction)s_01*J%(restriction)s_20 - J%(restriction)s_00*J%(restriction)s_21;

const double d20 = J%(restriction)s_01*J%(restriction)s_12 - J%(restriction)s_02*J%(restriction)s_11;
const double d21 = J%(restriction)s_02*J%(restriction)s_10 - J%(restriction)s_00*J%(restriction)s_12;
const double d22 = J%(restriction)s_00*J%(restriction)s_11 - J%(restriction)s_01*J%(restriction)s_10;
  
// Compute determinant of Jacobian
double detJ%(restriction)s = J%(restriction)s_00*d00 + J%(restriction)s_10*d10 + J%(restriction)s_20*d20;
  
// Compute inverse of Jacobian
const double Jinv%(restriction)s_00 = d00 / detJ%(restriction)s;
const double Jinv%(restriction)s_01 = d10 / detJ%(restriction)s;
const double Jinv%(restriction)s_02 = d20 / detJ%(restriction)s;
const double Jinv%(restriction)s_10 = d01 / detJ%(restriction)s;
const double Jinv%(restriction)s_11 = d11 / detJ%(restriction)s;
const double Jinv%(restriction)s_12 = d21 / detJ%(restriction)s;
const double Jinv%(restriction)s_20 = d02 / detJ%(restriction)s;
const double Jinv%(restriction)s_21 = d12 / detJ%(restriction)s;
const double Jinv%(restriction)s_22 = d22 / detJ%(restriction)s;

// Take absolute value of determinant
detJ%(restriction)s = std::abs(detJ%(restriction)s);
"""

# Code snippet for computing the scale factor (determinant)
scale_factor = """\
// Set scale factor
const double det = detJ;
"""

# Code snippet for computing the determinant of the facet mapping in 2D
facet_determinant_2D = """\
// Vertices on edges
static unsigned int edge_vertices[3][2] = {{1, 2}, {0, 2}, {0, 1}};

// Get vertices
const unsigned int v0 = edge_vertices[%(facet)s][0];
const unsigned int v1 = edge_vertices[%(facet)s][1];

// Compute scale factor (length of edge scaled by length of reference interval)
const double dx0 = x%(restriction)s[v1][0] - x%(restriction)s[v0][0];
const double dx1 = x%(restriction)s[v1][1] - x%(restriction)s[v0][1];
const double det = std::sqrt(dx0*dx0 + dx1*dx1);
"""

# Code snippet for computing the determinant of the facet mapping in 2D
facet_determinant_3D = """\
// Vertices on faces
static unsigned int face_vertices[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};

// Get vertices
const unsigned int v0 = face_vertices[%(facet)s][0];
const unsigned int v1 = face_vertices[%(facet)s][1];
const unsigned int v2 = face_vertices[&(facet)s][1];

// Compute scale factor (area of face scaled by area of reference triangle)
const double a0 = (x%(restriction)s[v0][1]*x%(restriction)s[v1][2] + x%(restriction)s[v0][2]*x%(restriction)s[v2][1] + x%(restriction)s[v1][1]*x%(restriction)s[v2][2])
              - (x%(restriction)s[v2][1]*x%(restriction)s[v1][2] + x%(restriction)s[v2][2]*x%(restriction)s[v0][1] + x%(restriction)s[v1][1]*x%(restriction)s[v0][2]);
const double a1 = (x%(restriction)s[v0][2]*x%(restriction)s[v1][0] + x%(restriction)s[v0][0]*x%(restriction)s[v2][2] + x%(restriction)s[v1][2]*x%(restriction)s[v2][0])
              - (x%(restriction)s[v2][2]*x%(restriction)s[v1][0] + x%(restriction)s[v2][0]*x%(restriction)s[v0][2] + x%(restriction)s[v1][2]*x%(restriction)s[v0][0]);
const double a2 = (x%(restriction)s[v0][0]*x%(restriction)s[v1][1] + x%(restriction)s[v0][1]*x%(restriction)s[v2][0] + x%(restriction)s[v1][0]*x%(restriction)s[v2][1])
              - (x%(restriction)s[v2][0]*x%(restriction)s[v1][1] + x%(restriction)s[v2][1]*x%(restriction)s[v0][0] + x%(restriction)s[v1][0]*x%(restriction)s[v0][1]);
const double det = std::sqrt(a0*a0 + a1*a1 + a2*a2);
"""

# Code snippet for evaluate_dof in 2D
evaluate_dof_2D = """\
static double values[%d];
static double coordinates[2];

// Nodal coordinates on reference cell
static double X[%d][2] = %s;

// Components for each dof
static unsigned int components[%d] = %s;

// Extract vertex coordinates
const double * const * x = c.coordinates;

// Evaluate basis functions for affine mapping
const double w0 = 1.0 - X[i][0] - X[i][1];
const double w1 = X[i][0];
const double w2 = X[i][1];

// Compute affine mapping x = F(X)
coordinates[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0];
coordinates[1] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1];

// Evaluate function at coordinates
f.evaluate(values, coordinates, c);

// Pick component for evaluation
return values[components[i]];"""

# Code snippet for evaluate_dof in 3D
evaluate_dof_3D = """\
static double values[%d];
static double coordinates[3];

// Nodal coordinates on reference cell
static double X[%d][3] = %s;

// Components for each dof
static unsigned int components[%d] = %s;

// Extract vertex coordinates
const double * const * x = c.coordinates;

// Evaluate basis functions for affine mapping
const double w0 = 1.0 - X[i][0] - X[i][1] - X[i][2];
const double w1 = X[i][0];
const double w2 = X[i][1];
const double w3 = X[i][2];

// Compute affine mapping x = F(X)
coordinates[0] = w0*x[0][0] + w1*x[1][0] + w2*x[2][0] + w3*x[3][0];
coordinates[0] = w0*x[0][1] + w1*x[1][1] + w2*x[2][1] + w3*x[3][1];
coordinates[0] = w0*x[0][2] + w1*x[1][2] + w2*x[2][2] + w3*x[3][2];

// Evaluate function at coordinates
f.evaluate(values, coordinates, c);

// Pick component for evaluation
return values[components[i]];"""
