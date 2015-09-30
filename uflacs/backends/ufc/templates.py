
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


domain_template = """\
/// A representation of a geometric domain parameterized by a local basis on each cell
class %(classname)s
{%(preamble)s

  /// Return domain signature string
  %(pre)sconst char * signature() const %(post)s
%(signature)s
  /// Create object of the same type
  %(pre)sdomain * create() const %(post)s
%(create)s
  /// Return geometric dimension of the domain
  %(pre)sstd::size_t geometric_dimension() const %(post)s
%(geometric_dimension)s
  /// Return topological dimension of the domain
  %(pre)sstd::size_t topological_dimension() const %(post)s
%(topological_dimension)s
  /// Return cell shape of the domain
  %(pre)sshape cell_shape() const %(post)s
%(cell_shape)s
  /// Create finite_element object representing the coordinate parameterization
  %(pre)sfinite_element * create_coordinate_finite_element() const %(post)s
%(create_coordinate_finite_element)s
  /// Create dofmap object representing the coordinate parameterization
  %(pre)sdofmap * create_coordinate_dofmap() const %(post)s
%(create_coordinate_dofmap)s
  /// Compute physical coordinates x from reference coordinates X, the inverse of compute_reference_coordinates
  %(pre)svoid compute_physical_coordinates(
      double * x, std::size_t num_points,
      const double * X, const double * K,
      const double * coordinate_dofs, int cell_orientation) const %(post)s
%(compute_physical_coordinates)s
  /// Compute reference coordinates X from physical coordinates x, the inverse of compute_physical_coordinates
  %(pre)svoid compute_reference_coordinates(
      double * X, std::size_t num_points,
      const double * x,
      const double * coordinate_dofs, int cell_orientation) const %(post)s
%(compute_reference_coordinates)s
  /// Compute Jacobian of coordinate mapping J = dx/dX at reference coordinates X
  %(pre)svoid compute_jacobians(
      double * J, std::size_t num_points,
      const double * X,
      const double * coordinate_dofs, int cell_orientation) const %(post)s
%(compute_jacobians)s
  /// Compute determinants of (pseudo-)Jacobians J
  %(pre)svoid compute_jacobian_determinants(
      double * detJ, std::size_t num_points,
      const double * J) const %(post)s
%(compute_jacobian_determinants)s
  /// Compute (pseudo-)inverses K of (pseudo-)Jacobians J
  %(pre)svoid compute_jacobian_inverses(
      double * K, std::size_t num_points,
      const double * J, const double * detJ) const %(post)s
%(compute_jacobian_inverses)s
  /// Combined (for convenience) computation of x, J, detJ, K from X and coordinate_dofs on a cell
  %(pre)svoid compute_geometry(
      double * x, double * J, double * detJ, double * K, std::size_t num_points,
      const double * X,
      const double * coordinate_dofs, int cell_orientation) const %(post)s
%(compute_geometry)s
};
"""


factory_template = """\
class %(namespace)s%(classname)s;

extern "C" %(namespace)s%(classname)s * create_%(classname)s()
{
  return new %(namespace)s%(classname)s();
}
"""


def extract_keywords(template):
    r = re.compile(r"%\(([a-zA-Z0-9_]*)\)")
    return set(r.findall(template))

def carry_key(key):
    return "%%(%s)s" % (key,)

def updated(src, **kwargs):
    r = src.copy()
    r.update(kwargs)
    return r

def generate_ihi(template, classname):
    empty_keys = dict.fromkeys(extract_keywords(template), "")
    carry_keys = {key: carry_key(key) for key in empty_keys}
    assert template == (template % carry_keys)

    interface_keys = updated(empty_keys,
        pre="virtual ",
        post="= 0;",
        classname=classname,
        preamble="\npublic:\n  virtual ~%s() {}" % (classname,),
        )

    implementation_keys = dict(
        pre="",
        post="final override",
        classname=carry_key("classname"),
        preamble="%s\npublic:" % (carry_key("preamble"),),
        )

    function_keys = set(empty_keys) - set(implementation_keys)

    implementation_keys.update({ key: "  {\n%%(%s)s\n  }\n" % (key,) for key in function_keys })

    interface      = template % interface_keys
    implementation = template % implementation_keys

    return interface, implementation


domain_interface, domain_implementation = generate_ihi(domain_template, classname="domain")


interface_prefix = """\
#ifndef _UFC_H_
#define _UFC_H_
/**
 * Licence header.
 */
"""

interface_postfix = """\
#endif
"""

"""
ufc_interface = "\n\n".join([
    interface_prefix,
    cell_interface,
    function_interface,
    domain_interface,
    finite_element_interface,
    dofmap_interface,
    form_interface,
    interface_postfix,
    ])
"""
