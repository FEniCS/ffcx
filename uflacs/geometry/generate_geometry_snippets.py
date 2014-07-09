from uflacs.codeutils.format_code import Indented, WithKeywords
from six.moves import xrange

# A dependency graph like this might be a way to automatically figure out which quantities to generate?
dependencies = {
   'J': ('vertex_coordinates',),
   'detJ': ('J',),
   'K': ('J', 'detJ'),
   'x': ('xi', 'J', 'vertex_coordinates'),
   'xi': ('x', 'K', 'vertex_coordinates'),
  }

# FIXME: This is now duplicated between here and ffc backend, get from one place always
restriction_postfix = { "+": "_0", "-": "_1", None: "" }

#def missing_geometry_snippets(existing, wanted):
#    ...
#def generate_missing_geometry_snippets(existing, wanted, restriction):
#    ...

# TODO: Get names from some central place? Or is it good enough to just rely on convention?

def get_r_keywords(r):
    return { "r": restriction_postfix[r] }


### Generating generic code snippets such as matrix-multiply

def generate_array_definition_snippets(name, expressions, d, typename="const double"):
    "Generate combined definition and declaration of name[] = expressions[] dimension d."
    decl = ['%s %s[%d] = {' % (typename, name, d),
            Indented([e+',' for e in expressions[:-1]]
                     + [expressions[-1], '};'])]
    return decl

def generate_z_Axpy_snippets(name_z, name_A, name_x, name_y, zd, xd):
    fmt_A = dict(((i, j), '%s[%d*%d+%d]' % (name_A, i, xd, j)) for i in xrange(zd) for j in xrange(xd))
    fmt_x = ['%s[%d]' % (name_x, j) for j in xrange(xd)]
    fmt_y = ['%s[%d]' % (name_y, i) for i in xrange(zd)]
    fmt_Ax = [' + '.join('%s * %s' % (fmt_A[(i, j)], fmt_x[j]) for j in xrange(xd)) for i in xrange(zd)]
    fmt_z = ['%s + %s' % (fmt_Ax[i], fmt_y[i]) for i in xrange(zd)]
    return generate_array_definition_snippets(name_z, fmt_z, zd)

def generate_z_Axmy_snippets(name_z, name_A, name_x, name_y, zd, xd):
    "Generate combined definition and declaration of z = A (x - y) with dimensions zd,xd."
    fmt_A = dict(((i, j), '%s[%d*%d+%d]' % (name_A, i, xd, j)) for i in xrange(zd) for j in xrange(xd))
    fmt_xmy = ['(%s[%d] - %s[%d])' % (name_x, j, name_y, j) for j in xrange(xd)]
    fmt_z = [' + '.join('%s * %s' % (fmt_A[(i, j)], fmt_xmy[j]) for j in xrange(xd)) for i in xrange(zd)]
    return generate_array_definition_snippets(name_z, fmt_z, zd)


### Generating inline computations of geometry

def generate_x_from_xi_snippets(cell, restriction):
    "Generate combined definition and declaration of x = J xi + v."
    gd = cell.geometric_dimension()
    td = cell.topological_dimension()
    rk = get_r_keywords(restriction)
    name_A = "J%(r)s" % rk
    name_x = "xi%(r)s" % rk
    name_y = "vertex_coordinates%(r)s" % rk
    name_z = "x%(r)s" % rk
    return generate_z_Axpy_snippets(name_z, name_A, name_x, name_y, gd, td)

def generate_xi_from_x_snippets(cell, restriction):
    "Generate combined definition and declaration of xi = K (x - v)."
    gd = cell.geometric_dimension()
    td = cell.topological_dimension()
    rk = get_r_keywords(restriction)
    name_A = "K%(r)s" % rk
    name_x = "x%(r)s" % rk
    name_y = "vertex_coordinates%(r)s" % rk
    name_z = "xi%(r)s" % rk
    return generate_z_Axmy_snippets(name_z, name_A, name_x, name_y, td, gd)


### Generating calls to functions from ufc_geometry.h

_code = {
    "interval": {
        1: [],
        2: [],
        3: [],
        },
    "triangle": {
        2: [],
        3: [],
        },
    "tetrahedron": {
        3: [],
        }
    }

jacobian_code = {
    "interval": {
        1: ["double J%(r)s[1];",
            "compute_jacobian_interval_1d(J%(r)s, vertex_coordinates%(r)s);"],
        2: ["double J%(r)s[2];",
            "compute_jacobian_interval_2d(J%(r)s, vertex_coordinates%(r)s);"],
        3: ["double J%(r)s[3];",
            "compute_jacobian_interval_3d(J%(r)s, vertex_coordinates%(r)s);"],
        },
    "triangle": {
        2: ["double J%(r)s[4];",
            "compute_jacobian_triangle_2d(J%(r)s, vertex_coordinates%(r)s);"],
        3: ["double J%(r)s[6];",
            "compute_jacobian_triangle_3d(J%(r)s, vertex_coordinates%(r)s);"],
        },
    "tetrahedron": {
        3: ["double J%(r)s[9];",
            "compute_jacobian_tetrahedron_3d(J%(r)s, vertex_coordinates%(r)s);"],
        }
    }

jacobian_determinants_code = {
    "interval": {
        1: ["double detJ%(r)s;",
            "compute_jacobian_determinants_interval_1d(detJ%(r)s, J%(r)s);"],
        2: ["double det2%(r)s;", "double detJ%(r)s;",
            "compute_jacobian_determinants_interval_2d(det2%(r)s, detJ%(r)s, J%(r)s);"],
        3: ["double det2%(r)s;", "double detJ%(r)s;",
            "compute_jacobian_determinants_interval_3d(det2%(r)s, detJ%(r)s, J%(r)s);"],
        },
    "triangle": {
        2: ["double detJ%(r)s;",
            "compute_jacobian_determinants_triangle_2d(detJ%(r)s, J%(r)s);"],
        3: ["double den%(r)s;", "double det2%(r)s;", "double detJ%(r)s;", "double detc%(r)s[3];",
            "compute_jacobian_determinants_triangle_3d(den%(r)s, det2%(r)s, detJ%(r)s, detc%(r)s, J%(r)s);"],
        },
    "tetrahedron": {
        3: ["double detJ%(r)s;", "double detd%(r)s[9];",
            "compute_jacobian_determinants_tetrahedron_3d(detJ%(r)s, detd%(r)s, J%(r)s);"],
        }
    }

jacobian_inverse_code = {
    "interval": {
        1: ["double K%(r)s[1];",
            "new_compute_jacobian_inverse_interval_1d(K%(r)s, detJ%(r)s);",],
        2: ["double K%(r)s[2];",
            "new_compute_jacobian_inverse_interval_2d(K%(r)s, det2%(r)s, J%(r)s);"],
        3: ["double K%(r)s[3];",
            "new_compute_jacobian_inverse_interval_3d(K%(r)s, det2%(r)s, J%(r)s);"],
        },
    "triangle": {
        2: ["double K%(r)s[4];",
            "new_compute_jacobian_inverse_triangle_2d(K%(r)s, detJ%(r)s, J%(r)s);"],
        3: ["double K%(r)s[6];",
            "new_compute_jacobian_inverse_triangle_3d(K%(r)s, den%(r)s, detc%(r)s, J%(r)s);"],
        },
    "tetrahedron": {
        3: ["double K%(r)s[9];",
            "new_compute_jacobian_inverse_tetrahedron_3d(K%(r)s, detJ%(r)s, detd%(r)s);"],
        }
    }

cell_volume_code = {
    "interval": {
        1: ["double volume%(r)s = std::fabs(detJ%(r)s);"],
        2: ["double volume%(r)s = std::fabs(detJ%(r)s);"],
        3: ["double volume%(r)s = std::fabs(detJ%(r)s);"],
        },
    "triangle": {
        2: ["double volume%(r)s = std::fabs(detJ%(r)s) / 2.0;"],
        3: ["double volume%(r)s = std::fabs(detJ%(r)s) / 2.0;"],
        },
    "tetrahedron": {
        3: ["double volume%(r)s = std::fabs(detJ%(r)s) / 6.0;"],
        }
    }

circumradius_code = {
    "interval": {
        1: ["double circumradius%(r)s = volume%(r)s / 2.0;"],
        2: ["double circumradius%(r)s = volume%(r)s / 2.0;"],
        3: ["double circumradius%(r)s = volume%(r)s / 2.0;"],
        },
    "triangle": {
        2: ["double circumradius%(r)s;",
            "compute_circumradius_triangle_2d(circumradius%(r)s, vertex_coordinates%(r)s, J%(r)s, volume%(r)s);"],
        3: ["double circumradius%(r)s;",
            "compute_circumradius_triangle_3d(circumradius%(r)s, vertex_coordinates%(r)s, J%(r)s, volume%(r)s);"],
        },
    "tetrahedron": {
        3: ["double circumradius%(r)s;",
            "compute_circumradius_tetrahedron_3d(circumradius%(r)s, vertex_coordinates%(r)s, J%(r)s, volume%(r)s);"],
        }
    }

facet_scaling_code = {
    "interval": {
        1: [],
        2: [],
        3: [],
        },
    "triangle": {
        2: ["double dx%(r)s[2];",
            "compute_edge_scaling_factors_triangle_2d(dx%(r)s, vertex_coordinates%(r)s, facet%(r)s);",
            "double det%(r)s;",
            "compute_facet_scaling_factor_triangle_2d(det%(r)s, dx%(r)s);"],
        3: ["double dx%(r)s[3];",
            "compute_edge_scaling_factors_triangle_3d(dx%(r)s, vertex_coordinates%(r)s, facet%(r)s);",
            "double det%(r)s;",
            "compute_facet_scaling_factor_triangle_3d(det%(r)s, dx%(r)s);"],
        },
    "tetrahedron": {
        3: ["double a%(r)s[3];",
            "compute_face_scaling_factors_tetrahedron_3d(a%(r)s, vertex_coordinates%(r)s, facet%(r)s);",
            "double det%(r)s;",
            "compute_facet_scaling_factor_tetrahedron_3d(det%(r)s, a%(r)s);"],
        }
    }

facet_area_code = {
    "interval": {
        1: ["double facetarea = 1.0;"],
        2: ["double facetarea = 1.0;"],
        3: ["double facetarea = 1.0;"],
        },
    "triangle": {
        2: ["double facetarea = det%(r)s;"],
        3: ["double facetarea = det%(r)s;"],
        },
    "tetrahedron": {
        3: ["double facetarea = det%(r)s / 2.0;"], # 'det' is scaled by area of reference triangle
        }
    }

facet_direction_code = {
    "interval": {
        1: ["bool direction%(r)s;",
            "compute_facet_normal_direction_interval_1d(direction%(r)s, vertex_coordinates%(r)s, facet%(r)s);"],
        2: [],
        3: [],
        },
    "triangle": {
        2: ["bool direction%(r)s;",
            "compute_facet_normal_direction_triangle_2d(direction%(r)s, vertex_coordinates%(r)s, dx%(r)s, facet%(r)s);"],
        3: [],
        },
    "tetrahedron": {
        3: ["bool direction%(r)s;",
            "compute_facet_normal_direction_tetrahedron_3d(direction%(r)s, vertex_coordinates%(r)s, a%(r)s, facet%(r)s);"],
        }
    }

# FIXME: Compute with + restriction and flip sign of the - direction?
facet_normal_code = {
    "interval": {
        1: ["double n[1];",
            "compute_facet_normal_interval_1d(n, direction%(r)s);"],
        2: ["double n[2];",
            "compute_facet_normal_interval_2d(n, vertex_coordinates%(r)s, facet%(r)s);"],
        3: ["double n[3];",
            "compute_facet_normal_interval_3d(n, vertex_coordinates%(r)s, facet%(r)s);"],
        },
    "triangle": {
        2: ["double n[2];",
            "compute_facet_normal_triangle_2d(n, dx%(r)s, det%(r)s, direction%(r)s);"],
        3: ["double n[3];",
            "compute_facet_normal_triangle_3d(n, vertex_coordinates%(r)s, facet%(r)s);"],
        },
    "tetrahedron": {
        3: ["double n[3];",
            "compute_facet_normal_tetrahedron_3d(n, a%(r)s, det%(r)s, direction%(r)s);"],
        }
    }

def template_with_r_keywords(templates, cell, r):
    na = cell.cellname()
    gd = cell.geometric_dimension()
    code = templates[na][gd]
    return WithKeywords(code, get_r_keywords(r))

def generate_jacobian_snippets(cell, restriction):
    return template_with_r_keywords(jacobian_code, cell, restriction)

def generate_jacobian_determinants_snippets(cell, restriction):
    return template_with_r_keywords(jacobian_determinants_code, cell, restriction)

def generate_jacobian_inverse_snippets(cell, restriction):
    return template_with_r_keywords(jacobian_inverse_code, cell, restriction)

def generate_cell_scaling_factor_snippets(cell):
    code = ['double det = std::fabs(detJ);']
    return code

def generate_cell_volume_snippets(cell, restriction):
    return template_with_r_keywords(cell_volume_code, cell, restriction)

def generate_circumradius_snippets(cell, restriction):
    return template_with_r_keywords(circumradius_code, cell, restriction)

def generate_facet_scaling_factor_snippets(cell, restriction):
    return template_with_r_keywords(facet_scaling_code, cell, restriction)

def generate_facet_area_snippets(cell, restriction):
    return template_with_r_keywords(facet_area_code, cell, restriction)

def generate_facet_direction_snippets(cell, restriction):
    return template_with_r_keywords(facet_direction_code, cell, restriction)

def generate_facet_normal_snippets(cell, restriction):
    return template_with_r_keywords(facet_normal_code, cell, restriction)
