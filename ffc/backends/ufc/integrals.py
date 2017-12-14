# -*- coding: utf-8 -*-
# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2017

cell_integral_combined = """
class %(classname)s: public ufc::cell_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s) : ufc::cell_integral()%(initializer_list)s
  {
%(constructor)s
  }

  ~%(classname)s() override
  {
%(destructor)s
  }

  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       int cell_orientation) const final override
  {
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
  }

};
"""

cell_integral_header = """
class %(classname)s: public ufc::cell_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  const std::vector<bool> & enabled_coefficients() const final override;

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       int cell_orientation) const final override;

};
"""

cell_integral_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::cell_integral()%(initializer_list)s
{
%(constructor)s
}

%(classname)s::~%(classname)s()
{
%(destructor)s
}

const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
                                    int cell_orientation) const
{
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
}
"""

exterior_facet_integral_combined = """
class %(classname)s: public ufc::exterior_facet_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s) : ufc::exterior_facet_integral()%(initializer_list)s
  {
%(constructor)s
  }

  ~%(classname)s() override
  {
%(destructor)s
  }

  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t facet,
                       int cell_orientation) const final override
  {
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
  }

};
"""

exterior_facet_integral_header = """
class %(classname)s: public ufc::exterior_facet_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  const std::vector<bool> & enabled_coefficients() const final override;

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t facet,
                       int cell_orientation) const final override;

};
"""

exterior_facet_integral_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::exterior_facet_integral()%(initializer_list)s
{
%(constructor)s
}

%(classname)s::~%(classname)s()
{
%(destructor)s
}

const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
                                    std::size_t facet,
                                    int cell_orientation) const
{
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
}
"""

interior_facet_integral_combined = """
class %(classname)s: public ufc::interior_facet_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s) : ufc::interior_facet_integral()%(initializer_list)s
  {
%(constructor)s
  }

  ~%(classname)s() override
  {
%(destructor)s
  }

  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs_0,
                       const double * coordinate_dofs_1,
                       std::size_t facet_0,
                       std::size_t facet_1,
                       int cell_orientation_0,
                       int cell_orientation_1) const final override
  {
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
  }

};
"""

interior_facet_integral_header = """
class %(classname)s: public ufc::interior_facet_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  const std::vector<bool> & enabled_coefficients() const final override;

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs_0,
                       const double * coordinate_dofs_1,
                       std::size_t facet_0,
                       std::size_t facet_1,
                       int cell_orientation_0,
                       int cell_orientation_1) const final override;

};
"""

interior_facet_integral_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::interior_facet_integral()%(initializer_list)s
{
%(constructor)s
}

%(classname)s::~%(classname)s()
{
%(destructor)s
}

const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs_0,
                                    const double * coordinate_dofs_1,
                                    std::size_t facet_0,
                                    std::size_t facet_1,
                                    int cell_orientation_0,
                                    int cell_orientation_1) const
{
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
}
"""

vertex_integral_combined = """
class %(classname)s: public ufc::vertex_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s) : ufc::vertex_integral()%(initializer_list)s
  {
%(constructor)s
  }

  ~%(classname)s() override
  {
%(destructor)s
  }

  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t vertex,
                       int cell_orientation) const final override
  {
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
  }

};
"""

vertex_integral_header = """
class %(classname)s: public ufc::vertex_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  const std::vector<bool> & enabled_coefficients() const final override;

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t vertex,
                       int cell_orientation) const final override;

};
"""

vertex_integral_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::vertex_integral()%(initializer_list)s
{
%(constructor)s
}

%(classname)s::~%(classname)s()
{
%(destructor)s
}

const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
                                    std::size_t vertex,
                                    int cell_orientation) const
{
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
}
"""

custom_integral_combined = """\
class %(classname)s: public ufc::custom_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s) : ufc::custom_integral()%(initializer_list)s
  {
%(constructor)s
  }

  ~%(classname)s() override
  {
%(destructor)s
  }

  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  std::size_t num_cells() const final override
  {
%(num_cells)s
  }

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       const double * facet_normals,
                       int cell_orientation) const final override
  {
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
  }

};
"""

custom_integral_header = """
class %(classname)s: public ufc::custom_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  const std::vector<bool> & enabled_coefficients() const final override;

  std::size_t num_cells() const final override;

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       const double * facet_normals,
                       int cell_orientation) const final override;

};
"""

custom_integral_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::custom_integral()%(initializer_list)s
{
%(constructor)s
}

%(classname)s::~%(classname)s()
{
%(destructor)s
}

const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

std::size_t %(classname)s::num_cells() const
{
%(num_cells)s
}

void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
                                    std::size_t num_quadrature_points,
                                    const double * quadrature_points,
                                    const double * quadrature_weights,
                                    const double * facet_normals,
                                    int cell_orientation) const
{
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
}
"""

cutcell_integral_combined = """
class %(classname)s: public ufc::cutcell_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s) : ufc::cutcell_integral()%(initializer_list)s
  {
%(constructor)s
  }

  ~%(classname)s() override
  {
%(destructor)s
  }

  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       int cell_orientation) const final override
  {
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
  }

};
"""

cutcell_integral_header = """
class %(classname)s: public ufc::cutcell_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  const std::vector<bool> & enabled_coefficients() const final override;

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       int cell_orientation) const final override;

};
"""

cutcell_integral_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::cutcell_integral()%(initializer_list)s
{
%(constructor)s
}

%(classname)s::~%(classname)s()
{
%(destructor)s
}

const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
                                    std::size_t num_quadrature_points,
                                    const double * quadrature_points,
                                    const double * quadrature_weights,
                                    int cell_orientation) const
{
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
}
"""

interface_integral_combined = """
class %(classname)s: public ufc::interface_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s) : ufc::interface_integral()%(initializer_list)s
  {
%(constructor)s
  }

  ~%(classname)s() override
  {
%(destructor)s
  }

  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       const double * facet_normals,
                       int cell_orientation) const final override
  {
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
  }

};
"""

interface_integral_header = """
class %(classname)s: public ufc::interface_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  const std::vector<bool> & enabled_coefficients() const final override;

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       const double * facet_normals,
                       int cell_orientation) const final override;

};
"""

interface_integral_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::interface_integral()%(initializer_list)s
{
%(constructor)s
}

%(classname)s::~%(classname)s()
{
%(destructor)s
}

const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
                                    std::size_t num_quadrature_points,
                                    const double * quadrature_points,
                                    const double * quadrature_weights,
                                    const double * facet_normals,
                                    int cell_orientation) const
{
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
}
"""

overlap_integral_combined = """
class %(classname)s: public ufc::overlap_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s) : ufc::overlap_integral()%(initializer_list)s
  {
%(constructor)s
  }

  ~%(classname)s() override
  {
%(destructor)s
  }

  const std::vector<bool> & enabled_coefficients() const final override
  {
%(enabled_coefficients)s
  }

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       int cell_orientation) const final override
  {
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
  }

};
"""

overlap_integral_header = """
class %(classname)s: public ufc::overlap_integral
{%(members)s
public:

  %(classname)s(%(constructor_arguments)s);

  ~%(classname)s() override;

  const std::vector<bool> & enabled_coefficients() const final override;

  void tabulate_tensor(double * A,
                       const double * const * w,
                       const double * coordinate_dofs,
                       std::size_t num_quadrature_points,
                       const double * quadrature_points,
                       const double * quadrature_weights,
                       int cell_orientation) const final override;

};
"""

overlap_integral_implementation = """
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::overlap_integral()%(initializer_list)s
{
%(constructor)s
}

%(classname)s::~%(classname)s()
{
%(destructor)s
}

const std::vector<bool> & %(classname)s::enabled_coefficients() const
{
%(enabled_coefficients)s
}

void %(classname)s::tabulate_tensor(double * A,
                                    const double * const * w,
                                    const double * coordinate_dofs,
                                    std::size_t num_quadrature_points,
                                    const double * quadrature_points,
                                    const double * quadrature_weights,
                                    int cell_orientation) const
{
%(tabulate_tensor_comment)s
%(tabulate_tensor)s
}
"""
