"Headers for the UFC 1.0 format."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-08 -- 2007-01-08"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

class dof_map:

    constructor = """\
dof_map()"""

    signature = """\
const char* signature()"""

    needs_mesh_entities = """\
bool needs_mesh_entities(unsigned int d)"""
    
    init_mesh = """\
bool init_mesh(const mesh& mesh)"""

    init_cell = """\
void init_cell(const mesh& mesh,
               const cell& cell)"""

    global_dimension = """\
unsigned int global_dimension()"""

    local_dimension = """\
unsigned int local_dimension()"""

    tabulate_dofs = """\
void tabulate_dofs(unsigned int* dofs,
                   const mesh& m,
                   const cell& c)"""

    tabulate_facet_dofs = """\
void tabulate_facet_dofs(unsigned int* dofs,
                         const mesh& m,
                         const cell& c,
                         unsigned int facet)"""
    
class finite_element:

    constructor = """\
finite_element()"""

    signature = """\
const char* signature()"""

    cell_shape = """\
shape cell_shape()"""

    space_dimension = """\
unsigned int space_dimension()"""

    value_rank = """\
unsigned int value_rank()"""

    value_dimension = """\
unsigned int value_dimension(unsigned int i)"""

    evaluate_basis = """\
void evaluate_basis(double* values,
                    const double* x,
                    unsigned int i,
                    const cell& c)"""

    evaluate_dof = """\
double evaluate_dof(unsigned int i,
                    const function& f,
                    const cell& c)"""

    interpolate_vertex = """\
void interpolate_vertex_values(double* vertex_values,
                               const double* dof_values)"""

    num_sub_elements = """\
unsigned int num_sub_elements(unsigned int i)"""

    sub_element = """\
const finite_element& sub_element(unsigned int i)"""

class cell_integral:

  constructor = """\
cell_integral()"""

  tabulate_tensor = """\
void tabulate_tensor(double* A,
                     const double * const * w,
                     const cell& c)"""

class exterior_facet_integral:

    constructor = """\
exterior_facet_integral()"""

    tabulate_tensor = """\
void tabulate_tensor(double* A,
                     const double * const * w,
                     const cell& c,
                     unsigned int facet)"""

class interior_facet_integral:

    constructor = """\
interior_facet_integral()"""

    tabulate_tensor = """\
void tabulate_tensor(double* A,
                     const double * const * w,
                     const cell& c0,
                     const cell& c1,
                     unsigned int facet0,
                     unsigned int facet1)"""

class form:

    constructor = """\
form()"""

    signature = """\
const char* signature()"""

    rank = """\
unsigned int rank()"""

    num_coefficients = """\
unsigned int num_coefficients()"""

    create_cell_integral = """\
cell_integral* create_cell_integral()"""

    create_interior_facet_integral = """\
interior_facet_integral* create_interior_facet_integral()"""

    create_exterior_facet_integral = """\
exterior_facet_integral* create_exterior_facet_integral()"""

    create_dof_map = """\
dof_map* create_dof_map(unsigned int i)"""

    create_finite_element = """\
finite_element* create_finite_element(unsigned int i)"""
