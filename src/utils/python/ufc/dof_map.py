# Code generation format strings for UFC (Unified Form-assembly Code) v. 1.0.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenics.org/) 2006-2007.

dof_map_combined = """\
/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class %(classname)s: public ufc::dof_map
{%(members)s
public:

  /// Constructor
  %(classname)s() : ufc::dof_map()
  {
%(constructor)s
  }

  /// Destructor
  virtual ~%(classname)s()
  {
%(destructor)s
  }

  /// Return a string identifying the dof map
  virtual const char* signature() const
  {
%(signature)s
  }

  /// Return true iff mesh entities of topological dimension d are needed
  virtual bool needs_mesh_entities(unsigned int d) const
  {
%(needs_mesh_entities)s
  }

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m)
  {
%(init_mesh)s
  }

  /// Initialize dof map for given cell
  virtual void init_cell(const ufc::mesh& m,
                         const ufc::cell& c)
  {
%(init_cell)s
  }

  /// Finish initialization of dof map for cells
  virtual void init_cell_finalize()
  {
%(init_cell_finalize)s
  }

  /// Return the dimension of the global finite element function space
  virtual unsigned int global_dimension() const
  {
%(global_dimension)s
  }

  /// Return the dimension of the local finite element function space
  virtual unsigned int local_dimension() const
  {
%(local_dimension)s
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
%(num_facet_dofs)s
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(unsigned int* dofs,
                             const ufc::mesh& m,
                             const ufc::cell& c) const
  {
%(tabulate_dofs)s
  }

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(unsigned int* dofs,
                                   const ufc::mesh& m,
                                   const ufc::cell& c,
                                   unsigned int facet) const
  {
%(tabulate_facet_dofs)s
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double **coordinates,
                                    const ufc::cell& c) const
  {
%(tabulate_coordinates)s
  }

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const
  {
%(num_sub_dof_maps)s
  }

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const
  {
%(create_sub_dof_map)s
  }

};
"""

dof_map_header = """\
/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class %(classname)s: public ufc::dof_map
{%(members)s
public:

  /// Constructor
  %(classname)s();

  /// Destructor
  virtual ~%(classname)s();

  /// Return a string identifying the dof map
  virtual const char* signature() const;

  /// Return true iff mesh entities of topological dimension d are needed
  virtual bool needs_mesh_entities(unsigned int d) const;

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m);

  /// Initialize dof map for given cell
  virtual void init_cell(const ufc::mesh& m,
                         const ufc::cell& c);

  /// Finish initialization of dof map for cells
  virtual void init_cell_finalize();

  /// Return the dimension of the global finite element function space
  virtual unsigned int global_dimension() const;

  /// Return the dimension of the local finite element function space
  virtual unsigned int local_dimension() const;

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(unsigned int* dofs,
                             const ufc::mesh& m,
                             const ufc::cell& c) const;

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(unsigned int* dofs,
                                   const ufc::mesh& m,
                                   const ufc::cell& c,
                                   unsigned int facet) const;

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double **coordinates,
                                    const ufc::cell& c) const;

  /// Return the number of sub dof maps (for a mixed element)
  virtual unsigned int num_sub_dof_maps() const;

  /// Create a new dof_map for sub dof map i (for a mixed element)
  virtual ufc::dof_map* create_sub_dof_map(unsigned int i) const;

};
"""

dof_map_implementation = """\
/// Constructor
%(classname)s::%(classname)s() : ufc::dof_map()
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s()
{
%(destructor)s
}

/// Return a string identifying the dof map
const char* %(classname)s::signature() const
{
%(signature)s
}

/// Return true iff mesh entities of topological dimension d are needed
bool %(classname)s::needs_mesh_entities(unsigned int d) const
{
%(needs_mesh_entities)s
}

/// Initialize dof map for mesh (return true iff init_cell() is needed)
bool %(classname)s::init_mesh(const ufc::mesh& m)
{
%(init_mesh)s
}

/// Initialize dof map for given cell
void %(classname)s::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
%(init_cell)s
}

/// Finish initialization of dof map for cells
void %(classname)s::init_cell_finalize()
{
%(init_cell_finalize)s
}

/// Return the dimension of the global finite element function space
unsigned int %(classname)s::global_dimension() const
{
%(global_dimension)s
}

/// Return the dimension of the local finite element function space
unsigned int %(classname)s::local_dimension() const
{
%(local_dimension)s
}

/// Return the number of dofs on each cell facet
unsigned int %(classname)s::num_facet_dofs() const
{
%(num_facet_dofs)s
}

/// Tabulate the local-to-global mapping of dofs on a cell
void %(classname)s::tabulate_dofs(unsigned int* dofs,
                                  const ufc::mesh& m,
                                  const ufc::cell& c) const
{
%(tabulate_dofs)s
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void %(classname)s::tabulate_facet_dofs(unsigned int* dofs,
                                        const ufc::mesh& m,
                                        const ufc::cell& c,
                                        unsigned int facet) const
{
%(tabulate_facet_dofs)s
}

/// Tabulate the coordinates of all dofs on a cell
void %(classname)s::tabulate_coordinates(double **coordinates,
                                         const ufc::cell& c) const
{
%(tabulate_coordinates)s
}

/// Return the number of sub dof maps (for a mixed element)
unsigned int %(classname)s::num_sub_dof_maps() const
{
%(num_sub_dof_maps)s
}

/// Create a new dof_map for sub dof map i (for a mixed element)
ufc::dof_map* %(classname)s::create_sub_dof_map(unsigned int i) const
{
%(create_sub_dof_map)s
}

"""
