# Code generation format strings for UFC (Unified Form-assembly Code) v. 2.0.5.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2011.

dofmap_combined = """\
/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class %(classname)s: public ufc::dofmap
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s) : ufc::dofmap()%(initializer_list)s
  {
%(constructor)s
  }

  /// Destructor
  virtual ~%(classname)s()
  {
%(destructor)s
  }

  /// Return a string identifying the dofmap
  virtual const char* signature() const
  {
%(signature)s
  }

  /// Return true iff mesh entities of topological dimension d are needed
  virtual bool needs_mesh_entities(unsigned int d) const
  {
%(needs_mesh_entities)s
  }

  /// Initialize dofmap for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m)
  {
%(init_mesh)s
  }

  /// Initialize dofmap for given cell
  virtual void init_cell(const ufc::mesh& m,
                         const ufc::cell& c)
  {
%(init_cell)s
  }

  /// Finish initialization of dofmap for cells
  virtual void init_cell_finalize()
  {
%(init_cell_finalize)s
  }

  /// Return the topological dimension of the associated cell shape
  virtual unsigned int topological_dimension() const
  {
%(topological_dimension)s
  }

  /// Return the geometric dimension of the associated cell shape
  virtual unsigned int geometric_dimension() const
  {
%(geometric_dimension)s
  }

  /// Return the dimension of the global finite element function space
  virtual unsigned int global_dimension() const
  {
%(global_dimension)s
  }

  /// Return the dimension of the local finite element function space for a cell
  virtual unsigned int local_dimension(const ufc::cell& c) const
  {
%(local_dimension)s
  }

  /// Return the maximum dimension of the local finite element function space
  virtual unsigned int max_local_dimension() const
  {
%(max_local_dimension)s
  }

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const
  {
%(num_facet_dofs)s
  }

  /// Return the number of dofs associated with each cell entity of dimension d
  virtual unsigned int num_entity_dofs(unsigned int d) const
  {
%(num_entity_dofs)s
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
                                   unsigned int facet) const
  {
%(tabulate_facet_dofs)s
  }

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  virtual void tabulate_entity_dofs(unsigned int* dofs,
                                    unsigned int d, unsigned int i) const
  {
%(tabulate_entity_dofs)s
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double** coordinates,
                                    const ufc::cell& c) const
  {
%(tabulate_coordinates)s
  }

  /// Return the number of sub dofmaps (for a mixed element)
  virtual unsigned int num_sub_dofmaps() const
  {
%(num_sub_dofmaps)s
  }

  /// Create a new dofmap for sub dofmap i (for a mixed element)
  virtual ufc::dofmap* create_sub_dofmap(unsigned int i) const
  {
%(create_sub_dofmap)s
  }

  /// Create a new class instance
  virtual ufc::dofmap* create() const
  {
%(create)s
  }

};
"""

dofmap_header = """\
/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class %(classname)s: public ufc::dofmap
{%(members)s
public:

  /// Constructor
  %(classname)s(%(constructor_arguments)s);

  /// Destructor
  virtual ~%(classname)s();

  /// Return a string identifying the dofmap
  virtual const char* signature() const;

  /// Return true iff mesh entities of topological dimension d are needed
  virtual bool needs_mesh_entities(unsigned int d) const;

  /// Initialize dofmap for mesh (return true iff init_cell() is needed)
  virtual bool init_mesh(const ufc::mesh& m);

  /// Initialize dofmap for given cell
  virtual void init_cell(const ufc::mesh& m,
                         const ufc::cell& c);

  /// Finish initialization of dofmap for cells
  virtual void init_cell_finalize();

  /// Return the topological dimension of the associated cell shape
  virtual unsigned int topological_dimension() const;

  /// Return the geometric dimension of the associated cell shape
  virtual unsigned int geometric_dimension() const;

  /// Return the dimension of the global finite element function space
  virtual unsigned int global_dimension() const;

  /// Return the dimension of the local finite element function space for a cell
  virtual unsigned int local_dimension(const ufc::cell& c) const;

  /// Return the maximum dimension of the local finite element function space
  virtual unsigned int max_local_dimension() const;

  /// Return the number of dofs on each cell facet
  virtual unsigned int num_facet_dofs() const;

  /// Return the number of dofs associated with each cell entity of dimension d
  virtual unsigned int num_entity_dofs(unsigned int d) const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(unsigned int* dofs,
                             const ufc::mesh& m,
                             const ufc::cell& c) const;

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(unsigned int* dofs,
                                   unsigned int facet) const;

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  virtual void tabulate_entity_dofs(unsigned int* dofs,
                                    unsigned int d, unsigned int i) const;

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double** coordinates,
                                    const ufc::cell& c) const;

  /// Return the number of sub dofmaps (for a mixed element)
  virtual unsigned int num_sub_dofmaps() const;

  /// Create a new dofmap for sub dofmap i (for a mixed element)
  virtual ufc::dofmap* create_sub_dofmap(unsigned int i) const;

  /// Create a new class instance
  virtual ufc::dofmap* create() const;

};
"""

dofmap_implementation = """\
/// Constructor
%(classname)s::%(classname)s(%(constructor_arguments)s) : ufc::dofmap()%(initializer_list)s
{
%(constructor)s
}

/// Destructor
%(classname)s::~%(classname)s()
{
%(destructor)s
}

/// Return a string identifying the dofmap
const char* %(classname)s::signature() const
{
%(signature)s
}

/// Return true iff mesh entities of topological dimension d are needed
bool %(classname)s::needs_mesh_entities(unsigned int d) const
{
%(needs_mesh_entities)s
}

/// Initialize dofmap for mesh (return true iff init_cell() is needed)
bool %(classname)s::init_mesh(const ufc::mesh& m)
{
%(init_mesh)s
}

/// Initialize dofmap for given cell
void %(classname)s::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
%(init_cell)s
}

/// Finish initialization of dofmap for cells
void %(classname)s::init_cell_finalize()
{
%(init_cell_finalize)s
}

/// Return the topological dimension of the associated cell shape
unsigned int %(classname)s::topological_dimension() const
{
%(topological_dimension)s
}

/// Return the geometric dimension of the associated cell shape
unsigned int %(classname)s::geometric_dimension() const
{
%(geometric_dimension)s
}

/// Return the dimension of the global finite element function space
unsigned int %(classname)s::global_dimension() const
{
%(global_dimension)s
}

/// Return the dimension of the local finite element function space for a cell
unsigned int %(classname)s::local_dimension(const ufc::cell& c) const
{
%(local_dimension)s
}

/// Return the maximum dimension of the local finite element function space
unsigned int %(classname)s::max_local_dimension() const
{
%(max_local_dimension)s
}

/// Return the number of dofs on each cell facet
unsigned int %(classname)s::num_facet_dofs() const
{
%(num_facet_dofs)s
}

/// Return the number of dofs associated with each cell entity of dimension d
unsigned int %(classname)s::num_entity_dofs(unsigned int d) const
{
%(num_entity_dofs)s
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
                                        unsigned int facet) const
{
%(tabulate_facet_dofs)s
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void %(classname)s::tabulate_entity_dofs(unsigned int* dofs,
                                  unsigned int d, unsigned int i) const
{
%(tabulate_entity_dofs)s
}

/// Tabulate the coordinates of all dofs on a cell
void %(classname)s::tabulate_coordinates(double** coordinates,
                                         const ufc::cell& c) const
{
%(tabulate_coordinates)s
}

/// Return the number of sub dofmaps (for a mixed element)
unsigned int %(classname)s::num_sub_dofmaps() const
{
%(num_sub_dofmaps)s
}

/// Create a new dofmap for sub dofmap i (for a mixed element)
ufc::dofmap* %(classname)s::create_sub_dofmap(unsigned int i) const
{
%(create_sub_dofmap)s
}

/// Create a new class instance
ufc::dofmap* %(classname)s::create() const
{
%(create)s
}
"""
