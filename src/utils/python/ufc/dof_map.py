
dof_map_combined = """\
/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class %(classname)s: public ufc::dof_map
{
%(members)s

public:

  /// Constructor
  %(classname)s()
  {
%(constructor)s
  }

  /// Destructor
  ~%(classname)s()
  {
%(destructor)s
  }

  /// Return a string identifying the dof map
  const char* signature() const
  {
%(signature)s
  }

  /// Return true iff mesh entities of topological dimension d are needed
  bool needs_mesh_entities(unsigned int d) const
  {
%(needs_mesh_entities)s
  }

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  bool init_mesh(const ufc::mesh& mesh)
  {
%(init_mesh)s
  }

  /// Initialize dof map for given cell
  void init_cell(const ufc::mesh& m,
                 const ufc::cell& c)
  {
%(init_cell)s
  }

  /// Return the dimension of the global finite element function space
  unsigned int global_dimension() const
  {
%(global_dimension)s
  }

  /// Return the dimension of the local finite element function space
  unsigned int local_dimension() const
  {
%(local_dimension)s
  }

  /// Return the number of dofs on a facets of a cell
  unsigned int num_facet_dofs() const
  {
%(num_facet_dofs)s
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs(unsigned int* dofs,
                     const ufc::mesh& m,
                     const ufc::cell& c) const
  {
%(tabulate_dofs)s
  }

  /// Tabulate the local-to-global mapping of dofs on a facet of a cell
  void tabulate_facet_dofs(unsigned int* dofs,
                           const ufc::mesh& m,
                           const ufc::cell& c,
                           unsigned int facet) const
  {
%(tabulate_facet_dofs)s
  }

};
"""


dof_map_header = """\
/// This class defines the interface for a local-to-global mapping of
/// degrees of freedom (dofs).

class %(classname)s: public ufc::dof_map
{
%(members)s
public:

  /// Constructor
  %(classname)s();

  /// Destructor
  ~%(classname)s();

  /// Return a string identifying the dof map
  const char* signature() const;

  /// Return true iff mesh entities of topological dimension d are needed
  bool needs_mesh_entities(unsigned int d) const;

  /// Initialize dof map for mesh (return true iff init_cell() is needed)
  bool init_mesh(const ufc::mesh& mesh);

  /// Initialize dof map for given cell
  void init_cell(const ufc::mesh& m,
                 const ufc::cell& c);

  /// Return the dimension of the global finite element function space
  unsigned int global_dimension() const;

  /// Return the dimension of the local finite element function space
  unsigned int local_dimension() const;

  /// Return the number of dofs on a facets of a cell
  unsigned int num_facet_dofs() const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs(unsigned int* dofs,
                     const ufc::mesh& m,
                     const ufc::cell& c) const;

  /// Tabulate the local-to-global mapping of dofs on a facet of a cell
  void tabulate_facet_dofs(unsigned int* dofs,
                           const ufc::mesh& m,
                           const ufc::cell& c,
                           unsigned int facet) const;

};
"""



dof_map_implementation = """\
/// Constructor
%(classname)s::%(classname)s()
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
bool %(classname)s::init_mesh(const ufc::mesh& mesh)
{
%(init_mesh)s
}

/// Initialize dof map for given cell
void %(classname)s::init_cell(const ufc::mesh& m,
                              const ufc::cell& c)
{
%(init_cell)s
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

/// Return the number of dofs on a facets of a cell
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

/// Tabulate the local-to-global mapping of dofs on a facet of a cell
void %(classname)s::tabulate_facet_dofs(unsigned int* dofs,
                                        const ufc::mesh& m,
                                        const ufc::cell& c,
                                        unsigned int facet) const
{
%(tabulate_facet_dofs)s
}
"""
