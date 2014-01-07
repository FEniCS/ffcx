# Code generation format strings for UFC (Unified Form-assembly Code) v. 2.3.0.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2014.

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
  virtual bool needs_mesh_entities(std::size_t d) const
  {
%(needs_mesh_entities)s
  }

  /// Return the topological dimension of the associated cell shape
  virtual std::size_t topological_dimension() const
  {
%(topological_dimension)s
  }

  /// Return the geometric dimension of the associated cell shape
  virtual std::size_t geometric_dimension() const
  {
%(geometric_dimension)s
  }

  /// Return the dimension of the global finite element function space
  virtual std::size_t global_dimension(const std::vector<std::size_t>&
                                       num_global_entities) const
  {
%(global_dimension)s
  }

  /// Return the dimension of the local finite element function space for a cell
  virtual std::size_t local_dimension() const
  {
%(local_dimension)s
  }

  /// Return the number of dofs on each cell facet
  virtual std::size_t num_facet_dofs() const
  {
%(num_facet_dofs)s
  }

  /// Return the number of dofs associated with each cell entity of dimension d
  virtual std::size_t num_entity_dofs(std::size_t d) const
  {
%(num_entity_dofs)s
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(std::size_t* dofs,
                             const std::vector<std::size_t>& num_global_entities,
                             const ufc::cell& c) const
  {
%(tabulate_dofs)s
  }

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(std::size_t* dofs,
                                   std::size_t facet) const
  {
%(tabulate_facet_dofs)s
  }

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  virtual void tabulate_entity_dofs(std::size_t* dofs,
                                    std::size_t d, std::size_t i) const
  {
%(tabulate_entity_dofs)s
  }

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double** dof_coordinates,
                                    const double* vertex_coordinates) const
  {
%(tabulate_coordinates)s
  }

  /// Return the number of sub dofmaps (for a mixed element)
  virtual std::size_t num_sub_dofmaps() const
  {
%(num_sub_dofmaps)s
  }

  /// Create a new dofmap for sub dofmap i (for a mixed element)
  virtual ufc::dofmap* create_sub_dofmap(std::size_t i) const
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
  virtual bool needs_mesh_entities(std::size_t d) const;

  /// Return the topological dimension of the associated cell shape
  virtual std::size_t topological_dimension() const;

  /// Return the geometric dimension of the associated cell shape
  virtual std::size_t geometric_dimension() const;

  /// Return the dimension of the global finite element function space
  virtual std::size_t global_dimension(const std::vector<std::size_t>&
                                       num_global_entities) const;

  /// Return the dimension of the local finite element function space for a cell
  virtual std::size_t local_dimension() const;

  /// Return the number of dofs on each cell facet
  virtual std::size_t num_facet_dofs() const;

  /// Return the number of dofs associated with each cell entity of dimension d
  virtual std::size_t num_entity_dofs(std::size_t d) const;

  /// Tabulate the local-to-global mapping of dofs on a cell
  virtual void tabulate_dofs(std::size_t* dofs,
                             const std::vector<std::size_t>& num_global_entities,
                             const ufc::cell& c) const;

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  virtual void tabulate_facet_dofs(std::size_t* dofs,
                                   std::size_t facet) const;

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  virtual void tabulate_entity_dofs(std::size_t* dofs,
                                    std::size_t d, std::size_t i) const;

  /// Tabulate the coordinates of all dofs on a cell
  virtual void tabulate_coordinates(double** coordinates,
                                    const double* vertex_coordinates) const;

  /// Return the number of sub dofmaps (for a mixed element)
  virtual std::size_t num_sub_dofmaps() const;

  /// Create a new dofmap for sub dofmap i (for a mixed element)
  virtual ufc::dofmap* create_sub_dofmap(std::size_t i) const;

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
bool %(classname)s::needs_mesh_entities(std::size_t d) const
{
%(needs_mesh_entities)s
}

/// Return the topological dimension of the associated cell shape
std::size_t %(classname)s::topological_dimension() const
{
%(topological_dimension)s
}

/// Return the geometric dimension of the associated cell shape
std::size_t %(classname)s::geometric_dimension() const
{
%(geometric_dimension)s
}

/// Return the dimension of the global finite element function space
std::size_t %(classname)s::global_dimension(const std::vector<std::size_t>&
                                            num_global_entities) const
{
%(global_dimension)s
}

/// Return the dimension of the local finite element function space for a cell
std::size_t %(classname)s::local_dimension() const
{
%(local_dimension)s
}

/// Return the number of dofs on each cell facet
std::size_t %(classname)s::num_facet_dofs() const
{
%(num_facet_dofs)s
}

/// Return the number of dofs associated with each cell entity of dimension d
std::size_t %(classname)s::num_entity_dofs(std::size_t d) const
{
%(num_entity_dofs)s
}

/// Tabulate the local-to-global mapping of dofs on a cell
void %(classname)s::tabulate_dofs(std::size_t* dofs,
                                  const std::vector<std::size_t>& num_global_entities,
                                  const ufc::cell& c) const
{
%(tabulate_dofs)s
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void %(classname)s::tabulate_facet_dofs(std::size_t* dofs,
                                        std::size_t facet) const
{
%(tabulate_facet_dofs)s
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void %(classname)s::tabulate_entity_dofs(std::size_t* dofs,
                                  std::size_t d, std::size_t i) const
{
%(tabulate_entity_dofs)s
}

/// Tabulate the coordinates of all dofs on a cell
void %(classname)s::tabulate_coordinates(double** dof_coordinates,
                                         const double* vertex_coordinates) const
{
%(tabulate_coordinates)s
}

/// Return the number of sub dofmaps (for a mixed element)
std::size_t %(classname)s::num_sub_dofmaps() const
{
%(num_sub_dofmaps)s
}

/// Create a new dofmap for sub dofmap i (for a mixed element)
ufc::dofmap* %(classname)s::create_sub_dofmap(std::size_t i) const
{
%(create_sub_dofmap)s
}

/// Create a new class instance
ufc::dofmap* %(classname)s::create() const
{
%(create)s
}
"""
