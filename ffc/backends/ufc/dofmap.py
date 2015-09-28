# Code generation format strings for UFC (Unified Form-assembly Code) v. 1.7.0dev.
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2006-2015.

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
  ~%(classname)s() override
  {
%(destructor)s
  }

  /// Return a string identifying the dofmap
  const char * signature() const final override
  {
%(signature)s
  }

  /// Return true iff mesh entities of topological dimension d are needed
  bool needs_mesh_entities(std::size_t d) const final override
  {
%(needs_mesh_entities)s
  }

  /// Return the topological dimension of the associated cell shape
  std::size_t topological_dimension() const final override
  {
%(topological_dimension)s
  }

  /// Return the dimension of the global finite element function space
  std::size_t global_dimension(const std::vector<std::size_t>&
                               num_global_entities) const final override
  {
%(global_dimension)s
  }

  /// Return the dimension of the local finite element function space for a cell
  std::size_t num_element_dofs() const final override
  {
%(num_element_dofs)s
  }

  /// Return the number of dofs on each cell facet
  std::size_t num_facet_dofs() const final override
  {
%(num_facet_dofs)s
  }

  /// Return the number of dofs associated with each cell entity of dimension d
  std::size_t num_entity_dofs(std::size_t d) const final override
  {
%(num_entity_dofs)s
  }

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs(std::size_t * dofs,
                     const std::vector<std::size_t>& num_global_entities,
                     const std::vector<std::vector<std::size_t>>& entity_indices) const final override
  {
%(tabulate_dofs)s
  }

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  void tabulate_facet_dofs(std::size_t * dofs,
                           std::size_t facet) const final override
  {
%(tabulate_facet_dofs)s
  }

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  void tabulate_entity_dofs(std::size_t * dofs,
                            std::size_t d, std::size_t i) const final override
  {
%(tabulate_entity_dofs)s
  }


  /// Return the number of sub dofmaps (for a mixed element)
  std::size_t num_sub_dofmaps() const final override
  {
%(num_sub_dofmaps)s
  }

  /// Create a new dofmap for sub dofmap i (for a mixed element)
  ufc::dofmap * create_sub_dofmap(std::size_t i) const final override
  {
%(create_sub_dofmap)s
  }

  /// Create a new class instance
  ufc::dofmap * create() const final override
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
  ~%(classname)s() override;

  /// Return a string identifying the dofmap
  const char * signature() const final override;

  /// Return true iff mesh entities of topological dimension d are needed
  bool needs_mesh_entities(std::size_t d) const final override;

  /// Return the topological dimension of the associated cell shape
  std::size_t topological_dimension() const final override;

  /// Return the dimension of the global finite element function space
  std::size_t global_dimension(const std::vector<std::size_t>&
                               num_global_entities) const final override;

  /// Return the dimension of the local finite element function space for a cell
  std::size_t num_element_dofs() const final override;

  /// Return the number of dofs on each cell facet
  std::size_t num_facet_dofs() const final override;

  /// Return the number of dofs associated with each cell entity of dimension d
  std::size_t num_entity_dofs(std::size_t d) const final override;

  /// Tabulate the local-to-global mapping of dofs on a cell
  void tabulate_dofs(std::size_t * dofs,
                     const std::vector<std::size_t>& num_global_entities,
                     const std::vector<std::vector<std::size_t>>& entity_indices) const final override;

  /// Tabulate the local-to-local mapping from facet dofs to cell dofs
  void tabulate_facet_dofs(std::size_t * dofs,
                           std::size_t facet) const final override;

  /// Tabulate the local-to-local mapping of dofs on entity (d, i)
  void tabulate_entity_dofs(std::size_t * dofs,
                            std::size_t d, std::size_t i) const final override;

  /// Return the number of sub dofmaps (for a mixed element)
  std::size_t num_sub_dofmaps() const final override;

  /// Create a new dofmap for sub dofmap i (for a mixed element)
  ufc::dofmap * create_sub_dofmap(std::size_t i) const final override;

  /// Create a new class instance
  ufc::dofmap * create() const final override;

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
const char * %(classname)s::signature() const
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

/// Return the dimension of the global finite element function space
std::size_t %(classname)s::global_dimension(const std::vector<std::size_t>&
                                            num_global_entities) const
{
%(global_dimension)s
}

/// Return the dimension of the local finite element function space for a cell
std::size_t %(classname)s::num_element_dofs() const
{
%(num_element_dofs)s
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
void %(classname)s::tabulate_dofs(std::size_t * dofs,
                                  const std::vector<std::size_t>& num_global_entities,
                                  const std::vector<std::vector<std::size_t>>& entity_indices) const
{
%(tabulate_dofs)s
}

/// Tabulate the local-to-local mapping from facet dofs to cell dofs
void %(classname)s::tabulate_facet_dofs(std::size_t * dofs,
                                        std::size_t facet) const
{
%(tabulate_facet_dofs)s
}

/// Tabulate the local-to-local mapping of dofs on entity (d, i)
void %(classname)s::tabulate_entity_dofs(std::size_t * dofs,
                                  std::size_t d, std::size_t i) const
{
%(tabulate_entity_dofs)s
}

/// Return the number of sub dofmaps (for a mixed element)
std::size_t %(classname)s::num_sub_dofmaps() const
{
%(num_sub_dofmaps)s
}

/// Create a new dofmap for sub dofmap i (for a mixed element)
ufc::dofmap * %(classname)s::create_sub_dofmap(std::size_t i) const
{
%(create_sub_dofmap)s
}

/// Create a new class instance
ufc::dofmap * %(classname)s::create() const
{
%(create)s
}
"""
