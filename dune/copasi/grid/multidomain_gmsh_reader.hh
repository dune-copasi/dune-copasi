#ifndef DUNE_COPASI_GMSH_READER_HH
#define DUNE_COPASI_GMSH_READER_HH

#include <dune/copasi/concepts/grid.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/multidomaingrid.hh>

#include <dune/assembler/quantity_grid/grid.hh>

#include <dune/logging/logging.hh>

#include <dune/common/typetraits.hh>

#include <algorithm>
#include <type_traits>

namespace Dune {

template<class HostGrid, class MDGTraits>
class GridFactory<Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>> : public GridFactoryInterface<Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>> {
  
  using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;
  using Base = GridFactoryInterface<Grid>;
public:

  GridFactory()
    : _host_grid_factory{}
  {}

  // use default implementation from base class
  using Base::insertBoundarySegment;

  void insertVertex(const FieldVector<typename Grid::ctype, Grid::dimensionworld>& pos) override {
    _host_grid_factory.insertVertex(pos);
  }

  void insertElement(const GeometryType& type, const std::vector<unsigned int>& vertices) override {
    _host_grid_factory.insertElement(type, vertices);
  }

  void insertBoundarySegment(const std::vector<unsigned int>& vertices) override {
    _host_grid_factory.insertBoundarySegment(vertices);
  }

  std::unique_ptr<Grid> createGrid() override  {
    if (not _grid_ptr)
      makeGrid();
    assert(_grid_ptr);
    return std::exchange(_grid_ptr, nullptr);
  }

  void makeGrid(int max_subdomains = 1) {
    assert(not _grid_ptr);

    MDGTraits* traits;
    if constexpr (std::is_default_constructible_v<MDGTraits>)
      traits = new MDGTraits; // TODO assert max_subdomains!
    else
      traits = new MDGTraits(max_subdomains);
    
    _host_grid_ptr = _host_grid_factory.createGrid();
    _grid_ptr = std::make_unique<Grid>(_host_grid_ptr, *traits);
    delete traits;
  }

  Grid& grid() {
    assert(_grid_ptr);
    return *_grid_ptr;
  }

  HostGrid& hostGrid() {
    assert(_host_grid_ptr);
    return *_host_grid_ptr;
  }

  const std::shared_ptr<HostGrid>& hostGridPtr() {
    return _host_grid_ptr;
  }


  Dune::GridFactory<HostGrid>& hostGridFactory() {
    return _host_grid_factory;
  }

private:
  Dune::GridFactory<HostGrid> _host_grid_factory;
  std::shared_ptr<HostGrid> _host_grid_ptr;
  std::unique_ptr<Grid> _grid_ptr;
};


template<class HostGrid, class Quantity>
class GridFactory<Dune::Assembler::QuantityGrid<HostGrid, Quantity>> : public GridFactoryInterface<Dune::Assembler::QuantityGrid<HostGrid, Quantity>> {
  
  using Grid = Dune::Assembler::QuantityGrid<HostGrid, Quantity>;
  using Base = GridFactoryInterface<Dune::Assembler::QuantityGrid<HostGrid, Quantity>>;
public:

  // quantity_cast define what is the quantitative interpretation of positions instantiated with dimensionless units (i.e. double)
  // by default, it is assumed that they have the same units as the grid.
  GridFactory(
    Dune::GridFactory<HostGrid>&& host_grid_factory = {},
    std::function<Quantity(double)> quantity_cast = [](double v){return Quantity{v};}
  ) : _host_grid_factory{std::move(host_grid_factory)}
    , _quantity_cast{quantity_cast}
  {}

  // use default implementation from base class
  using Base::insertBoundarySegment;

  void insertVertex(const FieldVector<double, Grid::dimensionworld>& dimensionless_pos) {
    // TODO: switch units from options
    FieldVector<Quantity, Grid::dimensionworld> pos;
    for (std::size_t i = 0; i != Grid::dimensionworld; ++i)
      pos[i] = _quantity_cast(dimensionless_pos[i]);
    this->insertVertex(pos);
  }


  void insertVertex(const FieldVector<Quantity, Grid::dimensionworld>& pos) override {
    auto host_pos = std::bit_cast<FieldVector<typename HostGrid::ctype, Grid::dimensionworld>>(pos);
    _host_grid_factory.insertVertex(host_pos);
  }

  void insertElement(const GeometryType& type, const std::vector<unsigned int>& vertices) override {
    _host_grid_factory.insertElement(type, vertices);
  }

  void insertBoundarySegment(const std::vector<unsigned int>& vertices) override {
    _host_grid_factory.insertBoundarySegment(vertices);
  }

  std::unique_ptr<Grid> createGrid() override  {
    return std::make_unique<Grid>(_host_grid_factory.createGrid());
  }

private:
  GridFactory<HostGrid> _host_grid_factory;
  std::function<Quantity(double)> _quantity_cast;
};

// template<class HostGrid, class Quantity>
// class GridFactory<Dune::Assembler::QuantityGrid<HostGrid, Quantity>> : public GridFactoryInterface<Dune::Assembler::QuantityGrid<HostGrid, Quantity>>{

// };


template<class HostGrid, class MDGTraits>
class GmshReaderParser<Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>> : public Dune::GmshReaderParser<HostGrid> {
public:

  using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

  GmshReaderParser(Dune::GridFactory<Grid>& factory, bool v, bool i)
    : Dune::GmshReaderParser<HostGrid>{factory.hostGridFactory(), v, i}
    , _factory{factory}
  {}
 

  void read (const std::string& f)
  {
    Dune::GmshReaderParser<HostGrid>::read(f);
    const auto& index_map = this->elementIndexMap();
    int max_subdomains = *std::max_element(begin(index_map), end(index_map));
    _factory.makeGrid(max_subdomains);
    auto& grid = _factory.grid();
    
    grid.startSubDomainMarking();
    unsigned int i = 0;
    for (const auto& cell : elements(grid.leafGridView())) {
      auto subdomain = index_map[i] - 1; // gmsh index starts at 1!
      grid.addToSubDomain(subdomain, cell), i++;
    }

    grid.preUpdateSubDomains();
    grid.updateSubDomains();
    grid.postUpdateSubDomains();
  }

private:
  Dune::GridFactory<Grid>& _factory;
};

}


namespace Dune::Copasi {
/**
 * @brief      This class describes a multi domain gmsh reader.
 *
 * @tparam     Grid  The grid
 */
template<class Grid>
class MultiDomainGmshReader
{
  static_assert(Dune::Copasi::Concept::isMultiDomainGrid<Grid>(),
                "Grid is not multidomain grid");

  // this is due the fact that multidomain grids do not export its traits type.
  // And, since it is needed for the some grids construction, the only way is by
  // template specialization.
  static_assert(AlwaysFalse<Grid>::value,
                "Not implemented: This class only accepts multidomain grids "
                "from the dune-multidomain grid module");
};

/**
 * @brief      This class describes a multi domain gmsh reader.
 *
 * @tparam     HostGrid   The host grid of the multidomain grid
 * @tparam     MDGTraits  Traits of the multidomain grid (this is different to
 *                        the grid interface traits!)
 */
template<class HostGrid, class MDGTraits>
class [[deprecated]] MultiDomainGmshReader<Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>>
{
public:
  using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

  /**
   * @brief      Reads a grid out of a file name and a configuration file
   *
   * @param[in]  fileName                The gmsh file name
   * @param[in]  insertBoundarySegments  Bool to include boundary segments while
   *                                     reading file
   *
   * @return     A pair containing pointers to the host and the multidomain grid
   */
  static auto read(const std::string& fileName,
                   bool insertBoundarySegments = true)
  {
    Dune::GridFactory<Grid> factory;
    Dune::GmshReader<Grid>::read(factory, fileName, false, insertBoundarySegments);
    std::shared_ptr<Grid> grid_ptr = factory.createGrid();
    return std::make_pair(grid_ptr, factory.hostGridPtr());
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GMSH_READER_HH
