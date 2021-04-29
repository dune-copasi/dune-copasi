#ifndef DUNE_COPASI_GMSH_READER_HH
#define DUNE_COPASI_GMSH_READER_HH

#include <dune/copasi/concepts/grid.hh>

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/multidomaingrid.hh>

#include <dune/logging/logging.hh>

#include <dune/common/typetraits.hh>

#include <algorithm>
#include <type_traits>

namespace Dune::Copasi {


template<class Container, class Mapper>
struct ContainerDataHandle
  : public Dune::CommDataHandleIF<ContainerDataHandle<Container,Mapper>,
                                  typename Container::value_type
                                  >
{
  ContainerDataHandle(Container& container, const Mapper& mapper)
    : _grid(mapper.gridView().grid())
    , _container(container)
    , _mapper(mapper)
  {
    // auto mapper = MultipleCodimMultipleGeomTypeMapper{grid_view};
    auto& id_set = _grid.localIdSet();
    typename Mapper::Index index;
    for (const auto& e : elements(_mapper.gridView()))
      if (_mapper.contains(e,index))
        _id_map[id_set.id(e)] = _container[index];
  }

  bool contains (int dim, int codim) const
  {
    assert(dim == Grid::dimension);
    return _grid.comm().max(_mapper.types(codim).size());
  }

  bool fixedsize(int dim, int codim) const
  {
    // TODO: fix this
    return true;
  }

  template<typename Entity>
  std::size_t size(const Entity& e) const
  {
    return _mapper.size(e.type());
  }

  template<typename MessageBufferImp, typename Entity>
  void gather(MessageBufferImp& buf, const Entity& e) const
  {
    auto it = _id_map.find(_grid.localIdSet().id(e));
    assert(it != _id_map.end());
    buf.write(it->second);
  }

  template<typename MessageBufferImp, typename Entity>
  void scatter(MessageBufferImp& buf, const Entity& e, std::size_t n)
  {
    auto id = _grid.localIdSet().id(e);
    ValueType tmp;
    buf.read(tmp);
    _id_map.insert({id,tmp});
  }

  void update()
  {
    _mapper.update();
    auto& id_set = _grid.localIdSet();
    typename Mapper::Index index;
    _container.resize(_mapper.size());
    for (const auto& e : elements(_mapper.gridView()))
      if (_mapper.contains(e,index))
      {
        auto it = _id_map.find(id_set.id(e));
        assert(it != _id_map.end());
        _container[index] = it->second;
      }
  }

  using Grid = typename Mapper::GridView::Grid;
  using IdType = typename Grid::LocalIdSet::IdType;
  using ValueType = typename Container::value_type;

  const Grid& _grid;
  Container& _container;
  Mapper _mapper;
  std::map<IdType,ValueType> _id_map;
  std::map<std::array<std::size_t,2>,std::size_t> _sizes;
};


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
class MultiDomainGmshReader<Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>>
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
    auto logger = Dune::Logging::Logging::componentLogger({}, "grid");
    const bool cout_redirected = Dune::Logging::Logging::isCoutRedirected();
    const bool verbose = logger.level() > Logging::LogLevel::detail;

    if (not cout_redirected)
      Dune::Logging::Logging::redirectCout(logger.name(),
                                           Logging::LogLevel::detail);

    // make a grid factory
    Dune::GridFactory<HostGrid> factory;

    // create parse object
    GmshReaderParser<HostGrid> parser(factory, verbose, insertBoundarySegments);
    logger.info("Reading grid file: '{}'"_fmt,fileName);
    parser.read(fileName);

    auto index_map = parser.elementIndexMap();

    std::shared_ptr<HostGrid> host_grid = factory.createGrid();

    auto mapper = LeafMultipleCodimMultipleGeomTypeMapper{*host_grid,mcmgElementLayout()};
    if (host_grid->comm().size() > 0)
    {
      auto data_handle = ContainerDataHandle{index_map,mapper};
      logger.debug("Load balance grid"_fmt);
      host_grid->loadBalance(data_handle);
      data_handle.update();
    }

    // auto writer = VTKWriter{host_grid->leafGridView()};
    // writer.addCellData(index_map,"subdomain");
    // writer.write("test_multidomain_gmsh_reader");

    auto [min_subdomain_it,max_subdomain_it] = std::minmax_element(index_map.begin(), index_map.end());
    auto [min_subdomain,max_subdomain] = std::tie(*min_subdomain_it,*max_subdomain_it);
    // Subdoman GMSH indices should start at one
    assert(max_subdomain > 0);
    MDGTraits* traits;
    if constexpr (std::is_default_constructible_v<MDGTraits>)
      traits = new MDGTraits;
    else
      traits = new MDGTraits(max_subdomain);

    logger.detail("Creating multidomain grid from host grid"_fmt);
    std::shared_ptr<Grid> grid = std::make_shared<Grid>(*host_grid, *traits);
    delete traits;

    logger.trace("Creating multidomain sub-domains"_fmt);
    grid->startSubDomainMarking();

    for (const auto& cell : elements(mapper.gridView(),Partitions::interior)) {
      typename std::decay_t<decltype(mapper)>::Index index;
      if (mapper.contains(cell,index)) {
        // gmsh index starts at 1!
        auto subdomain = index_map[ index ] - 1;
        logger.trace( 1, "Cell {} added to subdomain {}"_fmt, cell.type(), subdomain);
        grid->addToSubDomain(subdomain, grid->wrapHostEntity(cell));
      }
    }

    grid->preUpdateSubDomains();
    grid->updateSubDomains();
    grid->postUpdateSubDomains();

    if (logger.level() >= Logging::LogLevel::detail) {
      logger.detail("Grid info"_fmt);
      gridinfo(*grid,"  ");
      for (int i = 0; i < max_subdomain; ++i) {
        logger.detail("Subdomain {} info"_fmt, i);
        gridinfo(grid->subDomain(i), "  ");

        // auto writer = VTKWriter{grid->subDomain(i).leafGridView()};
        // std::vector<int> values(grid->subDomain(i).leafGridView().size(0),0);
        // writer.addCellData(values,"subdomain");
        // writer.write("test_multidomain_gmsh_reader_" + std::to_string(i));
      }
    }

    if (not cout_redirected)
      Dune::Logging::Logging::restoreCout();

    return std::make_pair(grid, host_grid);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GMSH_READER_HH
