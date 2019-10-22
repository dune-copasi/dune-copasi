#ifndef DUNE_COPASI_GMSH_READER_HH
#define DUNE_COPASI_GMSH_READER_HH

#include <dune/copasi/concepts/grid.hh>

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/multidomaingrid.hh>

#include <dune/logging/logging.hh>

#include <dune/common/typetraits.hh>

#include <algorithm>
#include <type_traits>

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
class MultiDomainGmshReader<Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>>
{
public:
  using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

  /**
   * @brief      Reads a grid out of a file name and a configuration file
   *
   * @param[in]  fileName                The gmsh file name
   * @param[in]  config                  The configuration file
   * @param[in]  insertBoundarySegments  Bool to include boundary segments while
   *                                     reading file
   *
   * @return     A pair containing pointers to the host and the multidomain grid
   */
  static auto read(const std::string& fileName,
                   ParameterTree config,
                   bool insertBoundarySegments = true)
  {
    Logging::Logger _logger =
      Logging::Logging::componentLogger(config, "default");
    const bool cout_redirected = Dune::Logging::Logging::isCoutRedirected();
    const bool verbose = _logger.level() > Logging::LogLevel::info;

    if (not cout_redirected)
      Dune::Logging::Logging::redirectCout(_logger.name());

    // make a grid factory
    Dune::GridFactory<HostGrid> factory;

    // create parse object
    GmshReaderParser<HostGrid> parser(factory, verbose, insertBoundarySegments);
    parser.read(fileName);

    auto index_map = parser.elementIndexMap();

    std::shared_ptr<HostGrid> host_grid = factory.createGrid();

    int max_subdomains = *std::max_element(index_map.begin(), index_map.end());
    MDGTraits* traits;
    if constexpr (std::is_default_constructible_v<MDGTraits>)
      traits = new MDGTraits;
    else
      traits = new MDGTraits(max_subdomains);

    std::shared_ptr<Grid> grid = std::make_shared<Grid>(*host_grid, *traits);
    delete traits;

    grid->startSubDomainMarking();
    unsigned int i = 0;
    for (const auto& cell : elements(grid->leafGridView())) {
      // gmsh index starts at 1!
      grid->addToSubDomain(index_map[i] - 1, cell), i++;
    }

    grid->preUpdateSubDomains();
    grid->updateSubDomains();
    grid->postUpdateSubDomains();

    // _logger.debug("Load balance grid"_fmt);
    // grid->loadBalance();

    if (_logger.level() >= Logging::LogLevel::trace) {
      _logger.debug("Grid info"_fmt);
      gridinfo(*grid);
      for (int i = 0; i < max_subdomains; ++i) {
        _logger.debug("Subdomain {} info"_fmt, i);
        gridinfo(grid->subDomain(i));
      }
    }

    if (not cout_redirected)
      Dune::Logging::Logging::restoreCout();

    return std::make_pair(grid, host_grid);
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GMSH_READER_HH