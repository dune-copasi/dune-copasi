#ifndef DUNE_COPASI_GMSH_READER_HH
#define DUNE_COPASI_GMSH_READER_HH

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/multidomaingrid.hh>

#include <dune/logging/logging.hh>

#include <algorithm>
#include <type_traits>

namespace Dune::Copasi {

template<class Grid>
class GmshReader;

template<class HostGrid, class MDGTraits>
class GmshReader<Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>>
{
public:
  using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

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