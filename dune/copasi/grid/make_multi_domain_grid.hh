#ifndef DUNE_COPASI_GRID_MAKE_MULTI_DOMAIN_GRID_HH
#define DUNE_COPASI_GRID_MAKE_MULTI_DOMAIN_GRID_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/ostream_to_spdlog.hh>
#include <dune/copasi/common/axis_names.hh>
#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/parser/context.hh>
#include <dune/copasi/grid/cell_data.hh>
#include <dune/copasi/grid/cell_data_parser.hh>

#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid/uggridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/common/parametertree.hh>

#include <memory>

#include <fmt/ranges.h>

namespace Dune::Copasi {

/**
 * @brief Creates a multi domain grid from a config file
 * From given configuration file with grid and compartment section as described
 * in the program documentation (--grid.* and --compartment.*), this funcition
 * will read or create a grid and its sub-domains (or compartments).
 * It will also assign a unique subdomain 'id' to each compartment in the
 * configuration file. Finally, it will create a cell data with the contents of the
 * `gmsh_id` if the file was read by gmsh, and the data files parsed with `cell_data_parser`. 
 *
 * @tparam MDGrid                   The type of the multidomain-grid
* @tparam CellDataType              The type to store values in the cell data object
 * @param config                    The configuration file
 * @param parser_context            A parser context to interpret parsed math expression
 * @return std::pair<std::unique_ptr<MDGrid>, std::unique_ptr<CellData<typename MDGrid::LevelGridView, CellDataType>>> Pointer to the resulting grid and cell data
 */
template<Concept::MultiDomainGrid MDGrid, class CellDataType = double>
std::tuple<std::unique_ptr<MDGrid>, std::unique_ptr<CellData<typename MDGrid::LevelGridView, CellDataType>>>
make_multi_domain_grid_with_cell_data(Dune::ParameterTree& config, std::shared_ptr<const ParserContext> parser_context = {})
{
  using HostGrid = typename MDGrid::HostGrid;
  using MDEntity = typename MDGrid::template Codim<0>::Entity;

  const auto& grid_config = config.sub("grid");
  std::size_t const dim = grid_config.get("dimension", std::size_t{ MDGrid::dimensionworld });
  if (dim != MDGrid::dimensionworld) {
    throw format_exception(IOError{},
                           "Executable grid dimension does not match with input grid dimension!");
  }

  auto out_guard = ostream2spdlog();
  auto& compartments_config = config.sub("compartments");

  std::unique_ptr<HostGrid> host_grid_ptr;
  std::unique_ptr<MDGrid> md_grid_ptr;
  std::vector<int> gmsh_physical_entity;

  std::vector<std::pair<std::string,std::function<bool(const MDEntity&)>>> compartments;

  if (grid_config.hasKey("path")) {
    auto grid_path = grid_config.template get<std::string>("path");
    spdlog::info("Reading grid file '{}'", grid_path);
    auto host_grid_factory = Dune::GridFactory<HostGrid>{};
    auto host_grid_parser = Dune::GmshReaderParser<HostGrid>{ host_grid_factory, false, false };

    host_grid_parser.read(grid_path);
    gmsh_physical_entity = host_grid_parser.elementIndexMap();

    host_grid_ptr = host_grid_factory.createGrid();
  } else {
    std::array<unsigned int, MDGrid::dimensionworld> cells{};
    std::fill_n(begin(cells), MDGrid::dimensionworld, 1);
    cells = grid_config.get("cells", cells);
    auto upper_right =
      grid_config.get("extensions", FieldVector<double, MDGrid::dimensionworld>(1.));
    auto lower_left = grid_config.get("origin", FieldVector<double, MDGrid::dimensionworld>(0.));
    upper_right += lower_left;
    if (MDGrid::dimensionworld == 1) {
      host_grid_ptr =
        StructuredGridFactory<HostGrid>::createCubeGrid(lower_left, upper_right, cells);
    } else {
      host_grid_ptr =
        StructuredGridFactory<HostGrid>::createSimplexGrid(lower_left, upper_right, cells);
    }
  }

  typename MDGrid::MDGridTraits const md_grid_traits{};
  md_grid_ptr = std::make_unique<MDGrid>(std::move(host_grid_ptr), md_grid_traits);

  auto level = grid_config.get("refinement_level", int{ 0 });
  if (level > 0) {
    spdlog::info("Applying refinement of level: {}", level);
    md_grid_ptr->globalRefine(level);
  }

  std::cout << fmt::format("MultiDomainGrid: ");
  gridinfo(*md_grid_ptr, "    ");

  auto cell_data = std::make_unique<CellData<typename MDGrid::LevelGridView, CellDataType>>(md_grid_ptr->levelGridView(0));
  
  // inject gmsh id into cell_data
  if (not gmsh_physical_entity.empty())
    cell_data->addData("gmsh_id", gmsh_physical_entity);

  // parse and inject data-file contents into cell_data
  if (grid_config.hasSub("cell_data"))
    cell_data_parser(grid_config, *cell_data);

  std::vector<CellDataType> cell_values(cell_data->size());
  std::vector<bool> cell_mask(cell_data->size());

  // assign ids to dynamically generated compartments
  for (const auto& compartment : compartments_config.getSubKeys()) {
    auto& compartment_config = compartments_config.sub(compartment);
    if (not compartment_config.hasKey("id")) {
      // assign a unique id to the config file
      std::size_t id = compartments.size();
      compartment_config["id"] = std::to_string(id);
      const auto& type = compartment_config.get("type", "expression");
      if (type == "expression") {
        auto parser_type =
          string2parser.at(compartment_config.get("parser_type", default_parser_str));
        std::shared_ptr const parser_ptr = make_parser(parser_type);
        for (std::size_t i = 0; i != cell_data->size(); ++i)
          parser_ptr->define_variable(cell_data->keys()[i], &cell_values[i]);
        parser_ptr->set_expression(compartment_config["expression"]);
        if (parser_context)
          parser_context->add_context(*parser_ptr);
        auto position = std::make_shared<FieldVector<double, MDGrid::dimensionworld>>();
        for (std::size_t i = 0; i != axis_names.size(); ++i) {
          auto pos_arg = fmt::format("position_{}", axis_names[i]);
          if (i < MDGrid::dimensionworld) {
            parser_ptr->define_variable(pos_arg, &(*position)[i]);
          } else {
            parser_ptr->define_constant(pos_arg, 0.);
          }
        }
        parser_ptr->compile();
        compartments.emplace_back(compartment, [position, parser_ptr, &cell_data, &cell_values, &cell_mask](const MDEntity& entity) {
          (*position) = entity.geometry().center();
          cell_data->getData(entity, cell_values, cell_mask);
          return FloatCmp::ne(std::invoke(*parser_ptr), 0.);
        });
      } else {
        throw format_exception(NotImplemented{}, "Not know type '{}' of compartment", type);
      }
    }
  }

  // check that all compartments have an id
  if (compartments_config.getSubKeys().size() != compartments.size())  {
    throw format_exception(InvalidStateException{}, "Compartment ids were setup with wrong number of sub-domains");
  }
  for (const auto& compartment : compartments_config.getSubKeys()) {
    if (not compartments_config.sub(compartment).hasKey("id")) {
      spdlog::warn("Compartment '{}' were not assigned an 'id'!", compartment);
    }
  }

  spdlog::info("Applying sub-domain partition");
  md_grid_ptr->startSubDomainMarking();
  std::size_t max_domains_per_entity = 0;
  std::set<std::size_t> used_sub_domains;
  for (const auto& cell : elements(md_grid_ptr->leafGridView())) {
    for (std::size_t id = max_domains_per_entity = 0; id != compartments.size(); ++id) {
      if (compartments[id].second(cell)) {
        used_sub_domains.insert(id);
        ++max_domains_per_entity;
        md_grid_ptr->addToSubDomain(id, cell);
      }
    }
    if (max_domains_per_entity > static_cast<std::size_t>(md_grid_traits.maxSubDomainIndex()) + 1u) {
      throw format_exception(GridError{},
                             "This version of dune-copasi does not support to"
                             " have more than {} domains per entity!",
                             md_grid_traits.maxSubDomainIndex() + 1);
    }
  }

  if (used_sub_domains.size() != compartments.size())  {
    spdlog::warn("Grid partition uses sub-domain{} {} but config file has {} compartment{}",
      used_sub_domains.size() == 1 ? "" : "s", used_sub_domains, 
      compartments.size(), compartments.size() == 1 ? "" : "s");
  }

  md_grid_ptr->preUpdateSubDomains();
  md_grid_ptr->updateSubDomains();
  md_grid_ptr->postUpdateSubDomains();

  for (std::size_t id = 0; id != compartments.size(); ++id) {
    std::cout << fmt::format("  SubDomain {{{}: {}}}", id, compartments[id].first);
    gridinfo(md_grid_ptr->subDomain(id), "      ");
  }

  return std::make_tuple(std::move(md_grid_ptr), std::move(cell_data));
}


/**
 * @copydoc make_multi_domain_grid_with_cell_data
 * 
 * @return std::unique_ptr<MDGrid> Pointer to the resulting grid
 */
template<Concept::MultiDomainGrid MDGrid>
std::unique_ptr<MDGrid>
make_multi_domain_grid(Dune::ParameterTree& config,
                       std::shared_ptr<const ParserContext> parser_context = {})
{
  return std::get<0>(make_multi_domain_grid_with_cell_data<MDGrid>(config, parser_context));
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_MAKE_MULTI_DOMAIN_GRID_HH
