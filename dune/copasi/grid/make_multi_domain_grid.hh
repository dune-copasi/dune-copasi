#ifndef DUNE_COPASI_GRID_MAKE_MULTI_DOMAIN_GRID_HH
#define DUNE_COPASI_GRID_MAKE_MULTI_DOMAIN_GRID_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/ostream_to_spdlog.hh>
#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/parser/context.hh>
#include <dune/copasi/parser/grid_context.hh>

#include <dune/grid/common/exceptions.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/uggrid/uggridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/common/parametertree.hh>

#include <memory>

namespace Dune::Copasi {

/**
 * @brief Creates a multi domain grid from a config file
 * From given configuration file with grid and compartment section as described
 * in the program documentation (--grid.* and --compartment.*), this funcition
 * will read or create a grid and its sub-domains (or compartments).
 * Additionally, it will assign a unique subdomain 'id' to each compartment in the
 * configuration file.
 *
 * @tparam MDGrid                   The type of the multidomain-grid
 * @param config                    The configuration file
 * @param parser_context            A parser context to interpret parsed math expression
 * @return std::unique_ptr<MDGrid>  Pointer to the resulting grid
 */
template<Concept::MultiDomainGrid MDGrid>
std::shared_ptr<MDGrid>
make_multi_domain_grid(Dune::ParameterTree& config,
                       std::shared_ptr<const ParserContext> parser_context = {},
                       std::shared_ptr<ParserGridContext<MDGrid>> parser_grid_context = {})
{
  using HostGrid = typename MDGrid::HostGrid;
  using HostEntity = typename HostGrid::template Codim<0>::Entity;

  auto& grid_config = config.sub("grid");
  std::size_t const dim = grid_config.get("dimension", std::size_t{ MDGrid::dimensionworld });
  if (dim != MDGrid::dimensionworld) {
    throw format_exception(IOError{},
                           "Executable grid dimension does not match with input grid dimension!");
  }

  auto out_guard = ostream2spdlog();
  auto& compartments_config = config.sub("compartments");

  std::shared_ptr<HostGrid> host_grid_ptr;
  std::shared_ptr<MDGrid> md_grid_ptr;

  std::vector<std::pair<std::string,std::function<bool(const HostEntity&)>>> compartments;

  double gmsh_id = std::numeric_limits<double>::max();
  std::function<void(HostEntity)> update_gmsh_id;
  std::function<double(HostEntity)> get_gmsh_id;

  if (grid_config.hasKey("path")) {
    auto grid_path = grid_config.template get<std::string>("path");
    spdlog::info("Reading grid file '{}'", grid_path);
    auto host_grid_factory = Dune::GridFactory<HostGrid>{};
    auto host_grid_parser = Dune::GmshReaderParser<HostGrid>{ host_grid_factory, false, false };

    host_grid_parser.read(grid_path);
    auto& gmsh_physical_entity = host_grid_parser.elementIndexMap();

    host_grid_ptr = host_grid_factory.createGrid();
    auto mcmg = MultipleCodimMultipleGeomTypeMapper<typename HostGrid::LeafGridView>(host_grid_ptr->leafGridView(), mcmgElementLayout());

    update_gmsh_id = [gmsh_physical_entity, mcmg, &gmsh_id](const HostEntity& entity) {
      HostEntity _entity = entity;
      while (_entity.hasFather()) {
        _entity = _entity.father();
      }
      assert(_entity.level() == 0);
      gmsh_id = gmsh_physical_entity[mcmg.index(_entity)];
    };

    // Create a get grid gmsh_id function to add to the parserContext (TO DO: @Dylan)
    get_gmsh_id = [gmsh_physical_entity, mcmg](const HostEntity& entity) {
      HostEntity _entity = entity;
      while (_entity.hasFather()) {
        _entity = _entity.father();
      }
      assert(_entity.level() == 0);
      return gmsh_physical_entity[mcmg.index(_entity)];
    };

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
        if (update_gmsh_id)
          parser_ptr->define_variable("gmsh_id", &gmsh_id);
        parser_ptr->set_expression(compartment_config["expression"]);
        if (parser_context)
          parser_context->add_context(*parser_ptr);
        auto position = std::make_shared<FieldVector<double, MDGrid::dimensionworld>>();
        const std::vector<std::string> dim_name = { "x", "y", "z" };
        for (std::size_t i = 0; i != 3; ++i) {
          auto pos_arg = fmt::format("position_{}", dim_name.at(i));
          if (i < MDGrid::dimensionworld) {
            parser_ptr->define_variable(pos_arg, &(*position)[i]);
          } else {
            parser_ptr->define_constant(pos_arg, 0.);
          }
        }
        parser_ptr->compile();
        compartments.emplace_back(compartment, [position, parser_ptr, update_gmsh_id](const HostEntity& entity) {
          (*position) = entity.geometry().center();
          if (update_gmsh_id)
            update_gmsh_id(entity);
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

  typename MDGrid::MDGridTraits const md_grid_traits(compartments.size());
  md_grid_ptr = std::make_unique<MDGrid>(std::move(host_grid_ptr), md_grid_traits);

  auto level = grid_config.get("refinement_level", int{ 0 });
  if (level > 0) {
    spdlog::info("Applying refinement of level: {}", level);
    md_grid_ptr->globalRefine(level);
  }

  spdlog::info("Applying sub-domain partition");
  md_grid_ptr->startSubDomainMarking();
  std::size_t max_domains = 0;
  std::size_t max_domains_per_entity = 0;
  for (const auto& cell : elements(md_grid_ptr->leafGridView())) {
    for (std::size_t id = max_domains_per_entity = 0; id != compartments.size(); ++id) {
      if (compartments[id].second(md_grid_ptr->hostEntity(cell))) {
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
    max_domains = std::max(max_domains, max_domains_per_entity);
  }

  if (max_domains != compartments.size())  {
    spdlog::warn("Grid has {} sub-domain{} but config has {} compartment{}",
      max_domains, max_domains == 1 ? "" : "s",
      compartments.size(), compartments.size() == 1 ? "" : "s");
  }

  md_grid_ptr->preUpdateSubDomains();
  md_grid_ptr->updateSubDomains();
  md_grid_ptr->postUpdateSubDomains();

  std::cout << fmt::format("MultiDomainGrid: ");
  gridinfo(*md_grid_ptr, "    ");
  for (std::size_t id = 0; id != compartments.size(); ++id) {
    std::cout << fmt::format("  SubDomain {{{}: {}}}", id, compartments[id].first);
    gridinfo(md_grid_ptr->subDomain(id), "      ");
  }

  parser_grid_context->configure(grid_config, md_grid_ptr, host_grid_ptr);

  return md_grid_ptr;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_MAKE_MULTI_DOMAIN_GRID_HH
