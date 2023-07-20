#ifndef DUNE_COPASI_GRID_MAKE_MULTI_DOMAIN_GRID_HH
#define DUNE_COPASI_GRID_MAKE_MULTI_DOMAIN_GRID_HH

#include <dune-copasi-config.h>

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/ostream_to_spdlog.hh>
#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/parser/context.hh>

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

template<Concept::MultiDomainGrid MDGrid>
std::unique_ptr<MDGrid>
make_multi_domain_grid(Dune::ParameterTree& config,
                       std::shared_ptr<const ParserContext> parser_context = {})
{
  using HostGrid = typename MDGrid::HostGrid;
  using Entity = typename MDGrid::template Codim<0>::Entity;

  const auto& grid_config = config.sub("grid");
  std::size_t const dim = grid_config.get("dimension", std::size_t{ MDGrid::dimensionworld });
  if (dim != MDGrid::dimensionworld) {
    throw format_exception(IOError{},
                           "Executable grid dimension does not match with input grid dimension!");
  }

  auto out_guard = ostream2spdlog();

  std::unique_ptr<HostGrid> host_grid_ptr;
  std::unique_ptr<MDGrid> md_grid_ptr;

  std::unique_ptr<Dune::MultipleCodimMultipleGeomTypeMapper<typename MDGrid::LeafGridView>> mcmg;
  std::map<std::size_t, std::string> compartment_names;
  std::vector<std::function<bool(const Entity&)>> compartment_fncs;

  if (grid_config.hasKey("path")) {
    auto grid_path = grid_config.template get<std::string>("path");
    spdlog::info("Reading grid file '{}'", grid_path);
    auto host_grid_factory = Dune::GridFactory<HostGrid>{};
    auto host_grid_parser = Dune::GmshReaderParser<HostGrid>{ host_grid_factory, false, false };

    host_grid_parser.read(grid_path);
    const auto& gmsh_physical_entity = host_grid_parser.elementIndexMap();

    host_grid_ptr = host_grid_factory.createGrid();

    // get compartments with fixed id
    for (const auto& compartment : config.sub("compartments").getSubKeys()) {
      const auto& compartment_config = config.sub("compartments").sub(compartment, true);
      if (compartment_config.hasKey("id")) {
        auto const id = compartment_config.template get<std::size_t>("id");
        if (compartment_names.contains(id)) {
          compartment_names[id] += fmt::format(", {}", compartment);
        } else {
          compartment_names[id] = compartment;
        }
        compartment_fncs.resize(id + 1);
        compartment_fncs[id] = [&gmsh_physical_entity, &mcmg, id](const Entity& entity) {
          Entity _entity = entity;
          while (_entity.hasFather()) {
            _entity = _entity.father();
          }
          assert(_entity.level() == 0);
          auto gmhs_id = gmsh_physical_entity[mcmg->index(_entity)];
          assert(gmhs_id > 0);
          return id + 1 == static_cast<std::size_t>(gmhs_id);
        };
      }
    }

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

  // fix id of other compartments
  for (const auto& compartment : config.sub("compartments").getSubKeys()) {
    auto& compartment_config = config.sub("compartments").sub(compartment);
    if (not compartment_config.hasKey("id")) {
      std::size_t id = 0;
      std::size_t count = 0;
      while (compartment_names.contains(id = count++)) {
      }
      compartment_names[id] = compartment;
      compartment_config["id"] = std::to_string(id);
      compartment_fncs.resize(id + 1);
      const auto& type = compartment_config["type"];
      if (type == "expression") {
        auto parser_type =
          string2parser.at(compartment_config.get("parser_type", default_parser_str));
        std::shared_ptr const parser_ptr = make_parser(parser_type);
        if (parser_context)
          parser_context->add_context(*parser_ptr);
        auto position = std::make_shared<FieldVector<double, MDGrid::dimensionworld>>();
        std::vector<std::string> dim_name = { "x", "y", "z" };
        for (std::size_t i = 0; i != MDGrid::dimensionworld; ++i) {
          parser_ptr->define_variable(fmt::format("position_{}", dim_name.at(i)), &(*position)[i]);
        }
        parser_ptr->set_expression(compartment_config["expression"]);
        parser_ptr->compile();
        compartment_fncs[id] = [position, parser_ptr](const Entity& entity) {
          (*position) = entity.geometry().center();
          return FloatCmp::ne(std::invoke(*parser_ptr), 0.);
        };
      } else {
        throw format_exception(NotImplemented{}, "Not know type '{}' of compartment", type);
      }
    }
  }

  typename MDGrid::MDGridTraits const md_grid_traits(compartment_fncs.size());
  md_grid_ptr = std::make_unique<MDGrid>(std::move(host_grid_ptr), md_grid_traits);
  mcmg = std::make_unique<MultipleCodimMultipleGeomTypeMapper<typename MDGrid::LeafGridView>>(
    md_grid_ptr->leafGridView(), mcmgElementLayout());

  auto level = grid_config.get("refinement_level", int{ 0 });
  if (level > 0) {
    spdlog::info("Applying refinement of level: {}", level);
    md_grid_ptr->globalRefine(level);
  }

  spdlog::info("Applying sub-domain partition");
  md_grid_ptr->startSubDomainMarking();
  std::size_t max_domains_per_entity = 0;
  for (const auto& cell : elements(md_grid_ptr->leafGridView())) {
    for (std::size_t id = max_domains_per_entity = 0; id != compartment_fncs.size(); ++id) {
      if (compartment_fncs[id](cell)) {
        ++max_domains_per_entity;
        md_grid_ptr->addToSubDomain(id, cell);
      }
    }
    if (max_domains_per_entity > static_cast<std::size_t>(md_grid_traits.maxSubDomainIndex() + 1)) {
      throw format_exception(GridError{},
                             "This version of dune-copasi does not support to"
                             " have more than {} domains per entity!",
                             md_grid_traits.maxSubDomainIndex() + 1);
    }
  }

  md_grid_ptr->preUpdateSubDomains();
  md_grid_ptr->updateSubDomains();
  md_grid_ptr->postUpdateSubDomains();

  std::cout << fmt::format("MultiDomainGrid: ");
  gridinfo(*md_grid_ptr, "    ");
  for (auto [id, compartment_name] : compartment_names) {
    std::cout << fmt::format("  SubDomain {{{}: {}}}", id, compartment_name);
    gridinfo(md_grid_ptr->subDomain(id), "      ");
  }

  return md_grid_ptr;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_MAKE_MULTI_DOMAIN_GRID_HH
