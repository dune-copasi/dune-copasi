#ifndef DUNE_COPASI_MODEL_MAKE_INITIAL_HH
#define DUNE_COPASI_MODEL_MAKE_INITIAL_HH

#include <dune/copasi/concepts/grid.hh>

#include <dune/copasi/model/local_equations/functor_factory.hh>

#include <dune/pdelab/common/trace.hh>

#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dune/grid/concepts/grid.hh>

#include <dune/common/exceptions.hh>

#include <spdlog/spdlog.h>

#include <optional>
#include <string>
#include <unordered_map>

namespace Dune::Copasi {

template<class GridFunction, Dune::Concept::Grid Grid>
[[nodiscard]] inline static std::unordered_map<std::string, GridFunction>
make_initial(const Grid& grid,
             const ParameterTree& config,
             const FunctorFactory<Grid::dimensionworld>& functor_factory)
{
  TRACE_EVENT("dune", "InitialCondition");

  std::unordered_map<std::string, GridFunction> functions;

  auto time = config.get("time_step_operator.time_begin", double{ 0. });

  const auto& compartments_config = config.sub("compartments", true);

  for (const auto& compartment : compartments_config.getSubKeys()) {

    // find grid view to use
    auto sub_grid_view = [&]() {
      const auto& compartment_config = compartments_config.sub(compartment);
      if constexpr (Concept::MultiDomainGrid<Grid>) {
        using SubDomainIndex = typename Grid::SubDomainIndex;
        auto const domain_id = compartment_config.template get<SubDomainIndex>("id");
        if (domain_id > grid.maxAssignedSubDomainIndex()) {
          DUNE_THROW(IOError, "\tCompartment ID does not exist in multi-domain grid");
        }
        return grid.subDomain(domain_id).leafGridView();
      } else {
        auto const domain_id = compartment_config.template get<std::size_t>("id");
        if (domain_id != 0) {
          DUNE_THROW(IOError, "\tCompartment ID does not exist in grid");
        }
        return grid.leafGridView();
      }
    }();

    auto components = config.sub("scalar_field", true).getSubKeys();
    for (const auto& component : components) {
      auto prefix = fmt::format("scalar_field.{}.initial", component);
      if (not config.hasSub(prefix)) {
        continue;
      }
      auto local_domain = std::make_shared<LocalDomain<Grid::dimensionworld>>();
      local_domain->time = time;
      auto fcn = functor_factory.make_scalar(prefix, config.sub(prefix), *local_domain, false);
      if (not fcn) {
        continue;
      }
      auto shared_fnc = std::make_shared<decltype(fcn)>(std::move(fcn));
      using GridDomain =
        typename decltype(sub_grid_view)::template Codim<0>::Geometry::GlobalCoordinate;
      std::function<double(GridDomain)> analytic_fcn =
        [lcl_domain = std::move(local_domain),
         _fnc = std::move(shared_fnc)](GridDomain global_position) {
          lcl_domain->position = global_position;
          return std::invoke(*_fnc);
        };

      auto grid_function =
        Dune::Functions::makeAnalyticGridViewFunction(std::move(analytic_fcn), sub_grid_view);

      std::string const component_id = fmt::format("{}.{}", compartment, component);
      auto [_, inserted] = functions.try_emplace(component_id, std::move(grid_function));
      if (not inserted) {
        DUNE_THROW(IOError, "\tTwo or more components share the same name");
      }
    }
  }

  return functions;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MAKE_INITIAL_HH
