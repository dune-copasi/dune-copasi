#ifndef DUNE_COPASI_MODEL_FACTORY_HH
#define DUNE_COPASI_MODEL_FACTORY_HH

#include <dune/copasi/concepts/grid.hh>

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/model/diffusion_reaction_mc.hh>
#include <dune/copasi/model/diffusion_reaction_mc_traits.hh>
#include <dune/copasi/model/diffusion_reaction_sc.hh>
#include <dune/copasi/model/diffusion_reaction_sc_traits.hh>

#include <dune/copasi/parser/context.hh>

// comma separated list of fem orders to compile for each dimension
#ifndef DUNE_COPASI_1D_FEM_ORDERS
#define DUNE_COPASI_1D_FEM_ORDERS 1
#endif

#ifndef DUNE_COPASI_2D_FEM_ORDERS
#define DUNE_COPASI_2D_FEM_ORDERS 1
#endif

#ifndef DUNE_COPASI_3D_FEM_ORDERS
#define DUNE_COPASI_3D_FEM_ORDERS 1
#endif

namespace Dune::Copasi {

namespace Impl {
template<class Model, std::size_t Order, bool SpeciesBlocked>
using SingleCompartmentTraits = ModelDiffusionPkReactionTraits<typename Model::Grid,
                                                               typename Model::GridView,
                                                               Order,
                                                               typename Model::RangeQuatinty,
                                                               typename Model::TimeQuantity,
                                                               SpeciesBlocked>;

template<class Model, std::size_t Order, bool SpeciesBlocked, bool CompartmentBlocked>
using MultiCompartmentTraits = ModelMultiCompartmentDiffusionReactionPkTraits<
  SingleCompartmentTraits<Model, Order, SpeciesBlocked>,
  CompartmentBlocked>;
}

// multi-domaingrid case!
template<class Model>
  requires Concept::MultiDomainGrid<typename Model::Grid> &&
           Concept::SubDomainGrid<typename Model::GridView::Grid>
std::unique_ptr<Model>
make_model(
  const typename Model::Grid& md_grid,
  const ParameterTree& config,
  std::shared_ptr<ParserContext> parser_context  )
{

  const auto fem_orders = []() {
    if constexpr (Model::Grid::dimensionworld == 1) {
      return std::index_sequence<DUNE_COPASI_1D_FEM_ORDERS>{};
    } else if constexpr (Model::Grid::dimensionworld == 2) {
      return std::index_sequence<DUNE_COPASI_2D_FEM_ORDERS>{};
    } else if constexpr (Model::Grid::dimensionworld == 3) {
      return std::index_sequence<DUNE_COPASI_3D_FEM_ORDERS>{};
    }
    return std::index_sequence<1>{};
  }();

  auto model_config = config.sub("model");
  auto config_fem_order = model_config.get("order", std::size_t{ 1 });

  auto field_blocked = model_config.get("blocked_layout.scalar_fields", false);
  auto compartments_blocked = model_config.get("blocked_layout.compartments", false);

  std::set<std::string> compartments;
  for (const auto& component : model_config.sub("scalar_field").getSubKeys()) {
    compartments.insert(model_config[fmt::format("scalar_field.{}.compartment", component)]);
  }

  // default case when order is unknown
  auto not_know_order = [config_fem_order]() {
    throw format_exception(
      IOError{}, "Polynomail degree '{}' not available for this simulation", config_fem_order);
  };

  std::unique_ptr<Model> model;

  if (compartments.size() == 1) {
    // unroll static switch case for dynamic order case
    Dune::Hybrid::switchCases(
      fem_orders,
      config_fem_order,
      [&](auto fem_order) {
        if (field_blocked) {
          model = std::make_unique<
            ModelDiffusionReaction<Impl::SingleCompartmentTraits<Model, fem_order, true> > >(md_grid, config, std::move(parser_context));
        } else {
          model = std::make_unique<
            ModelDiffusionReaction<Impl::SingleCompartmentTraits<Model, fem_order, false> > >(md_grid, config, std::move(parser_context));
        }
      },
      not_know_order );
  } else {
    // unroll static switch case for dynamic order case
    Dune::Hybrid::switchCases(
      fem_orders,
      config_fem_order,
      [&](auto fem_order) {
        if (compartments_blocked) {
          if (field_blocked) {
            model = std::make_unique<ModelMultiCompartmentDiffusionReaction<
              Impl::MultiCompartmentTraits<Model, fem_order, true, true>>>(md_grid, config, std::move(parser_context));
          } else {
            model = std::make_unique<ModelMultiCompartmentDiffusionReaction<
              Impl::MultiCompartmentTraits<Model, fem_order, false, true>>>(md_grid, config, std::move(parser_context));
          }
        } else {
          if (field_blocked) {
            model = std::make_unique<ModelMultiCompartmentDiffusionReaction<
              Impl::MultiCompartmentTraits<Model, fem_order, true, false>>>(md_grid, config, std::move(parser_context));
          } else {
            model = std::make_unique<ModelMultiCompartmentDiffusionReaction<
              Impl::MultiCompartmentTraits<Model, fem_order, false, false>>>(md_grid, config, std::move(parser_context));
          }
        }
      } ,
      not_know_order);
  }

  assert(model); // this should not be reachable if parameters are invalid
  return model;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_FACTORY_HH
