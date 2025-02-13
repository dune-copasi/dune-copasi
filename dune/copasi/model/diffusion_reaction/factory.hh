#ifndef DUNE_COPASI_MODEL_FACTORY_HH
#define DUNE_COPASI_MODEL_FACTORY_HH

#include <dune/copasi/concepts/grid.hh>

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/grid/cell_data.hh>
#include <dune/copasi/model/diffusion_reaction/model_multi_compartment.hh>
#include <dune/copasi/model/diffusion_reaction/model_multi_compartment_traits.hh>
#include <dune/copasi/model/diffusion_reaction/model_single_compartment.hh>
#include <dune/copasi/model/diffusion_reaction/model_single_compartment_traits.hh>

#include <dune/copasi/model/functor_factory_parser.hh>

// comma separated list of fem orders to compile for each dimension
#ifndef DUNE_COPASI_1D_DIFFUSION_REACTION_FEM_ORDERS
#define DUNE_COPASI_1D_DIFFUSION_REACTION_FEM_ORDERS 1
#endif

#ifndef DUNE_COPASI_2D_DIFFUSION_REACTION_FEM_ORDERS
#define DUNE_COPASI_2D_DIFFUSION_REACTION_FEM_ORDERS 1
#endif

#ifndef DUNE_COPASI_3D_DIFFUSION_REACTION_FEM_ORDERS
#define DUNE_COPASI_3D_DIFFUSION_REACTION_FEM_ORDERS 1
#endif

namespace Dune::Copasi::DiffusionReaction {

namespace Impl {
template<class Model, std::size_t Order, bool SpeciesBlocked>
using SingleCompartmentTraits = ModelSingleCompartmentPkTraits<typename Model::Grid,
                                                               typename Model::GridView,
                                                               Order,
                                                               typename Model::RangeQuatinty,
                                                               typename Model::TimeQuantity,
                                                               SpeciesBlocked>;

template<class Model, std::size_t Order, bool SpeciesBlocked, bool CompartmentBlocked>
using MultiCompartmentTraits = ModelMultiCompartmentPkTraits<
  SingleCompartmentTraits<Model, Order, SpeciesBlocked>,
  CompartmentBlocked>;
}

// multi-domaingrid case!
template<class Model>
  requires Concept::MultiDomainGrid<typename Model::Grid> &&
           Concept::SubDomainGrid<typename Model::GridView::Grid>
std::unique_ptr<Model>
make_model(
  const ParameterTree& config,
  std::shared_ptr<const FunctorFactory<Model::Grid::dimensionworld>> functor_factory = nullptr,
  std::shared_ptr<const CellData<typename Model::Grid::LeafGridView, typename Model::RangeQuatinty>> cell_data = nullptr
)
{
  if (not functor_factory) {
    functor_factory = std::make_shared<FunctorFactoryParser<Model::Grid::dimensionworld>>();
  }

  const auto fem_orders = []() {
    if constexpr (Model::Grid::dimensionworld == 1) {
      return std::index_sequence<DUNE_COPASI_1D_DIFFUSION_REACTION_FEM_ORDERS>{};
    } else if constexpr (Model::Grid::dimensionworld == 2) {
      return std::index_sequence<DUNE_COPASI_2D_DIFFUSION_REACTION_FEM_ORDERS>{};
    } else if constexpr (Model::Grid::dimensionworld == 3) {
      return std::index_sequence<DUNE_COPASI_3D_DIFFUSION_REACTION_FEM_ORDERS>{};
    } else {
      return std::index_sequence<1>{};
    }
  }();

  auto config_fem_order = config.get("order", std::size_t{ 1 });

  auto field_blocked = config.get("blocked_layout.scalar_fields", false);
  auto compartments_blocked = config.get("blocked_layout.compartments", false);

  std::set<std::string> compartments;
  for (const auto& component : config.sub("scalar_field").getSubKeys()) {
    compartments.insert(config[fmt::format("scalar_field.{}.compartment", component)]);
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
            ModelSingleCompartment<Impl::SingleCompartmentTraits<Model, fem_order, true>>>(
            functor_factory, cell_data);
        } else {
          model = std::make_unique<
            ModelSingleCompartment<Impl::SingleCompartmentTraits<Model, fem_order, false>>>(
            functor_factory, cell_data);
        }
      },
      not_know_order);
  } else {
    // unroll static switch case for dynamic order case
    Dune::Hybrid::switchCases(
      fem_orders,
      config_fem_order,
      [&](auto fem_order) {
        if (compartments_blocked) {
          if (field_blocked) {
            model = std::make_unique<ModelMultiCompartment<
              Impl::MultiCompartmentTraits<Model, fem_order, true, true>>>(functor_factory, cell_data);
          } else {
            model = std::make_unique<ModelMultiCompartment<
              Impl::MultiCompartmentTraits<Model, fem_order, false, true>>>(functor_factory, cell_data);
          }
        } else {
          if (field_blocked) {
            model = std::make_unique<ModelMultiCompartment<
              Impl::MultiCompartmentTraits<Model, fem_order, true, false>>>(functor_factory, cell_data);
          } else {
            model = std::make_unique<ModelMultiCompartment<
              Impl::MultiCompartmentTraits<Model, fem_order, false, false>>>(functor_factory, cell_data);
          }
        }
      },
      not_know_order);
  }

  assert(model); // this should not be reachable if parameters are invalid
  return model;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_FACTORY_HH
