#ifndef DUNE_COPASI_MODEL_FACTORY_HH
#define DUNE_COPASI_MODEL_FACTORY_HH

#include <dune/copasi/concepts/grid.hh>

#include <dune/copasi/model/diffusion_reaction_mc.hh>
#include <dune/copasi/model/diffusion_reaction_mc_traits.hh>
#include <dune/copasi/model/diffusion_reaction_sc.hh>
#include <dune/copasi/model/diffusion_reaction_sc_traits.hh>

#include <dune/copasi/model/local_equations/functor_factory_parser.hh>

#ifndef DUNE_COPASI_MAX_FEM_ORDER
#define DUNE_COPASI_MAX_FEM_ORDER 1
#endif

namespace Dune::Copasi {

namespace Impl {
template<class Model, std::size_t Order, bool SpeciesBlocked>
using SingleCompartmentTraits = ModelDiffusionPkReactionTraits<typename Model::Grid,
                                                               typename Model::GridView,
                                                               Order,
                                                               typename Model::RangeQuatinty,
                                                               typename Model::TimeQuantity,
                                                               false,
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
  const ParameterTree& config,
  std::shared_ptr<const FunctorFactory<Model::Grid::dimensionworld>> functor_factory = nullptr)
{

  if (not functor_factory) {
    auto parser_type = string2parser.at(config.get("parser_type", default_parser_str));
    functor_factory =
      std::make_shared<FunctorFactoryParser<Model::Grid::dimensionworld>>(parser_type);
  }

  auto dynamic_fem_order = config.get("order", std::size_t{ 1 });
  auto order_range =
    Dune::range(Indices::_1, Dune::index_constant<DUNE_COPASI_MAX_FEM_ORDER + 1>{});

  auto field_blocked = config.get("blocked_layout.scalar_fields", false);
  auto compartments_blocked = config.get("blocked_layout.compartments", false);

  std::set<std::string> compartments;
  if (not config.hasSub("scalar_field")) {
    DUNE_THROW(IOError, "A model must have an 'scalar_field' section");
  }
  for (const auto& component : config.sub("scalar_field", true).getSubKeys()) {
    compartments.insert(config[fmt::format("scalar_field.{}.compartment", component)]);
  }

  std::unique_ptr<Model> model;
  if (compartments.size() == 1) {
    Dune::Hybrid::switchCases(order_range, dynamic_fem_order, [&](auto static_fem_order) {
      if (field_blocked) {
        model = std::make_unique<
          ModelDiffusionReaction<Impl::SingleCompartmentTraits<Model, static_fem_order, true>>>(
          functor_factory);
      } else {
        model = std::make_unique<
          ModelDiffusionReaction<Impl::SingleCompartmentTraits<Model, static_fem_order, false>>>(
          functor_factory);
      }
    });
  } else {
    Dune::Hybrid::switchCases(order_range, dynamic_fem_order, [&](auto static_fem_order) {
      if (compartments_blocked) {
        if (field_blocked) {
          model = std::make_unique<ModelMultiCompartmentDiffusionReaction<
            Impl::MultiCompartmentTraits<Model, static_fem_order, true, true>>>(functor_factory);
        } else {
          model = std::make_unique<ModelMultiCompartmentDiffusionReaction<
            Impl::MultiCompartmentTraits<Model, static_fem_order, false, true>>>(functor_factory);
        }
      } else {
        if (field_blocked) {
          model = std::make_unique<ModelMultiCompartmentDiffusionReaction<
            Impl::MultiCompartmentTraits<Model, static_fem_order, true, false>>>(functor_factory);
        } else {
          model = std::make_unique<ModelMultiCompartmentDiffusionReaction<
            Impl::MultiCompartmentTraits<Model, static_fem_order, false, false>>>(functor_factory);
        }
      }
    });
  }

  if (not model) {
    DUNE_THROW(NotImplemented, "Not known model for the given configuration");
  }
  return model;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_FACTORY_HH
