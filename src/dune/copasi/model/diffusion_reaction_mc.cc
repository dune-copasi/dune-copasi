#include <dune-copasi-config.h>

#ifndef DUNE_COPASI_GRID_DIMENSION
#define DUNE_COPASI_GRID_DIMENSION 2
#endif

#ifndef DUNE_COPASI_FEM_ORDER
#define DUNE_COPASI_FEM_ORDER 1
#endif

#include <dune/copasi/model/diffusion_reaction_mc.hh>
#include <dune/copasi/model/diffusion_reaction_mc_traits.hh>
#include <dune/copasi/model/diffusion_reaction_sc.hh>
#include <dune/copasi/model/diffusion_reaction_sc_traits.hh>
#include <dune/copasi/model/local_equations/functor_factory_parser.hh>

#if DUNE_COPASI_GRID_DIMENSION < 2
#include <dune/grid/yaspgrid.hh>
#else
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/multidomaingrid.hh>

namespace Dune {

#if DUNE_COPASI_GRID_DIMENSION < 2
using HostGrid = Dune::YaspGrid<DUNE_COPASI_GRID_DIMENSION>;
#else
using HostGrid = Dune::UGGrid<DUNE_COPASI_GRID_DIMENSION>;
#endif

using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<DUNE_COPASI_GRID_DIMENSION, 10>;
using MDGrid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;
using SDGridView = typename MDGrid::SubDomainGrid::LeafGridView;

namespace Copasi {

template<bool SpeciesBlocked>
using SingleCompartmentTraits = ModelDiffusionPkReactionTraits<MDGrid,
                                                               SDGridView,
                                                               DUNE_COPASI_FEM_ORDER,
                                                               double,
                                                               double,
                                                               false,
                                                               SpeciesBlocked>;

extern template class ModelDiffusionReaction<SingleCompartmentTraits<true>>;
extern template class ModelDiffusionReaction<SingleCompartmentTraits<false>>;

template<bool SpeciesBlocked, bool CompartmentBlocked>
using MultiCompartmentTraits =
  ModelMultiCompartmentDiffusionReactionPkTraits<SingleCompartmentTraits<SpeciesBlocked>,
                                                 CompartmentBlocked>;

template class ModelMultiCompartmentDiffusionReaction<MultiCompartmentTraits<true, true>>;
template class ModelMultiCompartmentDiffusionReaction<MultiCompartmentTraits<false, true>>;
template class ModelMultiCompartmentDiffusionReaction<MultiCompartmentTraits<true, false>>;
template class ModelMultiCompartmentDiffusionReaction<MultiCompartmentTraits<false, false>>;

} // namespace Copasi
} // namespace Dune
