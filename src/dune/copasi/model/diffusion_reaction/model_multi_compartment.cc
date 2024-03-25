#ifndef DUNE_COPASI_GRID_DIMENSION
#define DUNE_COPASI_GRID_DIMENSION 2
#endif

#ifndef DUNE_COPASI_FEM_ORDER
#define DUNE_COPASI_FEM_ORDER 1
#endif

#include <dune/copasi/model/diffusion_reaction/model_multi_compartment.hh>
#include <dune/copasi/model/diffusion_reaction/model_multi_compartment_traits.hh>
#include <dune/copasi/model/diffusion_reaction/model_single_compartment_traits.hh>

#if DUNE_COPASI_GRID_DIMENSION < 2
#include <dune/grid/yaspgrid.hh>
#else
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/multidomaingrid/mdgridtraits.hh>
#include <dune/grid/multidomaingrid/multidomaingrid.hh>

namespace Dune {

#if DUNE_COPASI_GRID_DIMENSION < 2
using HostGrid =
  Dune::YaspGrid<DUNE_COPASI_GRID_DIMENSION,
                 Dune::EquidistantOffsetCoordinates<double, DUNE_COPASI_GRID_DIMENSION>>;
#else
using HostGrid = Dune::UGGrid<DUNE_COPASI_GRID_DIMENSION>;
#endif

using MDGTraits = Dune::mdgrid::FewSubDomainsTraits<DUNE_COPASI_GRID_DIMENSION, 64>;
using MDGrid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;
using SDGridView = typename MDGrid::SubDomainGrid::LeafGridView;

namespace Copasi::DiffusionReaction {

template<bool SpeciesBlocked>
using SingleCompartmentTraits = ModelSingleCompartmentPkTraits<MDGrid,
                                                               SDGridView,
                                                               DUNE_COPASI_FEM_ORDER,
                                                               double,
                                                               double,
                                                               SpeciesBlocked>;

template<bool SpeciesBlocked, bool CompartmentBlocked>
using MultiCompartmentTraits =
  ModelMultiCompartmentPkTraits<SingleCompartmentTraits<SpeciesBlocked>, CompartmentBlocked>;

template class ModelMultiCompartment<MultiCompartmentTraits<true, true>>;
template class ModelMultiCompartment<MultiCompartmentTraits<false, true>>;
template class ModelMultiCompartment<MultiCompartmentTraits<true, false>>;
template class ModelMultiCompartment<MultiCompartmentTraits<false, false>>;

} // namespace Copasi::DiffusionReaction
} // namespace Dune
