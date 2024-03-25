

#ifndef DUNE_COPASI_GRID_DIMENSION
#define DUNE_COPASI_GRID_DIMENSION 2
#endif

#ifndef DUNE_COPASI_FEM_ORDER
#define DUNE_COPASI_FEM_ORDER 1
#endif

#include <dune/copasi/model/diffusion_reaction/model_single_compartment.hh>
#include <dune/copasi/model/diffusion_reaction/model_single_compartment_traits.hh>
#include <dune/copasi/model/functor_factory_parser.hh>

#if DUNE_COPASI_GRID_DIMENSION < 2
#include <dune/grid/yaspgrid.hh>
#else
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/multidomaingrid/mdgridtraits.hh>
#include <dune/grid/multidomaingrid/multidomaingrid.hh>

namespace Dune {

#if DUNE_COPASI_GRID_DIMENSION < 2
using HostGrid = Dune::YaspGrid<DUNE_COPASI_GRID_DIMENSION, Dune::EquidistantOffsetCoordinates<double,DUNE_COPASI_GRID_DIMENSION>>;
#else
using HostGrid = Dune::UGGrid<DUNE_COPASI_GRID_DIMENSION>;
#endif
using HostGridView = typename HostGrid::Traits::LeafGridView;

using MDGTraits = Dune::mdgrid::FewSubDomainsTraits<DUNE_COPASI_GRID_DIMENSION, 64>;
using MDGrid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;
using SDGridView = typename MDGrid::SubDomainGrid::Traits::LeafGridView;

namespace Copasi {

template class FunctorFactoryParser<DUNE_COPASI_GRID_DIMENSION>;

namespace DiffusionReaction {

template class ModelSingleCompartment<ModelSingleCompartmentPkTraits<HostGrid,
                                                                     HostGridView,
                                                                     DUNE_COPASI_FEM_ORDER,
                                                                     double,
                                                                     double,
                                                                     false>>;
template class ModelSingleCompartment<ModelSingleCompartmentPkTraits<HostGrid,
                                                                     HostGridView,
                                                                     DUNE_COPASI_FEM_ORDER,
                                                                     double,
                                                                     double,
                                                                     true>>;

template class ModelSingleCompartment<ModelSingleCompartmentPkTraits<MDGrid,
                                                                     SDGridView,
                                                                     DUNE_COPASI_FEM_ORDER,
                                                                     double,
                                                                     double,
                                                                     false>>;
template class ModelSingleCompartment<ModelSingleCompartmentPkTraits<MDGrid,
                                                                     SDGridView,
                                                                     DUNE_COPASI_FEM_ORDER,
                                                                     double,
                                                                     double,
                                                                     true>>;

} // namespace DiffusionReaction
} // namespace Copasi
} // namespace Dune
