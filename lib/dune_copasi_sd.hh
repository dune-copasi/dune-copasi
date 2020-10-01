#ifndef DUNE_COPASI_SD_HH
#define DUNE_COPASI_SD_HH

#include <dune/copasi/model/diffusion_reaction.hh>
#include <dune/copasi/model/diffusion_reaction.cc>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/uggrid.hh>

namespace Dune {
namespace Copasi {

using HostGrid2D = Dune::UGGrid<2>;
using MDGTraits2D = Dune::mdgrid::DynamicSubDomainCountTraits<2, 1>;
using MDGrid2D = Dune::mdgrid::MultiDomainGrid<HostGrid2D, MDGTraits2D>;

using Grid2D = typename MDGrid2D::SubDomainGrid;
using GridView2D = typename Grid2D::Traits::LeafGridView;

#ifdef DUNE_COPASI_COMPILE_3D

using HostGrid3D = Dune::UGGrid<3>;
using MDGTraits3D = Dune::mdgrid::DynamicSubDomainCountTraits<3, 1>;
using MDGrid3D = Dune::mdgrid::MultiDomainGrid<HostGrid3D, MDGTraits3D>;

using Grid3D = typename MDGrid3D::SubDomainGrid;
using GridView3D = typename Grid3D::Traits::LeafGridView;

#endif

} // namespace Dorie
} // namespace Dune

#endif // DUNE_COPASI_SD_HH
