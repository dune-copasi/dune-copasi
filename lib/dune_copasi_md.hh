#ifndef DUNE_COPASI_MD_HH
#define DUNE_COPASI_MD_HH

#include <dune/copasi/model/multidomain_diffusion_reaction.hh>
#include <dune/copasi/model/multidomain_diffusion_reaction.cc>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/uggrid.hh>

namespace Dune {
namespace Copasi {

using HostGrid2D = Dune::UGGrid<2>;
using MDGTraits2D = Dune::mdgrid::DynamicSubDomainCountTraits<2, 1>;
using Grid2D = Dune::mdgrid::MultiDomainGrid<HostGrid2D, MDGTraits2D>;

#ifdef DUNE_COPASI_COMPILE_3D

using HostGrid3D = Dune::UGGrid<3>;
using MDGTraits3D = Dune::mdgrid::DynamicSubDomainCountTraits<3, 1>;
using Grid3D = Dune::mdgrid::MultiDomainGrid<HostGrid3D, MDGTraits3D>;

#endif

} // namespace Dorie
} // namespace Dune

#endif // DUNE_COPASI_MD_HH
