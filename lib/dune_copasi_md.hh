#ifndef DUNE_COPASI_MD_HH
#define DUNE_COPASI_MD_HH

#include <dune/copasi/model/multidomain_diffusion_reaction.hh>
#include <dune/copasi/model/multidomain_diffusion_reaction.cc>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/uggrid.hh>

namespace Dune {
namespace Copasi {

constexpr int dim = 2;
using HostGrid = Dune::UGGrid<dim>;
using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

} // namespace Dorie
} // namespace Dune

#endif // DUNE_COPASI_MD_HH