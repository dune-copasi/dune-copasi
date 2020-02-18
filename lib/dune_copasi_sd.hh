#ifndef DUNE_COPASI_SD_HH
#define DUNE_COPASI_SD_HH

#include <dune/copasi/model/diffusion_reaction.hh>
#include <dune/copasi/model/diffusion_reaction.cc>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/uggrid.hh>

namespace Dune {
namespace Copasi {

constexpr int dim = 2;
using HostGrid = Dune::UGGrid<dim>;
using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
using MDGrid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

using Grid = typename MDGrid::SubDomainGrid;
using GridView = typename Grid::Traits::LeafGridView;

} // namespace Dorie
} // namespace Dune

#endif // DUNE_COPASI_SD_HH