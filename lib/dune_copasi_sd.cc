#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/model/diffusion_reaction.cc>
#include <dune/copasi/model/diffusion_reaction.hh>

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

using ModelTraits1 =
  Dune::Copasi::ModelPkDiffusionReactionTraits<Grid, GridView, 1>;
template class ModelDiffusionReaction<ModelTraits1>;

using ModelTraits2 =
  Dune::Copasi::ModelPkDiffusionReactionTraits<Grid, GridView, 2>;
template class ModelDiffusionReaction<ModelTraits2>;

} // namespace Dorie
} // namespace Dune
