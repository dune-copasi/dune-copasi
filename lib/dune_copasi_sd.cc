#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

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

// only FV

using ModelTraits0 =
  Dune::Copasi::ModelPkDiffusionReactionTraits<Grid, GridView, 0>;
template class ModelDiffusionReaction<ModelTraits0>;

// only CG

using ModelTraits1 =
  Dune::Copasi::ModelPkDiffusionReactionTraits<Grid, GridView, 1>;
template class ModelDiffusionReaction<ModelTraits1>;

using ModelTraits2 =
  Dune::Copasi::ModelPkDiffusionReactionTraits<Grid, GridView, 2>;
template class ModelDiffusionReaction<ModelTraits2>;

// mixed FV and CG

using ModelTraits01 =
  Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid, GridView, 1>;
template class ModelDiffusionReaction<ModelTraits01>;

using ModelTraits02 =
  Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid, GridView, 2>;
template class ModelDiffusionReaction<ModelTraits02>;

} // namespace Dorie
} // namespace Dune
