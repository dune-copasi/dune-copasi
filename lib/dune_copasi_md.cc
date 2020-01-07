#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/model/multidomain_diffusion_reaction.cc>
#include <dune/copasi/model/multidomain_diffusion_reaction.hh>

#include <dune/grid/multidomaingrid.hh>

#include <dune/grid/uggrid.hh>

namespace Dune {
namespace Copasi {

constexpr int dim = 2;
using HostGrid = Dune::UGGrid<dim>;
using MDGTraits = Dune::mdgrid::DynamicSubDomainCountTraits<dim, 1>;
using Grid = Dune::mdgrid::MultiDomainGrid<HostGrid, MDGTraits>;

// Only FV

using ModelTraits0 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid, 0>;
template class ModelMultiDomainDiffusionReaction<ModelTraits0>;

// Only CG

using ModelTraits1 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid, 1>;
template class ModelMultiDomainDiffusionReaction<ModelTraits1>;

using ModelTraits2 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid, 2>;
template class ModelMultiDomainDiffusionReaction<ModelTraits2>;

// mixed FV and CG

using ModelTraits01 =
  Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid, 1>;
template class ModelMultiDomainDiffusionReaction<ModelTraits01>;

using ModelTraits02 =
  Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid, 2>;
template class ModelMultiDomainDiffusionReaction<ModelTraits02>;


} // namespace Dorie
} // namespace Dune