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

using ModelTraits1 =
  Dune::Copasi::ModelMultiDomainDiffusionReactionTraits<Grid, 1>;
template class ModelMultiDomainDiffusionReaction<ModelTraits1>;

using ModelTraits2 =
  Dune::Copasi::ModelMultiDomainDiffusionReactionTraits<Grid, 2>;
template class ModelMultiDomainDiffusionReaction<ModelTraits2>;

using Ordering = Dune::PDELab::EntityBlockedOrderingTag;
constexpr Dune::Copasi::JacobianMethod Jac =
  Dune::Copasi::JacobianMethod::Numerical;
using ModelTraits1Num =
  Dune::Copasi::ModelMultiDomainDiffusionReactionTraits<Grid, 1, Ordering, Jac>;
template class ModelMultiDomainDiffusionReaction<ModelTraits1Num>;

} // namespace Dorie
} // namespace Dune
