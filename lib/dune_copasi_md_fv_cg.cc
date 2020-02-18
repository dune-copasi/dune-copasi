#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune_copasi_md.hh"

namespace Dune {
namespace Copasi {

using ModelTraits01 =
  Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid, 1>;
template class ModelMultiDomainDiffusionReaction<ModelTraits01>;

using ModelTraits02 =
  Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid, 2>;
template class ModelMultiDomainDiffusionReaction<ModelTraits02>;

} // namespace Dorie
} // namespace Dune