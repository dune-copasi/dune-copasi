#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune_copasi_md.hh"

namespace Dune {
namespace Copasi {

using ModelTraits1 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid, 1>;
template class ModelMultiDomainDiffusionReaction<ModelTraits1>;

using ModelTraits2 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid, 2>;
template class ModelMultiDomainDiffusionReaction<ModelTraits2>;

} // namespace Dorie
} // namespace Dune