#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune_copasi_md.hh"

namespace Dune {
namespace Copasi {

using ModelTraits0 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid, 0>;
template class ModelMultiDomainDiffusionReaction<ModelTraits0>;

} // namespace Dorie
} // namespace Dune