#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune_copasi_sd.hh"

namespace Dune {
namespace Copasi {

using ModelTraits0 =
  Dune::Copasi::ModelPkDiffusionReactionTraits<Grid, GridView, 0>;
template class ModelDiffusionReaction<ModelTraits0>;

} // namespace Dorie
} // namespace Dune