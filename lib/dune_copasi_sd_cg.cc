#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune_copasi_sd.hh"

namespace Dune {
namespace Copasi {

using ModelTraits1 =
  Dune::Copasi::ModelPkDiffusionReactionTraits<Grid, GridView, 1>;
template class ModelDiffusionReaction<ModelTraits1>;

using ModelTraits2 =
  Dune::Copasi::ModelPkDiffusionReactionTraits<Grid, GridView, 2>;
template class ModelDiffusionReaction<ModelTraits2>;

} // namespace Dorie
} // namespace Dune