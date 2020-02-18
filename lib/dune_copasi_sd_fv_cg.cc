#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune_copasi_sd.hh"

namespace Dune {
namespace Copasi {

using ModelTraits01 =
  Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid, GridView, 1>;
template class ModelDiffusionReaction<ModelTraits01>;

using ModelTraits02 =
  Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid, GridView, 2>;
template class ModelDiffusionReaction<ModelTraits02>;

} // namespace Dorie
} // namespace Dune