#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune_copasi_sd.hh"

namespace Dune {
namespace Copasi {

using ModelTraits2DPk0 =
  Dune::Copasi::ModelPkDiffusionReactionTraits<Grid2D, GridView2D, 0>;
template class ModelDiffusionReaction<ModelTraits2DPk0>;

#ifdef DUNE_COPASI_COMPILE_3D

using ModelTraits3DPk0 =
  Dune::Copasi::ModelPkDiffusionReactionTraits<Grid3D, GridView3D, 0>;
template class ModelDiffusionReaction<ModelTraits3DPk0>;

#endif

} // namespace Dorie
} // namespace Dune
