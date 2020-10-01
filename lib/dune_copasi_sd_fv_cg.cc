#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune_copasi_sd.hh"

namespace Dune {
namespace Copasi {


using ModelTraits2DPk0 =
  Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid2D, GridView2D, 0>;
template class ModelDiffusionReaction<ModelTraits2DPk0>;

using ModelTraits2DPk1 =
  Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid2D, GridView2D, 1>;
template class ModelDiffusionReaction<ModelTraits2DPk1>;

using ModelTraits2DPk2 =
  Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid2D, GridView2D, 2>;
template class ModelDiffusionReaction<ModelTraits2DPk2>;

#ifdef DUNE_COPASI_COMPILE_3D

using ModelTraits3DPk0 =
  Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid3D, GridView3D, 0>;
template class ModelDiffusionReaction<ModelTraits3DPk0>;

using ModelTraits3DPk1 =
  Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid3D, GridView3D, 1>;
template class ModelDiffusionReaction<ModelTraits3DPk1>;

using ModelTraits3DPk2 =
  Dune::Copasi::ModelP0PkDiffusionReactionTraits<Grid3D, GridView3D, 2>;
template class ModelDiffusionReaction<ModelTraits3DPk2>;

#endif

} // namespace Dorie
} // namespace Dune
