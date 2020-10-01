#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune_copasi_md.hh"

namespace Dune {
namespace Copasi {

using ModelTraits2DPk1 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid2D, 1>;
template class ModelMultiDomainDiffusionReaction<ModelTraits2DPk1>;

using ModelTraits2DPk2 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid2D, 2>;
template class ModelMultiDomainDiffusionReaction<ModelTraits2DPk2>;

#ifdef DUNE_COPASI_COMPILE_3D

using ModelTraits3DPk1 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid3D, 1>;
template class ModelMultiDomainDiffusionReaction<ModelTraits3DPk1>;

using ModelTraits3DPk2 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid3D, 2>;
template class ModelMultiDomainDiffusionReaction<ModelTraits3DPk2>;

#endif

} // namespace Dorie
} // namespace Dune
