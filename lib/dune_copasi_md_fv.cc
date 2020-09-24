#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune_copasi_md.hh"

namespace Dune {
namespace Copasi {

using ModelTraits2DPk0 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid2D, 0>;
template class ModelMultiDomainDiffusionReaction<ModelTraits2DPk0>;

using ModelTraits3DPk0 =
  Dune::Copasi::ModelMultiDomainPkDiffusionReactionTraits<Grid3D, 0>;
template class ModelMultiDomainDiffusionReaction<ModelTraits3DPk0>;

} // namespace Dorie
} // namespace Dune
