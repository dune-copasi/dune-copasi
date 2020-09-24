#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "dune_copasi_md.hh"

namespace Dune {
namespace Copasi {

using ModelTraits2DPk0 =
  Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid2D, 0>;
template class ModelMultiDomainDiffusionReaction<ModelTraits2DPk0>;

using ModelTraits2DPk1 =
  Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid2D, 1>;
template class ModelMultiDomainDiffusionReaction<ModelTraits2DPk1>;

using ModelTraits2DPk2 =
  Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid2D, 2>;
template class ModelMultiDomainDiffusionReaction<ModelTraits2DPk2>;

using ModelTraits3DPk0 =
  Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid3D, 0>;
template class ModelMultiDomainDiffusionReaction<ModelTraits3DPk0>;

using ModelTraits3DPk1 =
  Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid3D, 1>;
template class ModelMultiDomainDiffusionReaction<ModelTraits3DPk1>;

using ModelTraits3DPk2 =
  Dune::Copasi::ModelMultiDomainP0PkDiffusionReactionTraits<Grid3D, 2>;
template class ModelMultiDomainDiffusionReaction<ModelTraits3DPk2>;

} // namespace Dorie
} // namespace Dune
