#include <dune/copasi/model_diffusion_reaction.hh>

namespace Dune::Copasi {
  
template<int components>
ModelDiffusionReaction<components>::ModelDiffusionReaction(
  std::shared_ptr<Grid> grid,
  const Dune::ParameterTree& config
) : ModelBase(config)
  , _grid(grid)
  , _grid_view(_grid->leafGridView())
{}

template<int components>
void ModelDiffusionReaction<components>::step()
{}

template<int components>
void ModelDiffusionReaction<components>::suggest_timestep(double dt)
{}

template<int components>
void ModelDiffusionReaction<components>::operator_setup()
{}

}