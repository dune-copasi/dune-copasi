#include <dune/copasi/model_diffusion_reaction.hh>

#include <string>

namespace Dune::Copasi {
  
template<int components>
ModelDiffusionReaction<components>::ModelDiffusionReaction(
  std::shared_ptr<Grid> grid,
  const Dune::ParameterTree& config
) : ModelBase(config)
  , _grid_view(grid->leafGridView()) // TODO: change this to a more appropriated view
{
  _state.grid = grid;
  operator_setup();
}

template<int components>
void ModelDiffusionReaction<components>::step()
{}

template<int components>
void ModelDiffusionReaction<components>::suggest_timestep(double dt)
{}

template<int components>
void ModelDiffusionReaction<components>::operator_setup()
{
  // create a finite element map
  _finite_element_map = std::make_shared<FEM>(_grid_view);

  // create (one) leaf grid function space
  LGFS leaf_grid_function_space(_grid_view,_finite_element_map);

  // create a power grid function space (lgfs is internally copied) 
  // TODO: check whether fem is hard copied or if the pointer is still vaild
  _state.grid_function_space = std::make_shared<GFS>(leaf_grid_function_space);

  // name each component
  for (int i=0; i<components; i++)
    _state.grid_function_space->child(i).name("u_"+std::to_string(i));
  
  // _constraints = std::make_unique<CC>();
  // Dune::PDELab::constraints(*_grid_function_space,*_constraints,false);

}

}