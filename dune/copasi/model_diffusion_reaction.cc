#include <dune/copasi/model_diffusion_reaction.hh>

#include <string>

namespace Dune::Copasi {
  
template<int components, class Param>
ModelDiffusionReaction<components,Param>::ModelDiffusionReaction(
  std::shared_ptr<Grid> grid,
  const Dune::ParameterTree& config
) : ModelBase(config)
  , _grid_view(grid->leafGridView()) // TODO: change this to a more appropriated view
{
  _parameterization = std::make_shared<Param>();
  _state.grid = grid;
  operator_setup();
  // set_state(1.);

  _logger.debug("ModelDiffusionReaction constructed"_fmt);
}

template<int components, class Param>
ModelDiffusionReaction<components,Param>::~ModelDiffusionReaction()
{
  _logger.debug("ModelDiffusionReaction destroyed"_fmt); 
}

template<int components, class Param>
void ModelDiffusionReaction<components,Param>::step()
{}

template<int components, class Param>
void ModelDiffusionReaction<components,Param>::suggest_timestep(double dt)
{}

template<int components, class Param>
template<class T>
void ModelDiffusionReaction<components,Param>::set_state(const T& input_state)
{
  if constexpr (std::is_arithmetic<T>::value)
  {
    // convert state to a vector of components
    Dune::FieldVector<RF,components> component_state( RF{input_state} );

    // set state from constant vector of components
    set_state(component_state);
  }
  else if constexpr (std::is_same_v<T,Dune::FieldVector<RF,components>>)
  {
    // convert state to a callable
    auto callable = [&](const auto& e, const auto& x)
    {
      return input_state;
    };

    // set state from callable
    // set_state(callable);

    // convert callable to a grid function
    auto grid_function = PDELab::makeGridFunctionFromCallable(_grid_view, 
                                                              callable);

    // interpolate grid function to model coefficients
    Dune::PDELab::interpolate(grid_function,
                              *(_state.grid_function_space),
                              *(_state.coefficients));
  }
  // else if constexpr (is_pdelab_callable<GV,T>::value)
  // {
  //   // convert callable to a grid function
  //   auto grid_function = PDELab::makeGridFunctionFromCallable(_grid_view, 
  //                                                             input_state);

  //   // set state from grid function
  //   set_state(grid_function);
  // }
  // else if constexpr (is_grid_function<T>::value)
  // {
  //   // interpolate grid function to model coefficients
  //   Dune::PDELab::interpolate(grid_function,
  //                             *(_state.grid_function_space),
  //                             *(_state.coefficients));
  // }
  // else if constexpr (std::is_same<T,ModelState>::value)
  // {
  //   DUNE_THROW(NotImplemented,"...");
  // }
  // else if constexpr (std::is_same<T,ConstModelState>::value)
  // {
  //   DUNE_THROW(NotImplemented,"...");
  // }
  else
  {
    static_assert(Dune::AlwaysFalse<T>::value,"Not known input model state");
  } 
}

template<int components, class Param>
void ModelDiffusionReaction<components,Param>::operator_setup()
{
  // reference to grid function space pointer
  auto& gfs = _state.grid_function_space;

  _logger.trace("create a finite element map"_fmt);
  _finite_element_map = std::make_shared<FEM>(_grid_view);

  _logger.trace("create (one) leaf grid function space"_fmt);
  LGFS leaf_grid_function_space(_grid_view,_finite_element_map);

  _logger.trace("create a power grid function space"_fmt);
  gfs = std::make_shared<GFS>(leaf_grid_function_space);

  _logger.trace("name each component"_fmt);
  for (int i=0; i<components; i++)
  {
    gfs->child(i).name("u_"+std::to_string(i));
    _logger.trace("component name {}: {}"_fmt, i, gfs->child(i).name());
  }
 

  auto b0lambda = [&](const auto& i, const auto& x)
    {return _parameterization->b(i,x);};
  auto b0 = Dune::PDELab::makeBoundaryConditionFromCallable(_grid_view,b0lambda);
  using B = Dune::PDELab::PowerConstraintsParameters<decltype(b0),components>;
  B b(b0);

  _logger.trace("assemble constraints"_fmt);
  _constraints = std::make_unique<CC>();
  Dune::PDELab::constraints(b,*gfs,*_constraints);

  _logger.info("constrained dofs: {} of {}"_fmt, _constraints->size(), 
                                                 gfs->globalSize());

  _logger.trace("create spatial local operator"_fmt);
  auto entity_it = _grid_view.template begin<0>();
  auto finite_element = _finite_element_map->find(*entity_it);
  LOP lop(*_parameterization,finite_element);

  _logger.trace("create temporal local operator"_fmt);
  TLOP tlop(finite_element);

  MBE mbe((int)pow(3,dim));

  _logger.trace("create spatial grid operator"_fmt);
  GOS gos(*gfs,*_constraints,*gfs,*_constraints,lop,mbe);

  _logger.trace("create temporal grid operator"_fmt);
  GOT got(*gfs,*_constraints,*gfs,*_constraints,tlop,mbe);

  _logger.trace("create instationary grid operator"_fmt);
  GOI goi(gos,got);

}

}