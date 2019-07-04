#include <dune/copasi/model_diffusion_reaction.hh>
#include <dune/copasi/pdelab_callable_adapter.hh>
#include <dune/copasi/pdelab_expression_adapter.hh>

#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include <dune/common/dynvector.hh>

#include <string>

#include <sys/stat.h>

namespace Dune::Copasi {

template<class Grid, class GridView>
ModelDiffusionReaction<Grid,GridView>::ModelDiffusionReaction(
  std::shared_ptr<Grid> grid,
  GV grid_view,
  const Dune::ParameterTree& config)
  : ModelBase(config)
  , _components(config.sub("reaction").getValueKeys().size())
  , _config(config)
  , _grid_view(grid_view)
  , _x(_state.coefficients)
  , _gfs(_state.grid_function_space)
{
  _state.grid = grid;
  operator_setup();

  auto& intial_config = _config.sub("initial");
  ExpressionToGridFunctionAdapter<GV, RF> initial(_grid_view, intial_config);

  set_state(initial);

  _writer = std::make_shared<W>(_grid_view, Dune::VTK::conforming);
  std::string filename = config.get("name", "output");
  struct stat st;

  if (stat(filename.c_str(), &st) != 0) {
    int stat = 0;
    stat = mkdir(filename.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    if (stat != 0 && stat != -1)
      std::cout << "Error: Cannot create directory " << filename << std::endl;
  }

  auto names = intial_config.getValueKeys();
  _sequential_writer = std::make_shared<SW>(_writer, filename, filename, "");
  // add data field for all components of the space to the VTK writer
  for (int i = 0; i < names.size(); ++i) {
    auto dgf = std::make_shared<const DGF>(_gfs, _x, i * _dof_per_component, (i + 1) * _dof_per_component);
    _sequential_writer->addVertexData(dgf, names[i]);
  }
  _sequential_writer->write(current_time(), Dune::VTK::appendedraw);
  _sequential_writer->vtkWriter()->clear();

  _logger.debug("ModelDiffusionReaction constructed"_fmt);
}

template<class Grid, class GridView>
ModelDiffusionReaction<Grid,GridView>::~ModelDiffusionReaction()
{
  _logger.debug("ModelDiffusionReaction destroyed"_fmt);
}

template<class Grid, class GridView>
void
ModelDiffusionReaction<Grid,GridView>::step()
{
  double dt = 0.001;

  // do time step
  auto x_new = std::make_shared<X>(*_x);
  _one_step_method->apply(current_time(), dt, *_x, *x_new);

  // accept time step
  _x = x_new;
  current_time() += dt;

  auto names = _config.sub("initial").getValueKeys();
  for (int i = 0; i < names.size(); ++i) {
    auto dgf = std::make_shared<const DGF>(_gfs, _x, i * _dof_per_component, (i + 1) * _dof_per_component);
    _sequential_writer->addVertexData(dgf, names[i]);
  }
  _sequential_writer->write(current_time(), Dune::VTK::appendedraw);
  _sequential_writer->vtkWriter()->clear();
}

template<class Grid, class GridView>
void
ModelDiffusionReaction<Grid,GridView>::suggest_timestep(double dt)
{
  DUNE_THROW(NotImplemented, "");
}

template<class Grid, class GridView>
template<class T>
void
ModelDiffusionReaction<Grid,GridView>::set_state(const T& input_state)
{
  if constexpr (std::is_arithmetic<T>::value) {
    _logger.trace("convert state to a vector of components"_fmt);
    DynamicVector<RF> component_state(_components, input_state);

    _logger.trace("set state from constant vector of components"_fmt);
    set_state(component_state);
  } else if constexpr (std::is_same_v<T, Dune::DynamicVector<RF>>) {
    assert(input_state.size() == _components);

    _logger.trace("convert vector of components to a callable"_fmt);
    auto callable = [&](const auto& x) { return input_state; };

    _logger.trace("set state from callable"_fmt);
    set_state(callable);
  } else if constexpr (Concept::isPDELabCallable<GV, T>()) {
    _logger.trace("convert callable to a grid function"_fmt);
    auto grid_function = makeGridFunctionFromCallable(_grid_view, input_state);

    _logger.trace("set state from grid function"_fmt);
    set_state(grid_function);
  } else if constexpr (Concept::isPDELabGridFunction<T>()) {
    _logger.trace("interpolate grid function to model coefficients"_fmt);
    Dune::PDELab::interpolate(input_state, *_gfs, *_x);
  } else {
    static_assert(Dune::AlwaysFalse<T>::value, "Not known input model state");
  }
}

template<class Grid, class GridView>
void
ModelDiffusionReaction<Grid,GridView>::operator_setup()
{
  BaseFEM base_fem(_grid_view);
  _dof_per_component = base_fem.maxLocalSize(); // todo: fix this

  _logger.trace("create a finite element map"_fmt);
  _finite_element_map = std::make_shared<FEM>(base_fem, _components);

  _logger.trace("create (one) leaf grid function space"_fmt);
  LGFS leaf_grid_function_space(_grid_view, _finite_element_map);

  _logger.trace("create a power grid function space"_fmt);
  _gfs = std::make_shared<GFS>(leaf_grid_function_space);

  _logger.trace("name grid function space"_fmt);
  _gfs->name("u");

  _logger.trace("create vector backend"_fmt);
  if (not _x)
    _x = std::make_shared<X>(*_gfs);
  else
    _x = std::make_shared<X>(*_gfs, *(_x->storage()));

  auto b0lambda = [&](const auto& i, const auto& x) { return false; };
  auto b0 =
    Dune::PDELab::makeBoundaryConditionFromCallable(_grid_view, b0lambda);

  _logger.trace("assemble constraints"_fmt);
  _constraints = std::make_unique<CC>();
  Dune::PDELab::constraints(b0, *_gfs, *_constraints);

  _logger.info(
    "constrained dofs: {} of {}"_fmt, _constraints->size(), _gfs->globalSize());

  _logger.trace("create spatial local operator"_fmt);
  auto entity_it = _grid_view.template begin<0>();
  FE finite_element;
  _local_operator = std::make_shared<LOP>(_grid_view, _config, finite_element);

  _logger.trace("create temporal local operator"_fmt);
  _temporal_local_operator =
    std::make_shared<TLOP>(_grid_view, _config, finite_element);

  MBE mbe((int)pow(3, dim));

  _logger.trace("create spatial grid operator"_fmt);
  _spatial_grid_operator = std::make_shared<GOS>(
    *_gfs, *_constraints, *_gfs, *_constraints, *_local_operator, mbe);

  _temporal_grid_operator = std::make_shared<GOT>(
    *_gfs, *_constraints, *_gfs, *_constraints, *_temporal_local_operator, mbe);

  _logger.trace("create instationary grid operator"_fmt);
  _grid_operator =
    std::make_shared<GOI>(*_spatial_grid_operator, *_temporal_grid_operator);

  _logger.trace("create linear solver"_fmt);
  _linear_solver = std::make_shared<LS>(5000, false);

  _logger.trace("create nonlinear solver"_fmt);
  _nonlinear_solver =
    std::make_shared<NLS>(*_grid_operator, *_x, *_linear_solver);
  _nonlinear_solver->setReassembleThreshold(0.0);
  _nonlinear_solver->setVerbosityLevel(2);
  _nonlinear_solver->setReduction(1e-8);
  _nonlinear_solver->setMinLinearReduction(1e-10);
  _nonlinear_solver->setMaxIterations(25);
  _nonlinear_solver->setLineSearchMaxIterations(10);

  _logger.trace("select and prepare time-stepping scheme"_fmt);
  using AlexMethod = Dune::PDELab::Alexander2Parameter<double>;
  _time_stepping_method = std::make_shared<AlexMethod>();
  _one_step_method = std::make_shared<OSM>(
    *_time_stepping_method, *_grid_operator, *_nonlinear_solver);
  _one_step_method->setVerbosityLevel(2);
}

} // namespace Dune::Copasi