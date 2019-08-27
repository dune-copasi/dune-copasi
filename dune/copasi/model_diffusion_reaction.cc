#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC

/**
 * This is the implementation of the ModelDiffusionReaction,
 * particularly, you want to notice that this is not an normal .cc
 * file but a header which has to be included when compiling.
 */

#include <dune/copasi/model_diffusion_reaction.hh>
#include <dune/copasi/pdelab_callable_adapter.hh>
#include <dune/copasi/pdelab_expression_adapter.hh>

#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include <string>

#include <sys/stat.h>

namespace Dune::Copasi {

template<class Grid, class GridView, int FEMorder, class OrderingTag>
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::
  ModelDiffusionReaction(std::shared_ptr<Grid> grid,
                         const Dune::ParameterTree& config,
                         GV grid_view,
                         ModelSetupPolicy setup_policy)
  : ModelBase(config)
  , _solver_logger(Logging::Logging::componentLogger(config, "solver"))
  , _components(config.sub("reaction").getValueKeys().size())
  , _config(config)
  , _grid_view(grid_view)
  , _grid(grid)
{
  setup(setup_policy);

  if (setup_policy == ModelSetupPolicy::All) {
    // _logger.debug("set initial state"_fmt);
    // auto& intial_config = _config.sub("initial");
    // ExpressionToGridFunctionAdapter<GV, RF> initial(_grid_view,
    // intial_config); set_state(initial);

    for (std::size_t i = 0; i < _states.size(); i++) {
      auto& x = _states[i].coefficients;
      auto& gfs = _states[i].grid_function_space;

      auto etity_transformation = [&](auto e) {
        if constexpr (Concept::isMultiDomainGrid<Grid>() and
                      Concept::isSubDomainGrid<typename GridView::Grid>())
          return _grid->multiDomainEntity(e);
        else
          return e;
      };
      using ET = decltype(etity_transformation);

      using Predicate = PDELab::vtk::DefaultPredicate;
      using Data =
        PDELab::vtk::DGFTreeCommonData<GFS, X, Predicate, GridView, ET>;
      std::shared_ptr<Data> data =
        std::make_shared<Data>(*gfs, *x, _grid_view, etity_transformation);
      PDELab::vtk::OutputCollector<SW, Data> collector(*_sequential_writer,
                                                       data);
      for (std::size_t k = 0; k < data->_lfs.degree(); k++)
        collector.addSolution(data->_lfs.child(k),
                              PDELab::vtk::defaultNameScheme());
    }
    _sequential_writer->write(current_time(), Dune::VTK::appendedraw);
    _sequential_writer->vtkWriter()->clear();
  }

  _logger.debug("ModelDiffusionReaction constructed"_fmt);
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::
  ~ModelDiffusionReaction()
{
  _logger.debug("ModelDiffusionReaction deconstructed"_fmt);
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
void
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::step()
{
  double dt = _config.template get<double>("time_step");

  const bool cout_redirected = Dune::Logging::Logging::isCoutRedirected();
  if (not cout_redirected)
    Dune::Logging::Logging::redirectCout(_solver_logger.name(),
                                         Dune::Logging::LogLevel::detail);


  // do time step
  for (std::size_t i = 0; i < _states.size(); i++) {
    _local_operators[i]->update(const_states());
    _temporal_local_operators[i]->update(const_states());

    auto& x = _states[i].coefficients;
    auto& gfs = _states[i].grid_function_space;
    auto x_new = std::make_shared<X>(*x);
    _one_step_methods[i]->apply(current_time(), dt, *x, *x_new);

    // accept time step
    x = x_new;

    auto etity_transformation = [&](auto e) {
      if constexpr (Concept::isMultiDomainGrid<Grid>() and
                    Concept::isSubDomainGrid<typename GridView::Grid>())
        return _grid->multiDomainEntity(e);
      else
        return e;
    };
    using ET = decltype(etity_transformation);

    using Predicate = PDELab::vtk::DefaultPredicate;
    using Data =
      PDELab::vtk::DGFTreeCommonData<GFS, X, Predicate, GridView, ET>;
    std::shared_ptr<Data> data =
      std::make_shared<Data>(*gfs, *x, _grid_view, etity_transformation);
    PDELab::vtk::OutputCollector<SW, Data> collector(*_sequential_writer, data);
    for (std::size_t k = 0; k < data->_lfs.degree(); k++)
      collector.addSolution(data->_lfs.child(k),
                            PDELab::vtk::defaultNameScheme());
  }

  if (not cout_redirected)
    Dune::Logging::Logging::restoreCout();

  current_time() += dt;

  _sequential_writer->write(current_time(), Dune::VTK::appendedraw);
  _sequential_writer->vtkWriter()->clear();
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
void
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::suggest_timestep(
  double dt)
{
  DUNE_THROW(NotImplemented, "");
}

// template<class Grid, class GridView, int FEMorder, class OrderingTag>
// template<class T>
// void
// ModelDiffusionReaction<Grid, GridView, FEMorder,
// OrderingTag>::set_state(const T& input_state)
// {
//   if constexpr (std::is_arithmetic<T>::value) {
//     _logger.trace("convert state to a vector of components"_fmt);
//     DynamicVector<RF> component_state(_components, input_state);

//     _logger.trace("set state from constant vector of components"_fmt);
//     set_state(component_state);
//   } else if constexpr (std::is_same_v<T, Dune::DynamicVector<RF>>) {
//     assert(input_state.size() == _components);

//     _logger.trace("convert vector of components to a callable"_fmt);
//     auto callable = [&](const auto& x) { return input_state; };

//     _logger.trace("set state from callable"_fmt);
//     set_state(callable);
//   } else if constexpr (Concept::isPDELabCallable<GV, T>()) {
//     _logger.trace("convert callable to a grid function"_fmt);
//     auto grid_function = makeGridFunctionFromCallable(_grid_view,
//     input_state);

//     _logger.trace("set state from grid function"_fmt);
//     set_state(grid_function);
//   } else if constexpr (Concept::isPDELabGridFunction<T>()) {
//     _logger.trace("interpolate grid function to model coefficients"_fmt);
//     Dune::PDELab::interpolate(input_state, *_gfs, *_x);
//   } else {
//     static_assert(Dune::AlwaysFalse<T>::value, "Not known input model
//     state");
//   }
// }

template<class Grid, class GridView, int FEMorder, class OrderingTag>
auto
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::
  setup_component_grid_function_space(std::string name) const
{
  _logger.trace("create a finite element map"_fmt);
  BaseFEM base_fem(_grid->leafGridView());
  auto finite_element_map = std::make_shared<FEM>(_grid_view, base_fem);

  _logger.trace("setup grid function space for component {}"_fmt, name);
  auto comp_gfs =
    std::make_shared<LGFS>(_grid->leafGridView(), finite_element_map);
  comp_gfs->name(name);

  return comp_gfs;
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
auto
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::
  setup_domain_grid_function_space(std::vector<std::string> comp_names) const
{
  _logger.debug("setup domain grid function space"_fmt);

  typename GFS::NodeStorage comp_gfs_vec(comp_names.size());

  for (std::size_t k = 0; k < comp_names.size(); k++) {
    comp_gfs_vec[k] = setup_component_grid_function_space(comp_names[k]);
  }

  _logger.trace("setup power grid function space"_fmt);
  return std::make_shared<GFS>(comp_gfs_vec);
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
void
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::
  setup_grid_function_space()
{
  // get component names
  auto& operator_splitting_config = _config.sub("operator");
  auto comp_names = operator_splitting_config.getValueKeys();
  std::sort(comp_names.begin(), comp_names.end());

  _operator_splitting.clear();
  std::size_t max_op_i = 0;
  for (auto&& var : comp_names) {
    std::size_t i = operator_splitting_config.get(var, 0);
    _operator_splitting.insert(std::make_pair(i, var));
    max_op_i = std::max(max_op_i, i + 1);
  }

  _states.resize(max_op_i);

  for (std::size_t i = 0; i < _states.size(); i++) {
    auto op_range = _operator_splitting.equal_range(i);
    std::size_t op_size = std::distance(op_range.first, op_range.second);
    assert(op_size > 0);
    std::vector<std::string> op_comp_names(op_size);
    std::transform(op_range.first,
                   op_range.second,
                   op_comp_names.begin(),
                   [](const auto& i) { return i.second; });
    _states[i].grid_function_space =
      setup_domain_grid_function_space(op_comp_names);
    _states[i].grid = _grid;
  }
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
void
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::
  setup_coefficient_vectors()
{
  _logger.debug("setup coefficient vector"_fmt);
  for (std::size_t i = 0; i < _states.size(); i++) {
    auto& x = _states[i].coefficients;
    auto& gfs = _states[i].grid_function_space;
    if (x)
      x = std::make_shared<X>(*gfs, *(x->storage()));
    else
      x = std::make_shared<X>(*gfs);
  }
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
void
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::
  setup_constraints()
{
  _logger.debug("setup constraints"_fmt);

  auto b0lambda = [&](const auto& i, const auto& x) { return false; };
  auto b0 =
    Dune::PDELab::makeBoundaryConditionFromCallable(_grid_view, b0lambda);

  _logger.trace("assemble constraints"_fmt);
  _constraints = std::make_unique<CC>();
  for (std::size_t i = 0; i < _states.size(); i++) {
    auto& gfs = _states[i].grid_function_space;
    Dune::PDELab::constraints(b0, *gfs, *_constraints);

    _logger.info("constrained dofs: {} of {}"_fmt,
                 _constraints->size(),
                 gfs->globalSize());
  }
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
auto
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::
  setup_local_operator(std::size_t i) const
{
  _logger.trace("setup local operators {}"_fmt, i);

  _logger.trace("create spatial local operator {}"_fmt, i);
  FE finite_element;

  std::vector<std::size_t> gfs_components, map_operator;
  auto& operator_splitting_config = _config.sub("operator");
  auto comp_names = operator_splitting_config.getValueKeys();
  std::sort(comp_names.begin(), comp_names.end());

  for (std::size_t j = 0; j < comp_names.size(); j++) {
    std::size_t k =
      operator_splitting_config.template get<std::size_t>(comp_names[j]);
    map_operator.push_back(k);
    if (k == i)
      gfs_components.push_back(j);
  }

  CM coefficient_mapper(map_operator, i);

  auto local_operator = std::make_shared<LOP>(
    _grid_view, _config, finite_element, gfs_components, coefficient_mapper);

  _logger.trace("create temporal local operator"_fmt);
  auto temporal_local_operator = std::make_shared<TLOP>(
    _grid_view, _config, finite_element, gfs_components, coefficient_mapper);

  return std::make_pair(local_operator, temporal_local_operator);
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
void
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::
  setup_local_operators()
{
  _logger.trace("setup local operators"_fmt);

  _local_operators.resize(_states.size());
  _temporal_local_operators.resize(_states.size());

  for (std::size_t i = 0; i < _states.size(); i++) {
    auto operators = setup_local_operator(i);
    _local_operators[i] = operators.first;
    _temporal_local_operators[i] = operators.second;
  }
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
void
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::
  setup_grid_operators()
{
  _logger.debug("setup grid operators"_fmt);
  _spatial_grid_operators.resize(_states.size());
  _temporal_grid_operators.resize(_states.size());
  _grid_operators.resize(_states.size());

  for (std::size_t i = 0; i < _states.size(); i++) {
    auto& gfs = _states[i].grid_function_space;
    auto& lop = _local_operators[i];
    auto& tlop = _temporal_local_operators[i];

    MBE mbe((int)pow(3, dim));

    _logger.trace("create spatial grid operator {}"_fmt, i);
    _spatial_grid_operators[i] = std::make_shared<GOS>(
      *gfs, *_constraints, *gfs, *_constraints, *lop, mbe);

    _logger.trace("create temporal grid operator {}"_fmt, i);
    _temporal_grid_operators[i] = std::make_shared<GOT>(
      *gfs, *_constraints, *gfs, *_constraints, *tlop, mbe);

    _logger.trace("create instationary grid operator {}"_fmt, i);
    _grid_operators[i] = std::make_shared<GOI>(*_spatial_grid_operators[i],
                                               *_temporal_grid_operators[i]);
  }
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
void
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::setup_solvers()
{
  _logger.debug("setup solvers"_fmt);
  _linear_solvers.resize(_states.size());
  _nonlinear_solvers.resize(_states.size());
  _time_stepping_methods.resize(_states.size());
  _one_step_methods.resize(_states.size());

  for (std::size_t i = 0; i < _states.size(); i++) {
    auto& x = _states[i].coefficients;
    _logger.trace("create linear solver"_fmt);
    _linear_solvers[i] = std::make_shared<LS>(5000, false);

    _logger.trace("create nonlinear solver"_fmt);
    _nonlinear_solvers[i] =
      std::make_shared<NLS>(*_grid_operators[i], *x, *_linear_solvers[i]);
    _nonlinear_solvers[i]->setReassembleThreshold(0.0);
    _nonlinear_solvers[i]->setVerbosityLevel(2);
    _nonlinear_solvers[i]->setReduction(1e-8);
    _nonlinear_solvers[i]->setMinLinearReduction(1e-10);
    _nonlinear_solvers[i]->setMaxIterations(25);
    _nonlinear_solvers[i]->setLineSearchMaxIterations(10);

    _logger.trace("select and prepare time-stepping scheme"_fmt);
    using AlexMethod = Dune::PDELab::Alexander2Parameter<double>;
    _time_stepping_methods[i] = std::make_shared<AlexMethod>();
    _one_step_methods[i] = std::make_shared<OSM>(
      *_time_stepping_methods[i], *_grid_operators[i], *_nonlinear_solvers[i]);
    _one_step_methods[i]->setVerbosityLevel(2);
  }
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
void
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::
  setup_vtk_writer()
{
  _logger.debug("setup vtk writer"_fmt);

  _writer = std::make_shared<W>(_grid_view, Dune::VTK::conforming);
  auto config_writer = _config.sub("writer");
  std::string file_name = config_writer.template get<std::string>("file_name");
  std::string path = config_writer.get("path", file_name);
  struct stat st;

  if (stat(path.c_str(), &st) != 0) {
    int stat = 0;
    stat = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
    if (stat != 0 && stat != -1)
      std::cout << "Error: Cannot create directory " << path << std::endl;
  }

  _sequential_writer = std::make_shared<SW>(_writer, file_name, path, path);
}

template<class Grid, class GridView, int FEMorder, class OrderingTag>
void
ModelDiffusionReaction<Grid, GridView, FEMorder, OrderingTag>::setup(
  ModelSetupPolicy setup_policy)
{
  _logger.trace("setup operator started"_fmt);

  if (setup_policy >= ModelSetupPolicy::GridFunctionSpace)
    setup_grid_function_space();
  if (setup_policy >= ModelSetupPolicy::CoefficientVector)
    setup_coefficient_vectors();
  if (setup_policy >= ModelSetupPolicy::Constraints)
    setup_constraints();
  if (setup_policy >= ModelSetupPolicy::LocalOperator)
    setup_local_operators();
  if (setup_policy >= ModelSetupPolicy::GridOperator)
    setup_grid_operators();
  if (setup_policy >= ModelSetupPolicy::Solver)
    setup_solvers();
  if (setup_policy >= ModelSetupPolicy::Writer)
    setup_vtk_writer();
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC