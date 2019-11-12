#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC

/**
 * This is the implementation of the ModelDiffusionReaction,
 * particularly, you want to notice that this is not an normal .cc
 * file but a header which has to be included when compiling.
 */

#include <dune/copasi/common/muparser_data_handler.hh>
#include <dune/copasi/common/pdelab_expression_adapter.hh>
#include <dune/copasi/concepts/pdelab.hh>
#include <dune/copasi/model/diffusion_reaction.hh>

#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include <string>

#include <sys/stat.h>

namespace Dune::Copasi {

template<class Traits>
ModelDiffusionReaction<Traits>::ModelDiffusionReaction(
  std::shared_ptr<Grid> grid,
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

  if constexpr (not Traits::is_sub_model) // todo fix this conditional for
                                          // something else
    if (setup_policy == ModelSetupPolicy::All) {
      MuParserDataHandler<TIFFGrayscale<unsigned short>> mu_data_handler;
      if (_config.hasSub("data"))
        mu_data_handler.add_tiff_functions(_config.sub("data"));

      auto initial_muparser = get_muparser_initial(_config, _grid_view, false);

      for (auto&& mu_grid_function : initial_muparser) {
        mu_data_handler.set_functions(mu_grid_function.parser());
        mu_grid_function.set_time(current_time());
        mu_grid_function.compile_parser();
      }
      set_initial(initial_muparser);
      write_states();
    }

  _logger.debug("ModelDiffusionReaction constructed"_fmt);
}

template<class Traits>
ModelDiffusionReaction<Traits>::~ModelDiffusionReaction()
{
  _logger.debug("ModelDiffusionReaction deconstructed"_fmt);
}

template<class Traits>
template<class GFGridView>
auto
ModelDiffusionReaction<Traits>::get_muparser_initial(
  const ParameterTree& model_config,
  const GFGridView& gf_grid_view,
  bool compile)
{
  if (not model_config.hasSub("initial"))
    DUNE_THROW(IOError, "configuration file does not have initial section");

  const auto& initial_config = model_config.sub("initial");
  const auto& vars = initial_config.getValueKeys();

  using GridFunction = ExpressionToGridFunctionAdapter<GFGridView, RF>;
  std::vector<GridFunction> functions;

  for (std::size_t i = 0; i < vars.size(); i++) {
    functions.emplace_back(gf_grid_view, initial_config[vars[i]], compile);
    assert(functions.size() == i + 1);
    if (compile)
      functions[i].compile_parser();
  }

  return functions;
}

template<class Traits>
template<class GF>
void
ModelDiffusionReaction<Traits>::set_initial(std::vector<GF>& initial)
{
  static_assert(Concept::isPDELabGridFunction<GF>(),
                "GF is not a PDElab grid functions");
  static_assert(std::is_same_v<typename GF::Traits::GridViewType, GV>,
                "GF has to have the same grid view as the templated grid view");
  static_assert(std::is_same_v<typename GF::Traits::GridViewType,
                               typename Grid::Traits::LeafGridView>,
                "GF has to have the same grid view as the templated grid");
  static_assert((int)GF::Traits::dimDomain == (int)Grid::dimension,
                "GF has to have domain dimension equal to the grid");
  static_assert(GF::Traits::dimRange == 1,
                "GF has to have range dimension equal to 1");

  _logger.debug("set initial state from grid functions"_fmt);

  if (initial.size() != _config.sub("diffusion").getValueKeys().size())
    DUNE_THROW(RangeError, "Wrong number of grid functions");

  for (auto& [op, state] : _states) {
    _logger.trace("interpolation of operator {}"_fmt, op);
    std::size_t comp_size = state.grid_function_space->degree();
    std::vector<std::shared_ptr<GF>> functions(comp_size);

    auto& operator_config = _config.sub("operator");
    auto comp_names = operator_config.getValueKeys();
    std::sort(comp_names.begin(), comp_names.end());

    std::size_t count_i = 0; // index for variables in operator
    std::size_t count_j = 0; // index for all variables
    for (const auto& var : comp_names) {
      assert(count_i < comp_names.size());
      if (op != operator_config.template get<std::size_t>(var))
        continue;

      _logger.trace("creating grid function for variable: {}"_fmt, var);
      functions[count_i] = stackobject_to_shared_ptr<GF>(initial[count_j]);
      functions[count_i]->set_time(current_time());
      count_i++;
      count_j++;
    }

    PDELab::DynamicPowerGridFunction<GF> comp_initial(functions);
    Dune::PDELab::interpolate(
      comp_initial, *state.grid_function_space, *state.coefficients);
  }
}

template<class Traits>
auto
ModelDiffusionReaction<Traits>::setup_component_grid_function_space(
  std::string name) const
{
  _logger.trace("create a finite element map"_fmt);

  // @todo make factory for fem classes
  typename Traits::BaseFEM base_fem(_grid->leafGridView());

  std::shared_ptr<FEM> finite_element_map;
  if constexpr (Traits::is_sub_model)
    finite_element_map = std::make_shared<FEM>(_grid_view, base_fem);
  else
    finite_element_map = std::make_shared<FEM>(base_fem);

  _logger.trace("setup grid function space for component {}"_fmt, name);
  const ES entity_set(_grid->leafGridView());
  auto comp_gfs = std::make_shared<LGFS>(entity_set, finite_element_map);
  comp_gfs->name(name);

  return comp_gfs;
}

template<class Traits>
auto
ModelDiffusionReaction<Traits>::setup_domain_grid_function_space(
  std::vector<std::string> comp_names) const
{
  _logger.debug("setup domain grid function space"_fmt);

  typename GFS::NodeStorage comp_gfs_vec(comp_names.size());

  for (std::size_t k = 0; k < comp_names.size(); k++) {
    comp_gfs_vec[k] = setup_component_grid_function_space(comp_names[k]);
  }

  _logger.trace("setup power grid function space"_fmt);
  return std::make_shared<GFS>(comp_gfs_vec);
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_grid_function_space()
{
  // get component names
  auto& operator_splitting_config = _config.sub("operator");
  auto comp_names = operator_splitting_config.getValueKeys();
  std::sort(comp_names.begin(), comp_names.end());

  _states.clear();
  // map operator -> variable_name
  _operator_splitting.clear();
  for (auto&& var : comp_names) {
    std::size_t op = operator_splitting_config.template get<std::size_t>(var);
    _operator_splitting.insert(std::make_pair(op, var));
    _states[op] = {}; // initializate map with empty states
    _states[op].grid = _grid;
  }

  for (auto& [op, state] : _states) {
    _logger.trace("setup grid function space for operator {}"_fmt, op);
    auto op_range = _operator_splitting.equal_range(op);
    std::size_t op_size = std::distance(op_range.first, op_range.second);
    assert(op_size > 0);
    std::vector<std::string> op_comp_names(op_size);
    std::transform(op_range.first,
                   op_range.second,
                   op_comp_names.begin(),
                   [](const auto& i) { return i.second; });
    _states[op].grid_function_space =
      setup_domain_grid_function_space(op_comp_names);
  }
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_coefficient_vectors()
{
  for (auto& [op, state] : _states) {
    _logger.debug("setup coefficient vector for operator {}"_fmt, op);
    const auto& gfs = _states[op].grid_function_space;
    auto& x = _states[op].coefficients;
    if (x)
      x = std::make_shared<X>(*gfs, *(x->storage()));
    else
      x = std::make_shared<X>(*gfs);
  }
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_constraints()
{
  _logger.debug("setup constraints"_fmt);

  auto b0lambda = [&](const auto& i, const auto& x) { return false; };
  auto b0 =
    Dune::PDELab::makeBoundaryConditionFromCallable(_grid_view, b0lambda);

  _logger.trace("assemble constraints"_fmt);
  for (auto& [op, state] : _states) {
    _constraints[op] = std::make_unique<CC>();
    auto& gfs = _states[op].grid_function_space;
    Dune::PDELab::constraints(b0, *gfs, *_constraints[op]);

    _logger.info("constrained dofs in operator {}: {} of {}"_fmt,
                 op,
                 _constraints[op]->size(),
                 gfs->globalSize());
  }
}

template<class Traits>
auto
ModelDiffusionReaction<Traits>::setup_local_operator(std::size_t i) const
{
  _logger.trace("setup local operators {}"_fmt, i);

  _logger.trace("create spatial local operator {}"_fmt, i);
  typename Traits::BaseFEM::Traits::FiniteElement finite_element;

  auto local_operator =
    std::make_shared<LOP>(_grid_view, _config, finite_element, i);

  _logger.trace("create temporal local operator {}"_fmt, i);
  auto temporal_local_operator =
    std::make_shared<TLOP>(_grid_view, _config, finite_element, i);

  return std::make_pair(local_operator, temporal_local_operator);
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_local_operators()
{
  _logger.trace("setup local operators"_fmt);

  _local_operators.clear();
  _temporal_local_operators.clear();

  for (auto& [op, state] : _states) {
    auto operators = setup_local_operator(op);
    _local_operators[op] = operators.first;
    _temporal_local_operators[op] = operators.second;
  }
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_grid_operators()
{
  _logger.debug("setup grid operators"_fmt);
  _spatial_grid_operators.clear();
  _temporal_grid_operators.clear();
  _grid_operators.clear();

  for (auto& [op, state] : _states) {
    auto& gfs = _states[op].grid_function_space;
    auto& lop = _local_operators[op];
    auto& tlop = _temporal_local_operators[op];

    // @todo fix this estimate for something more accurate
    MBE mbe((int)Dune::power((int)3, (int)Grid::dimension));

    _logger.trace("create spatial grid operator {}"_fmt, op);
    _spatial_grid_operators[op] = std::make_shared<GOS>(
      *gfs, *_constraints[op], *gfs, *_constraints[op], *lop, mbe);

    _logger.trace("create temporal grid operator {}"_fmt, op);
    _temporal_grid_operators[op] = std::make_shared<GOT>(
      *gfs, *_constraints[op], *gfs, *_constraints[op], *tlop, mbe);

    _logger.trace("create instationary grid operator {}"_fmt, op);
    _grid_operators[op] = std::make_shared<GOI>(*_spatial_grid_operators[op],
                                                *_temporal_grid_operators[op]);
  }
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_solvers()
{
  _logger.debug("setup solvers"_fmt);
  _linear_solvers.clear();
  _nonlinear_solvers.clear();
  _time_stepping_methods.clear();
  _one_step_methods.clear();

  for (auto& [op, state] : _states) {
    auto& x = state.coefficients;
    _logger.trace("create linear solver"_fmt);
    _linear_solvers[op] = std::make_shared<LS>(*_grid_operators[op]);

    _logger.trace("create nonlinear solver"_fmt);
    _nonlinear_solvers[op] =
      std::make_shared<NLS>(*_grid_operators[op], *x, *_linear_solvers[op]);

    _logger.trace("select and prepare time-stepping scheme"_fmt);
    using AlexMethod = Dune::PDELab::Alexander2Parameter<double>;
    _time_stepping_methods[op] = std::make_shared<AlexMethod>();
    _one_step_methods[op] = std::make_shared<OSM>(*_time_stepping_methods[op],
                                                  *_grid_operators[op],
                                                  *_nonlinear_solvers[op]);
    _one_step_methods[op]->setVerbosityLevel(2);
  }
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_vtk_writer()
{
  _logger.debug("setup vtk writer"_fmt);

  _writer = std::make_shared<W>(_grid_view, Dune::VTK::conforming);
  auto config_writer = _config.sub("writer");
  std::string file_name = config_writer.template get<std::string>("file_name");
  std::string path = config_writer.get("path", file_name);
  struct stat st;

  if (stat(path.c_str(), &st) != 0) {
    int stat = 0;
#if defined(_WIN32)
    stat = mkdir(path.c_str());
#else
    stat = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#endif
    if (stat != 0 && stat != -1)
      std::cout << "Error: Cannot create directory " << path << std::endl;
  }

  _sequential_writer = std::make_shared<SW>(_writer, file_name, path, path);
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::suggest_timestep(double dt)
{
  // @todo do time addaptivity
  _config["time_step"] = dt;
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup(ModelSetupPolicy setup_policy)
{
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

template<class Traits>
void
ModelDiffusionReaction<Traits>::step()
{
  double dt = _config.template get<double>("time_step");

  const bool cout_redirected = Dune::Logging::Logging::isCoutRedirected();
  if (not cout_redirected)
    Dune::Logging::Logging::redirectCout(_solver_logger.name(),
                                         Dune::Logging::LogLevel::detail);

  // do time step
  for (auto& [op, state] : _states) {
    _local_operators[op]->update(const_states());

    auto& x = state.coefficients;
    auto x_new = std::make_shared<X>(*x);
    _one_step_methods[op]->apply(current_time(), dt, *x, *x_new);

    // accept time step
    x = x_new;
  }

  if (not cout_redirected)
    Dune::Logging::Logging::restoreCout();

  current_time() += dt;

  write_states();
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::write_states() const
{
  for (auto& [op, state] : _states) {
    auto& x = state.coefficients;
    auto& gfs = state.grid_function_space;

    auto etity_transformation = [&](auto e) {
      if constexpr (Concept::isMultiDomainGrid<Grid>() and
                    Concept::isSubDomainGrid<typename GV::Grid>())
        return _grid->multiDomainEntity(e);
      else
        return e;
    };
    using ET = decltype(etity_transformation);

    using Predicate = PDELab::vtk::DefaultPredicate;
    using Data = PDELab::vtk::DGFTreeCommonData<GFS, X, Predicate, GV, ET>;
    std::shared_ptr<Data> data =
      std::make_shared<Data>(*gfs, *x, _grid_view, etity_transformation);
    PDELab::vtk::OutputCollector<SW, Data> collector(*_sequential_writer, data);
    for (std::size_t k = 0; k < data->_lfs.degree(); k++)
      collector.addSolution(data->_lfs.child(k),
                            PDELab::vtk::defaultNameScheme());
  }
  _sequential_writer->write(current_time(), Dune::VTK::base64);
  _sequential_writer->vtkWriter()->clear();
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC