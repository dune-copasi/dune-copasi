#ifndef DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC
#define DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC

/**
 * This is the implementation of the ModelMultiDomainDiffusionReaction,
 * particularly, you want to notice that this is not an normal .cc
 * file but a header which has to be included when compiling.
 */

#include <dune/copasi/common/muparser_data_handler.hh>
#include <dune/copasi/model/multidomain_diffusion_reaction.hh>

#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>

namespace Dune::Copasi {

template<class Traits>
ModelMultiDomainDiffusionReaction<Traits>::ModelMultiDomainDiffusionReaction(
  std::shared_ptr<Grid> grid,
  const Dune::ParameterTree& config,
  ModelSetupPolicy setup_policy)
  : ModelBase(config)
  , _solver_logger(Logging::Logging::componentLogger(config, "solver"))
  , _config(config)
  , _grid_view(grid->leafGridView())
  , _grid(grid)
  , _domains(config.sub("compartments").getValueKeys().size())
{
  setup(setup_policy);

  if (setup_policy == ModelSetupPolicy::All) {

    MuParserDataHandler<TIFFGrayscale<unsigned short>> mu_data_handler;
    if (_config.hasSub("data"))
      mu_data_handler.add_tiff_functions(_config.sub("data"));

    auto initial_muparser =
      get_muparser_initial(_config, _grid->leafGridView(), false);

    for (auto&& sd_grid_function : initial_muparser) {
      for (auto&& mu_grid_function : sd_grid_function) {
        mu_data_handler.set_functions(mu_grid_function.parser());
        mu_grid_function.compile_parser();
        mu_grid_function.set_time(current_time());
      }
    }
    set_initial(initial_muparser);

    update_data_handler();
    write_states();
  }

  _logger.debug("ModelMultiDomainDiffusionReaction constructed"_fmt);
}

template<class Traits>
ModelMultiDomainDiffusionReaction<Traits>::~ModelMultiDomainDiffusionReaction()
{
  _logger.debug("ModelMultiDomainDiffusionReaction deconstructed"_fmt);
}

template<class Traits>
template<class GFGridView>
auto
ModelMultiDomainDiffusionReaction<Traits>::get_muparser_initial(
  const ParameterTree& model_config,
  const GFGridView& gf_grid_view,
  bool compile)
{
  const auto& compartments = model_config.sub("compartments").getValueKeys();

  using GridFunction = ExpressionToGridFunctionAdapter<GV, RF>;
  std::vector<std::vector<std::shared_ptr<GridFunction>>> functions;

  for (std::size_t domain = 0; domain < compartments.size(); domain++) {
    const std::string compartement = compartments[domain];
    auto sub_model_config = model_config.sub(compartement);
    auto sub_model_initial =
      SubModel::get_muparser_initial(sub_model_config, gf_grid_view, compile);
    functions.emplace_back(sub_model_initial);
  }

  return functions;
}

template<class Traits>
template<class GF>
void
ModelMultiDomainDiffusionReaction<Traits>::set_initial(
  const std::vector<std::vector<std::shared_ptr<GF>>>& initial)
{
  using GridFunction = std::decay_t<GF>;
  static_assert(Concept::isPDELabGridFunction<GridFunction>(),
                "GridFunction is not a PDElab grid functions");
  static_assert(
    std::is_same_v<typename GridFunction::Traits::GridViewType, GV>,
    "GridFunction has to have the same grid view as the templated grid");
  static_assert((int)GridFunction::Traits::dimDomain == (int)Grid::dimension,
                "GridFunction has to have domain dimension equal to the grid");
  static_assert(GridFunction::Traits::dimRange == 1,
                "GridFunction has to have range dimension equal to 1");

  _logger.debug("set initial state from grid functions"_fmt);

  const auto& compartments_config = _config.sub("compartments");
  if (initial.size() != compartments_config.getValueKeys().size())
    DUNE_THROW(RangeError, "Wrong number of grid functions");

  using CompartementGridFunction =
    PDELab::DynamicPowerGridFunction<GridFunction>;
  using MultiDomainGridFunction =
    PDELab::DynamicPowerGridFunction<CompartementGridFunction>;

  for (auto& [op, state] : _states) {
    _logger.trace("interpolation of operator {}"_fmt, op);

    std::vector<std::shared_ptr<CompartementGridFunction>> md_functions(
      _domains);

    for (std::size_t i = 0; i < initial.size(); ++i) {
      std::size_t comp_size = state.grid_function_space->child(i).degree();
      std::vector<std::shared_ptr<GridFunction>> sd_functions(comp_size);

      auto compartment = compartments_config.getValueKeys()[i];
      auto& sub_model_config = _config.sub(compartment);
      auto& operator_config = sub_model_config.sub("operator");
      auto comp_names = operator_config.getValueKeys();
      std::sort(comp_names.begin(), comp_names.end());

      std::size_t count_i = 0; // index for variables in operator
      std::size_t count_j = 0; // index for all variables
      for (const auto& var : comp_names) {
        assert(count_i < comp_names.size());
        if (op != operator_config.template get<std::size_t>(var))
          continue;

        _logger.trace("creating grid function for variable: {}"_fmt, var);
        sd_functions[count_i] = initial[i][count_j];
        sd_functions[count_i]->set_time(current_time());
        count_i++;
        count_j++;
      }

      md_functions[i] =
        std::make_shared<CompartementGridFunction>(sd_functions);
    }

    MultiDomainGridFunction md_initial(md_functions);
    Dune::PDELab::interpolate(
      md_initial, *state.grid_function_space, *state.coefficients);
  }
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_grid_function_spaces()
{
  _logger.debug("setup grid function space"_fmt);

  const auto& compartments = _config.sub("compartments").getValueKeys();
  std::map<std::size_t, typename GFS::NodeStorage> gfs_operator_vec;

  _states.clear();

  for (std::size_t domain = 0; domain < _domains; ++domain) {
    const std::string compartement = compartments[domain];
    auto model_config = _config.sub(compartement);
    model_config["begin_time"] = _config["begin_time"];
    model_config["end_time"] = _config["end_time"];
    model_config["time_step"] = _config["time_step"];
    SubDomainGridView sub_grid_view = _grid->subDomain(domain).leafGridView();

    _logger.trace("create a sub model for compartment {}"_fmt, domain);
    auto sub_model = std::make_shared<SubModel>(
      _grid, model_config, sub_grid_view, ModelSetupPolicy::GridFunctionSpace);
    // extract grid function space from the sub model
    auto states = sub_model->states();
    for (auto& [op, state] : states) {
      gfs_operator_vec[op].resize(_domains);
      gfs_operator_vec[op][domain] = state.grid_function_space;
    }
  }

  for (auto& [op, gfs_vector] : gfs_operator_vec) {
    _logger.trace("setup power grid function space for operator {}"_fmt, op);

    for (std::size_t domain = 0; domain < _domains; domain++) {
      if (not gfs_vector[domain]) {
        DUNE_THROW(RangeError,
                   "Operators have to have at least one variable per domain!");
      }
    }

    assert(gfs_vector.size() == _domains);
    _states[op].grid_function_space = std::make_shared<GFS>(gfs_vector);
  }
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_coefficient_vectors()
{
  for (auto& [op, state] : _states) {
    _logger.debug("setup coefficient vector for operator {}"_fmt, op);
    const auto& gfs = state.grid_function_space;
    auto& x = state.coefficients;
    if (x)
      x = std::make_shared<X>(*gfs, *(x->storage()));
    else
      x = std::make_shared<X>(*gfs);
  }
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_constraints()
{
  _logger.debug("setup base constraints"_fmt);

  auto b0lambda = [&](const auto& i, const auto& x) { return false; };
  auto b0 =
    Dune::PDELab::makeBoundaryConditionFromCallable(_grid_view, b0lambda);

  _logger.trace("setup power constraints"_fmt);
  using B = Dune::PDELab::DynamicPowerConstraintsParameters<decltype(b0)>;
  std::vector<std::shared_ptr<decltype(b0)>> b0_vec;
  for (std::size_t i = 0; i < _domains; i++)
    b0_vec.emplace_back(std::make_shared<decltype(b0)>(b0));
  B b(b0_vec);

  _logger.trace("assemble constraints"_fmt);
  for (auto& [op, state] : _states) {
    const auto& gfs = state.grid_function_space;
    _constraints[op] = std::make_unique<CC>();
    Dune::PDELab::constraints(b, *gfs, *_constraints[op]);

    _logger.info("constrained dofs in operator {}: {} of {}"_fmt,
                 op,
                 _constraints[op]->size(),
                 gfs->globalSize());
  }
}

template<class Traits>
auto
ModelMultiDomainDiffusionReaction<Traits>::setup_local_operator(
  std::size_t i) const
{
  _logger.trace("setup local operators {}"_fmt, i);

  _logger.trace("create spatial local operator {}"_fmt, i);
  typename Traits::SubModelTraits::BaseFEM::Traits::FiniteElement
    finite_element;

  auto local_operator =
    std::make_shared<LOP>(_grid, _config, finite_element, i);

  _logger.trace("create temporal local operator {}"_fmt, i);
  auto temporal_local_operator =
    std::make_shared<TLOP>(_grid, _config, finite_element, i);

  return std::make_pair(local_operator, temporal_local_operator);
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_local_operators()
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
ModelMultiDomainDiffusionReaction<Traits>::setup_grid_operators()
{
  _logger.debug("setup grid operators"_fmt);
  _spatial_grid_operators.clear();
  _temporal_grid_operators.clear();
  _grid_operators.clear();

  for (auto& [op, state] : _states) {
    auto& gfs = _states[op].grid_function_space;
    auto& lop = _local_operators[op];
    auto& tlop = _temporal_local_operators[op];

    std::size_t max_comps(0);
    for (std::size_t i = 0; i < gfs->degree(); i++)
      max_comps = std::max(max_comps, gfs->child(i).degree());

    // @todo fix this estimate for something more accurate
    MBE mbe((int)Dune::power((int)3, (int)Grid::dimension) * max_comps);

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
ModelMultiDomainDiffusionReaction<Traits>::setup_solvers()
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
ModelMultiDomainDiffusionReaction<Traits>::setup_vtk_writer()
{
  _logger.debug("setup vtk writer"_fmt);
  const auto& compartments = _config.sub("compartments").getValueKeys();

  _sequential_writer.resize(_domains);

  for (std::size_t i = 0; i < _domains; ++i) {

    const std::string compartement = compartments[i];
    auto& model_config = _config.sub(compartement);

    int sub_domain_id =
      _config.sub("compartments").template get<int>(compartement);
    SubDomainGridView sub_grid_view =
      _grid->subDomain(sub_domain_id).leafGridView();

    std::shared_ptr<W> writer =
      std::make_shared<W>(sub_grid_view, Dune::VTK::conforming);

    auto config_writer = model_config.sub("writer");
    std::string file_name =
      config_writer.template get<std::string>("file_name");
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

    _sequential_writer[i] = std::make_shared<SW>(writer, file_name, path, "");
  }
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::suggest_timestep(double dt)
{
  // @todo do time addaptivity
  _config["time_step"] = dt;
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup(ModelSetupPolicy setup_policy)
{
  _logger.trace("setup operator started"_fmt);

  if (setup_policy >= ModelSetupPolicy::GridFunctionSpace)
    setup_grid_function_spaces();
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
ModelMultiDomainDiffusionReaction<Traits>::step()
{
  double dt = _config.template get<double>("time_step");

  const bool cout_redirected = Dune::Logging::Logging::isCoutRedirected();
  if (not cout_redirected)
    Dune::Logging::Logging::redirectCout(_solver_logger.name(),
                                         Dune::Logging::LogLevel::detail);

  double max_error = std::numeric_limits<double>::max();
  auto states_after = _states;

  _logger.info("Time Step {:.2e} + {:.2e} -> {:.2e}"_fmt,
               current_time(),
               dt,
               current_time() + dt);
  std::size_t op_iter = 0;
  do {
    const auto states_before = states_after;
    _solver_logger.notice("Iteration {}"_fmt, op_iter);
    for (auto& [op, state] : states_after) {
      _solver_logger.notice("Operator {}"_fmt, op);
      _local_operators[op]->update(const_states(states_after));

      auto& x = state.coefficients;
      auto x_new = std::make_shared<X>(*x);
      _one_step_methods[op]->apply(
        current_time(), dt, *_states[op].coefficients, *x_new);
      x = x_new;
    }

    if (op_iter > 0) {
      max_error = 0.;
      const auto& compartments = _config.sub("compartments").getValueKeys();
      for (std::size_t i = 0; i < _domains; i++) {
        const auto& comp_names =
          _config.sub(compartments[i] + ".initial").getValueKeys();
        for (std::size_t j = 0; j < comp_names.size(); j++) {
          auto gf_before = get_grid_function(states_before, i, j);
          auto gf_after = get_grid_function(states_after, i, j);
          PDELab::DifferenceSquaredAdapter<ComponentGridFunction,
                                           ComponentGridFunction>
            err(gf_before, gf_after);
          FieldVector<double, 1> split_error;
          PDELab::integrateGridFunction(err, split_error, 1);
          split_error = err.getGridView().comm().sum(std::sqrt(split_error));
          _solver_logger.debug(
            "error in domain {}, component {} is: {}"_fmt, i, j, split_error);
          max_error = std::max(max_error, split_error[0]);
        }
      }
    }
    op_iter++;
  } while (max_error >= 1e-14 and _states.size() > 1);

  // @todo integrate each component and calculate error after iteration

  if (not cout_redirected)
    Dune::Logging::Logging::restoreCout();

  current_time() += dt;

  // update to new states
  _states = states_after;

  update_data_handler();
  write_states();
}

template<class Traits>
auto
ModelMultiDomainDiffusionReaction<Traits>::get_data_handler(
  std::map<std::size_t, State> states) const
{
  std::vector<std::map<std::size_t, std::shared_ptr<DataHandler>>> data(
    _domains);

  const auto& compartments = _config.sub("compartments").getValueKeys();
  EntityTransformation et(_grid);
  for (std::size_t i = 0; i < _domains; ++i) {
    const std::string compartement = compartments[i];
    int sub_domain_id =
      _config.sub("compartments").template get<int>(compartement);
    SubDomainGridView sub_grid_view =
      _grid->subDomain(sub_domain_id).leafGridView();
    for (auto& [op, state] : states) {
      data[i][op] = std::make_shared<DataHandler>(
        *state.grid_function_space, *state.coefficients, sub_grid_view, et);
    }
  }
  return data;
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::update_data_handler()
{
  _data = get_data_handler(_states);
}

template<class Traits>
auto
ModelMultiDomainDiffusionReaction<Traits>::get_grid_function(
  const std::map<std::size_t, State>& states,
  std::size_t domain,
  std::size_t comp) const
{
  auto data = get_data_handler(states);
  const auto& compartments = _config.sub("compartments").getValueKeys();
  const std::string compartement = compartments[domain];
  auto& operator_config = _config.sub(compartement + ".operator");
  auto op_keys = operator_config.getValueKeys();
  std::sort(op_keys.begin(), op_keys.end());
  std::size_t op = operator_config.template get<std::size_t>(op_keys[comp]);
  std::size_t gfs_comp = 0;

  for (std::size_t i = 0; i < comp; i++)
    if (op == operator_config.template get<std::size_t>(op_keys[i]))
      gfs_comp++;

  const auto& data_comp = data[domain].at(op);
  ComponentGridFunction gf(data_comp->_lfs.child(domain).child(gfs_comp),
                           data_comp);
  return gf;
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::write_states() const
{

  const auto& compartments = _config.sub("compartments").getValueKeys();

  for (std::size_t i = 0; i < _domains; ++i) {
    const std::string compartement = compartments[i];

    for (auto& [op, state] : _states) {
      const auto& data = _data[i].at(op);
      PDELab::vtk::OutputCollector<SW, DataHandler> collector(
        *_sequential_writer[i], data);
      for (std::size_t k = 0; k < data->_lfs.child(i).degree(); k++) {
        collector.addSolution(data->_lfs.child(i).child(k),
                              PDELab::vtk::defaultNameScheme());
      }
    }
    _sequential_writer[i]->write(current_time(), Dune::VTK::base64);
    _sequential_writer[i]->vtkWriter()->clear();
  }
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC