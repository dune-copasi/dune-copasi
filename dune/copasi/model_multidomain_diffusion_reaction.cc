#ifndef DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC
#define DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC

/**
 * This is the implementation of the ModelMultiDomainDiffusionReaction,
 * particularly, you want to notice that this is not an normal .cc
 * file but a header which has to be included when compiling.
 */

#include <dune/copasi/model_multidomain_diffusion_reaction.hh>

namespace Dune::Copasi {

template<class Grid, int FEMorder, class OrderingTag>
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::
  ModelMultiDomainDiffusionReaction(std::shared_ptr<Grid> grid,
                                    const Dune::ParameterTree& config)
  : ModelBase(config)
  , _solver_logger(Logging::Logging::componentLogger(config, "solver"))
  , _config(config)
  , _grid_view(grid->leafGridView())
  , _grid(grid)
  , _domains(config.sub("compartments").getValueKeys().size())
{
  setup();

  using GridFunction = ExpressionToGridFunctionAdapter<GridView, RF>;
  using CompartementGridFunction =
    PDELab::DynamicPowerGridFunction<GridFunction>;
  using MultiDomainGridFunction =
    PDELab::DynamicPowerGridFunction<CompartementGridFunction>;

  const auto& compartments = _config.sub("compartments").getValueKeys();

  std::vector<std::shared_ptr<CompartementGridFunction>> comp_functions(
    _domains);
  for (std::size_t i = 0; i < _domains; ++i) {

    std::size_t comp_size = _gfs->child(i).degree();
    _logger.trace("component size: {}"_fmt, comp_size);
    std::vector<std::shared_ptr<GridFunction>> functions(comp_size);

    const std::string compartement = compartments[i];
    auto& intial_config = _config.sub(compartement + ".initial");
    _logger.trace("initial size: {}"_fmt, intial_config.getValueKeys().size());
    assert(comp_size == intial_config.getValueKeys().size());

    for (std::size_t k = 0; k < comp_size; k++) {
      std::string var = intial_config.getValueKeys()[i];
      std::string eq = intial_config[var];
      functions[k] =
        std::make_shared<ExpressionToGridFunctionAdapter<GridView, RF>>(
          _grid_view, eq);
    }
    comp_functions[i] = std::make_shared<CompartementGridFunction>(functions);
  }

  MultiDomainGridFunction initial(comp_functions);
  Dune::PDELab::interpolate(initial, *_gfs, *_x);

  auto etity_transformation = [&](auto e) {
    return _grid->multiDomainEntity(e);
  };
  using ET = decltype(etity_transformation);
  using Predicate = PDELab::vtk::DefaultPredicate;
  using Data =
    PDELab::vtk::DGFTreeCommonData<GFS, X, Predicate, SubDomainGridView, ET>;

  for (std::size_t i = 0; i < _domains; ++i) {
    const std::string compartement = compartments[i];

    int sub_domain_id =
      _config.sub("compartments").template get<int>(compartement);
    SubDomainGridView sub_grid_view =
      _grid->subDomain(sub_domain_id).leafGridView();

    std::shared_ptr<Data> data =
      std::make_shared<Data>(*_gfs, *_x, sub_grid_view, etity_transformation);
    PDELab::vtk::OutputCollector<SW, Data> collector(*_sequential_writer[i],
                                                     data);
    collector.addSolution(data->_lfs.child(i),
                          PDELab::vtk::defaultNameScheme());

    _sequential_writer[i]->write(current_time(), Dune::VTK::appendedraw);
    _sequential_writer[i]->vtkWriter()->clear();
  }

  _logger.debug("ModelMultiDomainDiffusionReaction constructed"_fmt);
}

template<class Grid, int FEMorder, class OrderingTag>
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::
  ~ModelMultiDomainDiffusionReaction()
{
  _logger.debug("ModelMultiDomainDiffusionReaction deconstructed"_fmt);
}

template<class Grid, int FEMorder, class OrderingTag>
void
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::
  setup_grid_function_space()
{
  _logger.debug("setup grid function space"_fmt);

  const auto& compartments = _config.sub("compartments").getValueKeys();
  typename GFS::NodeStorage subdomain_gfs_vec(_domains);

  for (std::size_t i = 0; i < _domains; ++i) {
    const std::string compartement = compartments[i];
    auto& model_config = _config.sub(compartement);
    SubDomainGridView sub_grid_view = _grid->subDomain(i).leafGridView();

    _logger.trace("create a sub model for compartment {}"_fmt, i);
    auto sub_model = std::make_shared<SubModel>(
      _grid, model_config, sub_grid_view, ModelSetupPolicy::GridFunctionSpace);
    // extract grid function space from the sub model
    auto states = sub_model->states();
    // submodel should not be operator splitting!
    assert(states.size() == 1);
    subdomain_gfs_vec[i] = states[0].grid_function_space;
  }

  _logger.trace("setup power grid function space"_fmt);
  _gfs = std::make_shared<GFS>(subdomain_gfs_vec);
}

template<class Grid, int FEMorder, class OrderingTag>
void
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::
  setup_coefficient_vector()
{
  _logger.debug("setup coefficient vector"_fmt);
  if (_x)
    _x = std::make_shared<X>(*_gfs, *(_x->storage()));
  else
    _x = std::make_shared<X>(*_gfs);
}

template<class Grid, int FEMorder, class OrderingTag>
void
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::
  setup_constraints()
{
  _logger.debug("setup constraints"_fmt);

  auto b0lambda = [&](const auto& i, const auto& x) { return false; };
  auto b0 =
    Dune::PDELab::makeBoundaryConditionFromCallable(_grid_view, b0lambda);
  using B = Dune::PDELab::DynamicPowerConstraintsParameters<decltype(b0)>;
  std::vector<std::shared_ptr<decltype(b0)>> b0_vec;
  for (std::size_t i = 0; i < _domains; i++)
    b0_vec.emplace_back(std::make_shared<decltype(b0)>(b0));
  B b(b0_vec);

  _logger.trace("assemble constraints"_fmt);
  _constraints = std::make_unique<CC>();
  Dune::PDELab::constraints(b, *_gfs, *_constraints);

  _logger.info(
    "constrained dofs: {} of {}"_fmt, _constraints->size(), _gfs->globalSize());
}

template<class Grid, int FEMorder, class OrderingTag>
void
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::
  setup_local_operators()
{
  _logger.trace("setup local operators"_fmt);

  _logger.trace("create spatial local operator"_fmt);
  FE finite_element;
  _local_operator = std::make_shared<LOP>(_grid, _config, finite_element);

  _logger.trace("create temporal local operator"_fmt);
  _temporal_local_operator =
    std::make_shared<TLOP>(_grid, _config, finite_element);
}

template<class Grid, int FEMorder, class OrderingTag>
void
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::
  setup_grid_operators()
{
  _logger.debug("setup grid operators"_fmt);
  MBE mbe((int)pow(3, dim));

  _logger.trace("create spatial grid operator"_fmt);
  _spatial_grid_operator = std::make_shared<GOS>(
    *_gfs, *_constraints, *_gfs, *_constraints, *_local_operator, mbe);

  _logger.trace("create temporal grid operator"_fmt);
  _temporal_grid_operator = std::make_shared<GOT>(
    *_gfs, *_constraints, *_gfs, *_constraints, *_temporal_local_operator, mbe);

  _logger.trace("create instationary grid operator"_fmt);
  _grid_operator =
    std::make_shared<GOI>(*_spatial_grid_operator, *_temporal_grid_operator);
}

template<class Grid, int FEMorder, class OrderingTag>
void
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::setup_solvers()
{
  _logger.debug("setup solvers"_fmt);

  _logger.trace("create linear solver"_fmt);
  _linear_solver = std::make_shared<LS>(*_grid_operator);

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

template<class Grid, int FEMorder, class OrderingTag>
void
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::
  setup_vtk_writer()
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
      stat = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
      if (stat != 0 && stat != -1)
        std::cout << "Error: Cannot create directory " << path << std::endl;
    }

    _sequential_writer[i] = std::make_shared<SW>(writer, file_name, path, path);
  }
}

template<class Grid, int FEMorder, class OrderingTag>
void
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::
  suggest_timestep(double dt)
{
  DUNE_THROW(NotImplemented, "not implemented");
}

template<class Grid, int FEMorder, class OrderingTag>
void
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::setup(
  ModelSetupPolicy setup_policy)
{
  _logger.trace("setup operator started"_fmt);

  if (setup_policy >= ModelSetupPolicy::GridFunctionSpace)
    setup_grid_function_space();
  if (setup_policy >= ModelSetupPolicy::CoefficientVector)
    setup_coefficient_vector();
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

template<class Grid, int FEMorder, class OrderingTag>
void
ModelMultiDomainDiffusionReaction<Grid, FEMorder, OrderingTag>::step()
{
  double dt = _config.template get<double>("time_step");

  const bool cout_redirected = Dune::Logging::Logging::isCoutRedirected();
  if (not cout_redirected)
    Dune::Logging::Logging::redirectCout(_solver_logger.name(),
                                         Dune::Logging::LogLevel::detail);

  // do time step
  auto x_new = std::make_shared<X>(*_x);
  _one_step_method->apply(current_time(), dt, *_x, *x_new);

  if (not cout_redirected)
    Dune::Logging::Logging::restoreCout();

  // accept time step
  _x = x_new;
  current_time() += dt;

  auto etity_transformation = [&](auto e) {
    return _grid->multiDomainEntity(e);
  };
  using ET = decltype(etity_transformation);
  using Predicate = PDELab::vtk::DefaultPredicate;
  using Data =
    PDELab::vtk::DGFTreeCommonData<GFS, X, Predicate, SubDomainGridView, ET>;

  const auto& compartments = _config.sub("compartments").getValueKeys();
  for (std::size_t i = 0; i < _domains; ++i) {
    const std::string compartement = compartments[i];

    int sub_domain_id =
      _config.sub("compartments").template get<int>(compartement);
    SubDomainGridView sub_grid_view =
      _grid->subDomain(sub_domain_id).leafGridView();

    std::shared_ptr<Data> data =
      std::make_shared<Data>(*_gfs, *_x, sub_grid_view, etity_transformation);
    PDELab::vtk::OutputCollector<SW, Data> collector(*_sequential_writer[i],
                                                     data);
    collector.addSolution(data->_lfs.child(i),
                          PDELab::vtk::defaultNameScheme());

    _sequential_writer[i]->write(current_time(), Dune::VTK::appendedraw);
    _sequential_writer[i]->vtkWriter()->clear();
  }
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC