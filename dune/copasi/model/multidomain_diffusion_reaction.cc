#ifndef DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC
#define DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC

/**
 * This is the implementation of the ModelMultiDomainDiffusionReaction,
 * particularly, you want to notice that this is not an normal .cc
 * file but a header which has to be included when compiling.
 */

#include <dune/copasi/common/bit_flags.hh>

#include <dune/copasi/common/muparser_data_handler.hh>
#include <dune/copasi/model/multidomain_diffusion_reaction.hh>

#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>

namespace Dune::Copasi {

template<class Traits>
ModelMultiDomainDiffusionReaction<Traits>::ModelMultiDomainDiffusionReaction(
  std::shared_ptr<Grid> grid,
  const Dune::ParameterTree& config,
  BitFlags<ModelSetup::Stages> setup_policy)
  : ModelBase(config)
  , _config(config)
  , _grid_view(grid->leafGridView())
  , _grid(grid)
  , _domains(config.sub("compartments",true).getValueKeys().size())
{
  setup(setup_policy);

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
  const auto& compartments = model_config.sub("compartments",true).getValueKeys();

  using GridFunction = ExpressionToGridFunctionAdapter<GV, RF>;
  std::vector<std::vector<std::shared_ptr<GridFunction>>> functions;

  for (std::size_t domain = 0; domain < compartments.size(); domain++) {
    const std::string compartement = compartments[domain];
    auto sub_model_config = model_config.sub(compartement,true);
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

  _logger.debug("Set initial state from grid functions"_fmt);

  const auto& compartments_config = _config.sub("compartments",true);
  if (initial.size() != compartments_config.getValueKeys().size())
    DUNE_THROW(RangeError, "Wrong number of grid functions");

  using CompartementGridFunction =
    PDELab::DynamicPowerGridFunction<GridFunction>;
  using MultiDomainGridFunction =
    PDELab::DynamicPowerGridFunction<CompartementGridFunction>;

  std::vector<std::shared_ptr<CompartementGridFunction>> md_functions(_domains);

  for (std::size_t i = 0; i < initial.size(); ++i)
    md_functions[i] = std::make_shared<CompartementGridFunction>(initial[i]);

  MultiDomainGridFunction md_initial(md_functions);
  Dune::PDELab::interpolate(
    md_initial, *_state.grid_function_space, *_state.coefficients);
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_grid_function_space()
{
  _logger.debug("Setup grid function space"_fmt);

  const auto& compartments = _config.sub("compartments",true).getValueKeys();
  auto gfs_vec = typename GFS::NodeStorage{};
  gfs_vec.resize(_domains);

  if (not _state) {
    _state.grid = _grid;
    _state.time = _config.get("time_stepping.begin",0.);
  }

  for (std::size_t domain = 0; domain < _domains; ++domain) {
    const std::string compartement = compartments[domain];

    // remove other compartments from sub model config, this will make the sub
    // model to only see one compartment tree
    auto sub_model_config = _config;
    sub_model_config.sub("compartments") = ParameterTree{};
    auto compartment_sect = "compartments." + compartement;
    sub_model_config[compartment_sect] = _config[compartment_sect];
    SubDomainGridView sub_grid_view = _grid->subDomain(domain).leafGridView();

    _logger.trace("Create a sub model for compartment {}"_fmt, domain);
    auto sub_model =
      std::make_shared<SubModel>(_grid,
                                 sub_model_config,
                                 sub_grid_view,
                                 ModelSetup::setup_grid_function_space);
    // extract grid function space from the sub model
    auto sub_state = sub_model->state();
    gfs_vec[domain] = sub_state.grid_function_space;
  }

  _state.grid_function_space = std::make_shared<GFS>(gfs_vec);
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_coefficient_vector()
{
  _logger.debug("Setup coefficient vector"_fmt);
  const auto& gfs = _state.grid_function_space;
  if (not gfs)
    DUNE_THROW(InvalidStateException,"Grid function space is not setup");
  _state.coefficients = std::make_shared<X>(*gfs); // create new storage
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_initial_condition()
{
  // get TIFF data if available
  MuParserDataHandler<TIFFGrayscale<unsigned short>> mu_data_handler;
  if (_config.hasSub("data"))
    mu_data_handler.add_tiff_functions(_config.sub("data"));

  // configure math parsers for initial conditions on each component
  auto initial_muparser =
    get_muparser_initial(_config, _grid->leafGridView(), false);

  for (auto&& sd_grid_function : initial_muparser) {
    for (auto&& mu_grid_function : sd_grid_function) {
      // make TIFF expression available in the parser
      mu_data_handler.set_functions(mu_grid_function->parser());
      mu_grid_function->compile_parser();
      mu_grid_function->set_time(_state.time);
    }
  }
  // evaluate compiled expressions as initial conditions
  set_initial(initial_muparser);
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_constraints()
{
  _logger.debug("Setup base constraints"_fmt);

  auto b0lambda = [&](const auto& i, const auto& x) { return false; };
  auto b0 =
    Dune::PDELab::makeBoundaryConditionFromCallable(_grid_view, b0lambda);

  _logger.trace("Setup power constraints"_fmt);
  using B = Dune::PDELab::DynamicPowerConstraintsParameters<decltype(b0)>;
  std::vector<std::shared_ptr<decltype(b0)>> b0_vec;
  for (std::size_t i = 0; i < _domains; i++)
    b0_vec.emplace_back(std::make_shared<decltype(b0)>(b0));
  B b(b0_vec);

  _logger.trace("Assemble constraints"_fmt);
  const auto& gfs = _state.grid_function_space;
  _constraints = std::make_unique<CC>();
  Dune::PDELab::constraints(b, *gfs, *_constraints);

  _logger.detail(
    "Constrained dofs: {} of {}"_fmt, _constraints->size(), gfs->globalSize());
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_local_operator()
{
  _logger.debug("Setup local operator"_fmt);

  _logger.trace("Create spatial local operator"_fmt);

  _local_operator = std::make_shared<LOP>(_grid, _config);

  _logger.trace("Create temporal local operator"_fmt);
  _temporal_local_operator = std::make_shared<TLOP>(_grid, _config);
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_grid_operator()
{
  _logger.debug("Setup grid operator"_fmt);

  auto& gfs = _state.grid_function_space;
  auto& lop = _local_operator;
  auto& tlop = _temporal_local_operator;

  std::size_t max_comps(0);
  for (std::size_t i = 0; i < gfs->degree(); i++)
    max_comps = std::max(max_comps, gfs->child(i).degree());

  // @todo fix this estimate for something more accurate
  MBE mbe((int)Dune::power((int)3, (int)Grid::dimension) * max_comps);

  _logger.trace("Create spatial grid operator"_fmt);
  _spatial_grid_operator = std::make_shared<GOS>(
    *gfs, *_constraints, *gfs, *_constraints, *lop, mbe);

  _logger.trace("Create temporal grid operator"_fmt);
  _temporal_grid_operator = std::make_shared<GOT>(
    *gfs, *_constraints, *gfs, *_constraints, *tlop, mbe);

  _logger.trace("Create instationary grid operator"_fmt);
  _grid_operator = std::make_shared<GOI>(*_spatial_grid_operator,
                                              *_temporal_grid_operator);
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_jacobian_operator()
{
  _logger.debug("Create jacobian operator"_fmt);
  _jacobian_operator = std::make_shared<JacobianOperator>(*_grid_operator);
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup_vtk_writer()
{
  _logger.debug("Setup VTK writer"_fmt);
  // this map contains the written timesteps for a given path and is shared
  // among different `state`s copied from this one.
  auto mapu_out =
    std::make_shared<std::map<std::string, std::vector<double>>>();

  _state.writer = [=](
                    const auto& state, const auto& path, bool append) mutable {
    auto log_writer = Logging::Logging::componentLogger({}, "writer");
    auto& gfs = state.grid_function_space;

    // create directory if necessary
    auto path_entry = std::filesystem::directory_entry{ path };
    if (not path_entry.exists()) {
      log_writer.info("Creating output directory '{}'"_fmt,
                      path_entry.path().string());
      std::error_code ec{ 0, std::generic_category() };
      std::filesystem::create_directories(path_entry.path(), ec);
      if (ec)
        DUNE_THROW(IOError,
                   "\n Category: " << ec.category().name() << '\n'
                                   << "Value: " << ec.value() << '\n'
                                   << "Message: " << ec.message() << '\n');
    }

    // Recover old timestesps in case something was written before
    auto& timesteps = (*mapu_out)[path.string()];

    auto data = get_data_handler(state);

    // let's avoid throubles and write a different file for each sub-domain
    // https://public.kitware.com/pipermail/paraview/2014-July/031732.html

    for (size_t i = 0; i < data.size(); i++) {
      std::string name = fmt::format("{}-{}", path.filename().string(), gfs->child(i).name());
      if (not append) {
        timesteps.clear();
        log_writer.detail("Creating a time sequence file: '{}.pvd'"_fmt, name);
      } else {
        log_writer.trace("Overriding time sequence file: '{}.pvd'"_fmt, name);
      }

      // setup writer again with old timesteps if necessary
      SubDomainGridView sub_grid_view =
        gfs->gridView().grid().subDomain(i).leafGridView();
      auto writer = std::make_shared<VTKWriter<SubDomainGridView>>(
        sub_grid_view, Dune::VTK::conforming);

      auto sequential_writer = VTKSequenceWriter{ writer, name, path, path };
      sequential_writer.setTimeSteps(timesteps);

      auto collector =
        PDELab::vtk::OutputCollector{ sequential_writer, data[i] };
      for (std::size_t k = 0; k < data[i]->_lfs.child(i).degree(); k++) {
        collector.addSolution(data[i]->_lfs.child(i).child(k),
                              PDELab::vtk::defaultNameScheme());
      }
      sequential_writer.write(state.time, Dune::VTK::base64);
      sequential_writer.vtkWriter()->clear();
    }
    timesteps.push_back(state.time);
  };
}

template<class Traits>
void
ModelMultiDomainDiffusionReaction<Traits>::setup(
  BitFlags<ModelSetup::Stages> setup_policy)
{
  _logger.info("Setting up multi-compartment diffusion-reaction model"_fmt);

  try
  {
    if (setup_policy.test(ModelSetup::Stages::GridFunctionSpace))
      setup_grid_function_space();
    if (setup_policy.test(ModelSetup::Stages::CoefficientVector))
      setup_coefficient_vector();
    if (setup_policy.test(ModelSetup::Stages::InitialCondition))
      setup_initial_condition();
    if (setup_policy.test(ModelSetup::Stages::Constraints))
      setup_constraints();
    if (setup_policy.test(ModelSetup::Stages::LocalOperator))
      setup_local_operator();
    if (setup_policy.test(ModelSetup::Stages::GridOperator))
      setup_grid_operator();
    if (setup_policy.test(ModelSetup::Stages::JacobianOperator))
      setup_jacobian_operator();
    if (setup_policy.test(ModelSetup::Stages::Writer))
      setup_vtk_writer();
  }
  catch(...)
  {
    // detailed report of the configuration paramter tree
    std::stringstream ss;
    _config.report(ss);
    _logger.error(
      2, "---- Multidomain Diffusion-Reaction model parameter tree"_fmt);
    for (std::string line; std::getline(ss, line);)
      _logger.error(2, "{}"_fmt, line);
    _logger.error(2, "----"_fmt);
    throw;
  }
}

template<class Traits>
auto
ModelMultiDomainDiffusionReaction<Traits>::get_data_handler(
  const ConstState& state) const
{
  std::vector<std::shared_ptr<DataHandler>> data(_domains);

  const auto& compartments = _config.sub("compartments",true).getValueKeys();
  EntityTransformation et(_grid);
  for (std::size_t i = 0; i < _domains; ++i) {
    const std::string compartement = compartments[i];
    int sub_domain_id =
      _config.sub("compartments",true).template get<int>(compartement);
    SubDomainGridView sub_grid_view =
      _grid->subDomain(sub_domain_id).leafGridView();
    data[i] = std::make_shared<DataHandler>(
      *state.grid_function_space, *state.coefficients, sub_grid_view, et);
  }
  return data;
}

template<class Traits>
auto
ModelMultiDomainDiffusionReaction<Traits>::get_grid_function(
  const ConstState& state,
  std::size_t domain,
  std::size_t comp) const -> std::shared_ptr<ComponentGridFunction>
{
  auto data = get_data_handler(state);
  return std::make_shared<ComponentGridFunction>(
    data[domain]->_lfs.child(domain).child(comp), data[domain]);
}

template<class Traits>
auto
ModelMultiDomainDiffusionReaction<Traits>::get_grid_function(
  std::size_t domain,
  std::size_t comp) const -> std::shared_ptr<ComponentGridFunction>
{
  return get_grid_function(const_state(), domain, comp);
}

template<class Traits>
auto
ModelMultiDomainDiffusionReaction<Traits>::get_grid_functions(
  const ConstState& state) const
  -> std::vector<std::vector<std::shared_ptr<ComponentGridFunction>>>
{
  const auto& compartments = _config.sub("compartments",true).getValueKeys();
  std::vector<std::vector<std::shared_ptr<ComponentGridFunction>>>
    grid_functions(_domains);

  for (std::size_t domain_i = 0; domain_i < _domains; domain_i++) {
    const std::string compartement = compartments[domain_i];
    std::size_t domain =
      _config.sub("compartments",true).template get<std::size_t>(compartement);
    const auto& model_config = _config.sub(compartments[domain_i],true);
    const auto& vars = model_config.sub("diffusion",true);
    grid_functions[domain_i].resize(vars.getValueKeys().size());
    for (std::size_t var_i = 0; var_i < vars.getValueKeys().size(); var_i++) {
      grid_functions[domain_i][var_i] =
        get_grid_function(state, domain, var_i);
    }
  }
  return grid_functions;
}

template<class Traits>
auto
ModelMultiDomainDiffusionReaction<Traits>::get_grid_functions() const
  -> std::vector<std::vector<std::shared_ptr<ComponentGridFunction>>>
{
  return get_grid_functions(const_state());
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_CC
