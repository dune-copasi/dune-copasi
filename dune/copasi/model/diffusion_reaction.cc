#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC

/**
 * This is the implementation of the ModelDiffusionReaction,
 * particularly, you want to notice that this is not an normal .cc
 * file but a header which has to be included when compiling.
 */

#include <dune/copasi/common/bit_flags.hh>
#include <dune/copasi/common/data_context.hh>
#include <dune/copasi/common/muparser_data_handler.hh>
#include <dune/copasi/common/pdelab_expression_adapter.hh>
#include <dune/copasi/common/tiff_grayscale.hh>
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
  BitFlags<ModelSetup::Stages> setup_policy)
  : ModelBase(config)
  , _config(config)
  , _compartment_name(_config.sub("compartments").getValueKeys().front())
  , _grid_view(grid_view)
  , _grid(grid)
{
  if (_config.sub("compartments",true).getValueKeys().size() != 1)
    DUNE_THROW(IOError, "'compartments' section must contain one entry");
  setup(setup_policy);

  _logger.trace("ModelDiffusionReaction constructed"_fmt);
}

template<class Traits>
ModelDiffusionReaction<Traits>::~ModelDiffusionReaction()
{
  _logger.trace("ModelDiffusionReaction deconstructed"_fmt);
}

template<class Traits>
template<class GFGridView>
auto
ModelDiffusionReaction<Traits>::get_muparser_initial(
  const ParameterTree& compartment_config,
  const GFGridView& gf_grid_view,
  bool compile)
{
  return get_muparser_expressions(
    compartment_config.sub("initial",true), gf_grid_view, compile);
}

template<class Traits>
template<class GF>
void
ModelDiffusionReaction<Traits>::set_initial(
  const std::vector<std::shared_ptr<GF>>& initial)
{
  using GridFunction = std::decay_t<GF>;
  static_assert(Concept::isPDELabGridFunction<GridFunction>(),
                "GridFunction is not a PDElab grid functions");
  static_assert(
    std::is_same_v<typename GridFunction::Traits::GridViewType, GV>,
    "GridFunction has to have the same grid view as the templated grid view");
  static_assert(
    std::is_same_v<typename GridFunction::Traits::GridViewType,
                   typename Grid::Traits::LeafGridView>,
    "GridFunction has to have the same grid view as the templated grid");
  static_assert((int)GridFunction::Traits::dimDomain == (int)Grid::dimension,
                "GridFunction has to have domain dimension equal to the grid");
  static_assert(GridFunction::Traits::dimRange == 1,
                "GridFunction has to have range dimension equal to 1");

  _logger.debug("Set initial state from grid functions"_fmt);

  if (initial.size() !=
      _config.sub(_compartment_name + ".diffusion",true).getValueKeys().size())
    DUNE_THROW(RangeError, "Wrong number of grid functions");

  PDELab::DynamicPowerGridFunction<GridFunction> comp_initial(initial);
  Dune::PDELab::interpolate(
    comp_initial, *_state.grid_function_space, *_state.coefficients);
}

template<class Traits>
auto
ModelDiffusionReaction<Traits>::setup_component_grid_function_space(
  const std::string& name) const
{
  _logger.trace("Create a finite element map"_fmt);

  // create entity mapper context
  using Entity = typename Traits::Grid::LeafGridView::template Codim<0>::Entity;
  using Index = std::size_t;
  using EntityMapper = std::function<Index(Entity)>;
  EntityMapper&& em = [](const auto entity) {
    return entity.geometry().type().isCube() ? 0 : 1;
  };

  auto create_fem = [&](GeometryType gt){
    // create data context with entity mapper, geometry type and grid view
    auto&& ctx = Context::data_context(em, gt, _grid_view);

    // create fem from factory
    return Factory<FEM>::create(ctx);
  };

  // if the grid view has no entities, we cannot know which geometry type the
  // local finite element is execpting in the constructor. So we just try with
  // simplicies and cubes and expect that it accepts one. Since the grid view
  // has no entities, it doesn't really matter what kind of mapping we create...
  std::shared_ptr<FEM> finite_element_map;

  constexpr auto partition = PartitionIteratorType::Interior_Partition;
  auto begin = _grid_view.template begin<0,partition>();
  auto end = _grid_view.template end<0,partition>();

  if (begin == end) {
    try {
      finite_element_map = create_fem(GeometryTypes::simplex(Grid::dimension));
    } catch (...) {
      finite_element_map = create_fem(GeometryTypes::cube(Grid::dimension));
    }
  } else {
    // get common geometry type on gridview
    if (not has_single_geometry_type(_grid_view))
      DUNE_THROW(InvalidStateException,
                "Grid view has to have only one geometry type");

    // in this case we know the geometry type for the finite element
    GeometryType gt = begin->geometry().type();
    finite_element_map = create_fem(gt);
  }

  _logger.trace("Setup grid function space for component {}"_fmt, name);
  const ES entity_set(_grid->leafGridView());
  auto comp_gfs = std::make_shared<LGFS>(entity_set, finite_element_map);
  comp_gfs->name(name);

  if (begin != end) {
    auto entity_converter = [&](const auto& e){
      if constexpr (Concept::isMultiDomainGrid<typename Traits::Grid>())
        return _grid->multiDomainEntity(e);
      else
        return e;
    };
    std::size_t order =
      finite_element_map->find(entity_converter(*begin)).localBasis().order();

    if (order == 0 and _grid->comm().size() > 1)
      DUNE_THROW(NotImplemented,
                 "Parallel programs on finite volume spaces are not supported");

    comp_gfs->setDataSetType(
      order == 0 ? PDELab::GridFunctionOutputParameters::Output::cellData
                : PDELab::GridFunctionOutputParameters::Output::vertexData);
  }
  return comp_gfs;
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_grid_function_space()
{
  _logger.debug("Setup domain grid function space"_fmt);
  auto components = _config.sub(_compartment_name + ".reaction",true).getValueKeys();

  if (not _state) {
    _state.grid = _grid;
    _state.time = _config.get("time_stepping.begin",0.);
  }
  auto comp_gfs_vec = typename GFS::NodeStorage{};
  for (const auto& name : components)
    comp_gfs_vec.push_back(setup_component_grid_function_space(name));

  _logger.trace("Setup domian power grid function space"_fmt);
  _logger.info("No. of components {}"_fmt, comp_gfs_vec.size());
  _state.grid_function_space = std::make_shared<GFS>(comp_gfs_vec);
  _state.grid_function_space->name(_compartment_name);

  if (_state.grid_function_space->degree() < 1)
    DUNE_THROW(InvalidStateException,
               "Grid function space is not correctly setup");
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_coefficient_vector()
{
  _logger.debug("Setup coefficient vector"_fmt);
  const auto& gfs = _state.grid_function_space;
  if (not gfs)
    DUNE_THROW(InvalidStateException,"Grid function space is not setup");
  _state.coefficients = std::make_shared<X>(*gfs); // create new storage
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_initial_condition()
{
  // if this is a sub model, instantiation of the following will fail
  if constexpr (not Traits::is_sub_model) {
    // get TIFF data if available
    MuParserDataHandler<TIFFGrayscale<unsigned short>> mu_data_handler;
    if (_config.hasSub("data"))
      mu_data_handler.add_tiff_functions(_config.sub("data"));

    // configure math parsers for initial conditions on each component
    auto initial_muparser =
      get_muparser_initial(_config.sub(_compartment_name), _grid_view, false);

    for (auto&& mu_grid_function : initial_muparser) {
      // make TIFF expression available in the parser
      mu_data_handler.set_functions(mu_grid_function->parser());
      mu_grid_function->set_time(_state.time);
      mu_grid_function->compile_parser();
    }
    // evaluate compiled expressions as initial conditions
    set_initial(initial_muparser);
  } else {
    DUNE_THROW(NotImplemented,
               "Initial condition setup is only available for complete models");
  }
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_constraints()
{
  _logger.debug("Setup constraints"_fmt);

  auto b0lambda = [&](const auto& i, const auto& x) { return false; };
  auto b0 =
    Dune::PDELab::makeBoundaryConditionFromCallable(_grid_view, b0lambda);

  _logger.trace("Assemble constraints"_fmt);
  _constraints = std::make_unique<CC>();
  auto& gfs = _state.grid_function_space;
  Dune::PDELab::constraints(b0, *gfs, *_constraints);

  _logger.detail(
    "Constrained dofs: {} of {}"_fmt, _constraints->size(), gfs->globalSize());
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_local_operator()
{
  _logger.debug("Setup local operator"_fmt);

  _logger.trace("Create spatial local operator"_fmt);
  _local_operator =
    std::make_shared<LOP>(_grid_view, _config.sub(_compartment_name));

  _logger.trace("Create temporal local operator"_fmt);
  _temporal_local_operator =
    std::make_shared<TLOP>(_grid_view, _config.sub(_compartment_name));
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_grid_operator()
{
  _logger.debug("Create grid operator"_fmt);

  auto& gfs = _state.grid_function_space;
  auto& lop = _local_operator;
  auto& tlop = _temporal_local_operator;

  // @todo fix this estimate for something more accurate
  MBE mbe((int)Dune::power((int)3, (int)Grid::dimension));

  _logger.trace("Create spatial grid operator"_fmt);
  _spatial_grid_operator =
    std::make_shared<GOS>(*gfs, *_constraints, *gfs, *_constraints, *lop, mbe);

  _logger.trace("Create temporal grid operator"_fmt);
  _temporal_grid_operator =
    std::make_shared<GOT>(*gfs, *_constraints, *gfs, *_constraints, *tlop, mbe);

  _logger.trace("Create instationary grid operator"_fmt);
  _grid_operator =
    std::make_shared<GOI>(*_spatial_grid_operator, *_temporal_grid_operator);
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup_vtk_writer()
{
  _logger.debug("Setup VTK writer"_fmt);
  // this map contains the written timesteps for a given path and is shared
  // among different `state`s copied from this one.
  auto mapu_out =
    std::make_shared<std::map<std::string, std::vector<double>>>();
  _state.writer = [=](
                    const auto& state, const auto& path, bool append) mutable {
    auto log_writer = Logging::Logging::componentLogger({}, "writer");
    auto& x = state.coefficients;
    auto& gfs = state.grid_function_space;

    // create directory if necessary
    auto path_entry = fs::directory_entry{ path };
    if (not path_entry.exists()) {
      log_writer.info("Creating output directory '{}'"_fmt,
                      path_entry.path().string());
      std::error_code ec{ 0, std::generic_category() };
      fs::create_directories(path_entry.path(), ec);
      if (ec)
        DUNE_THROW(IOError,
                   "\n Category: " << ec.category().name() << '\n'
                                   << "Value: " << ec.value() << '\n'
                                   << "Message: " << ec.message() << '\n');
    }

    // Recover old timestesps in case something was written before
    auto& timesteps = (*mapu_out)[path.string()];
    std::string name = fmt::format("{}-{}", path.filename().string(), gfs->name());
    if (not append) {
      timesteps.clear();
      log_writer.detail("Creating a time sequence file: '{}.pvd'"_fmt, name);
    } else {
      log_writer.trace("Overriding time sequence file: '{}.pvd'"_fmt, name);
    }

    // setup writer again with old timesteps if necessary
    using VTKGridView = typename GFS::Traits::GridView;
    auto writer =
      std::make_shared<VTKWriter<VTKGridView>>(gfs->gridView(), Dune::VTK::conforming);
    auto sequential_writer = VTKSequenceWriter{ writer, name, path.string(), "" };
    sequential_writer.setTimeSteps(timesteps);

    using Predicate = PDELab::vtk::DefaultPredicate;
    using Data = PDELab::vtk::DGFTreeCommonData<GFS, X, Predicate, VTKGridView>;
    auto data = std::make_shared<Data>(*gfs, *x, gfs->gridView());
    auto collector = PDELab::vtk::OutputCollector{sequential_writer, data};
    for (std::size_t k = 0; k < data->_lfs.degree(); k++)
      collector.addSolution(data->_lfs.child(k),
                            PDELab::vtk::defaultNameScheme());

    log_writer.detail("Writing solution for {:.2f}s time stamp"_fmt, state.time);
    log_writer.trace("Writing vtu file: '{0}/{0}-{1:0>5}.vtu'"_fmt,
                     name,
                     timesteps.size());
    sequential_writer.write(state.time, Dune::VTK::base64);
    sequential_writer.vtkWriter()->clear();
    timesteps = sequential_writer.getTimeSteps();
  };
}

template<class Traits>
void
ModelDiffusionReaction<Traits>::setup(BitFlags<ModelSetup::Stages> setup_policy)
{
  auto msg = "Setting up diffusion-reaction model for {} compartment"_fmt;
  if (Traits::is_sub_model)
    _logger.info(msg,_compartment_name);
  else
    _logger.notice(msg,_compartment_name);

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
    if (setup_policy.test(ModelSetup::Stages::Writer))
      setup_vtk_writer();
  }
  catch(...)
  {
    // detailed report of the input paramter tree
    std::stringstream ss;
    _config.report(ss);
    _logger.error(2, "---- Diffusion-Reaction model parameter tree"_fmt);
    for (std::string line; std::getline(ss, line);)
      _logger.error(2, "{}"_fmt, line);
    _logger.error(2, "----"_fmt);
    throw;
  }
}

template<class Traits>
auto
ModelDiffusionReaction<Traits>::get_data_handler(const ConstState& state)
{
  std::shared_ptr<DataHandler> data;
  data = std::make_shared<DataHandler>(*state.grid_function_space,
                                       *state.coefficients,
                                       state.grid_function_space->gridView());
  return data;
}

template<class Traits>
auto
ModelDiffusionReaction<Traits>::get_grid_function(const ConstState& state,
                                                  std::size_t comp)
  -> std::shared_ptr<ComponentGridFunction>
{
  auto data = get_data_handler(state);
  return std::make_shared<ComponentGridFunction>(data->_lfs.child(comp), data);
}

template<class Traits>
auto
ModelDiffusionReaction<Traits>::get_grid_function(std::size_t comp) const
  -> std::shared_ptr<ComponentGridFunction>
{
  return get_grid_function(const_state(), comp);
}

template<class Traits>
auto
ModelDiffusionReaction<Traits>::get_grid_functions(const ConstState& state)
  -> std::vector<std::shared_ptr<ComponentGridFunction>>
{
  std::vector<std::shared_ptr<ComponentGridFunction>> grid_functions;
  for (std::size_t i = 0; i < state.grid_function_space->degree(); i++)
    grid_functions.push_back(get_grid_function(state, i));
  return grid_functions;
}

template<class Traits>
auto
ModelDiffusionReaction<Traits>::get_grid_functions() const
  -> std::vector<std::shared_ptr<ComponentGridFunction>>
{
  return get_grid_functions(const_state());
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_CC
