#ifndef DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH

#ifndef DUNE_COPASI_MAX_SUBDOMAINS
#define DUNE_COPASI_MAX_SUBDOMAINS 3
#endif

#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/local_operator_multidomain.hh>
#include <dune/copasi/model_base.hh>
#include <dune/copasi/model_diffusion_reaction.cc>
#include <dune/copasi/model_diffusion_reaction.hh>

#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>

#include <dune/istl/io.hh>

#include <dune/common/parametertree.hh>

#include <array>
#include <memory>

namespace Dune::Copasi {

template<class Grid>
class ModelMultiDomainDiffusionReaction : public ModelBase
{
  // Check template
  static_assert(Concept::isMultiDomainGrid<Grid>(),
                "Provided grid type is not a multidomain grid");

  // just check concept for subdomain is working properly
  static_assert(Concept::isSubDomainGrid<typename Grid::SubDomainGrid>());

  static constexpr std::size_t max_subdomains = DUNE_COPASI_MAX_SUBDOMAINS;

  using SubDomainGridView = typename Grid::SubDomainGrid::LeafGridView;
  using GridView = typename Grid::LeafGridView;
  using Model = ModelDiffusionReaction<Grid, SubDomainGridView>;
  using ModelStorage = std::array<std::shared_ptr<Model>, max_subdomains>;

  //! World dimension
  static constexpr int dim = 2;
  //! Polynomial order
  static constexpr int order = 1;

  //! Domain field
  using DF = typename Grid::ctype;

  //! Range field
  using RF = double;

  //! Finite element
  using FE = Dune::PkLocalFiniteElement<DF, RF, dim, order>;

  //! Base finite element map
  using BaseFEM = PDELab::PkLocalFiniteElementMap<GridView, DF, RF, order>;

  //! Finite element map
  using FEM = DynamicPowerLocalFiniteElementMap<BaseFEM, SubDomainGridView>;

  //! Constraints builder
  using CON = PDELab::ConformingDirichletConstraints;

  //! Leaf vector backend
  using LVBE = PDELab::ISTL::VectorBackend<>;

  //! Leaf grid function space
  using LGFS = PDELab::GridFunctionSpace<GridView, FEM, CON, LVBE>;

  // using LGFS = typename Model::State::GridFunctionSpace;

  using VBE = typename LGFS::Traits::Backend;
  using GFS = PDELab::PowerGridFunctionSpace<LGFS, max_subdomains, VBE>;

  //! Coefficient vector
  using X = PDELab::Backend::Vector<GFS, RF>;

  //! Constraints container
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;

  //! Local operator
  using LOP =
    LocalOperatorMultiDomainDiffusionReaction<Grid, FE, max_subdomains>;

  //! Temporal local operator
  using TLOP =
    TemporalLocalOperatorMultiDomainDiffusionReaction<Grid, FE, max_subdomains>;

  //! Matrix backend
  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;

  //! Spatial grid operator
  using GOS =
    Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, RF, RF, RF, CC, CC>;

  //! Temporal grid operator
  using GOT =
    Dune::PDELab::GridOperator<GFS, GFS, TLOP, MBE, RF, RF, RF, CC, CC>;

  //! Instationary grid operator
  using GOI = Dune::PDELab::OneStepGridOperator<GOS, GOT>;

  //! Linear solver backend
  using LS = Dune::PDELab::ISTLBackend_SEQ_BCGS_SSOR;

  //! Nonlinear solver
  using NLS = Dune::PDELab::Newton<GOI, LS, X>;

  //! Time stepping parameter
  using TSP = Dune::PDELab::TimeSteppingParameterInterface<double>;

  //! One step method
  using OSM = Dune::PDELab::OneStepMethod<RF, GOI, NLS, X, X>;

  //! Writer
  using W = Dune::VTKWriter<SubDomainGridView>;

  //! Sequential writer
  using SW = Dune::Copasi::GridFunctionVTKSequenceWriter<SubDomainGridView>;

  //! Discrete grid function
  // using DGF = Dune::Copasi::DiscreteGridFunction<typename GFS::template
  // Child<0>::Type, X, SubDomainGridView, RF, -1, DynamicVector<RF>>;

  //! Model state structure
  using State = Dune::Copasi::ModelState<Grid, GFS, X>;

  //! Constant model state structure
  using ConstState = Dune::Copasi::ModelState<const Grid, const GFS, const X>;

public:
  ModelMultiDomainDiffusionReaction(std::shared_ptr<Grid> grid,
                                    const Dune::ParameterTree& config)
    : ModelBase(config)
    , _config(config)
    , _grid(grid)
  {

    // // create empty models
    // for (int i = _size; i < max_subdomains; ++i)
    // {
    //   auto& model_config = _config.sub(compartements[0]);
    //   SubDomainGridView sub_grid_view = grid->subDomain(i).leafGridView();
    //   _models[i] = std::make_shared<Model>(grid,sub_grid_view,model_config);
    // }

    operator_setup();
  }

  ~ModelMultiDomainDiffusionReaction()
  {
    _logger.debug("ModelMultiDomainDiffusionReaction destroyed"_fmt);
  }

  void operator_setup()
  {
    auto _grid_view = _grid->leafGridView();
    BaseFEM base_fem(_grid_view);
    _dof_per_component = base_fem.maxLocalSize();

    const auto& compartements = _config.sub("compartements").getValueKeys();
    _size = std::min(compartements.size(), max_subdomains);

    typename GFS::NodeStorage leaf_gfs;

    for (int i = 0; i < max_subdomains; ++i) {
      int components = 0;
      if (i < _size) {
        const std::string compartement = compartements[i];
        auto& model_config = _config.sub(compartement);
        auto& intial_config = model_config.sub("initial");
        auto names = intial_config.getValueKeys();
        components = names.size();
      }
      std::cout << components << std::endl;
      SubDomainGridView sub_grid_view = _grid->subDomain(i).leafGridView();

      auto finite_element_map =
        std::make_shared<FEM>(sub_grid_view, base_fem, components);

      leaf_gfs[i] = std::make_shared<LGFS>(_grid_view, finite_element_map);
    }

    _gfs = std::make_shared<GFS>(leaf_gfs);

    _logger.trace("create vector backend"_fmt);
    if (_x)
      _x = std::make_shared<X>(*_gfs, *(_x->storage()));
    else
      _x = std::make_shared<X>(*_gfs);

    _logger.trace("read constraints"_fmt);
    auto b0lambda = [&](const auto& i, const auto& x) { return false; };
    auto b0 =
      Dune::PDELab::makeBoundaryConditionFromCallable(_grid_view, b0lambda);

    _constraints = std::make_unique<CC>();
    Dune::PDELab::constraints(b0, *_gfs, *_constraints);

    _logger.info("constrained dofs: {} of {}"_fmt,
                 _constraints->size(),
                 _gfs->globalSize());

    auto entity_it = _grid_view.template begin<0>();
    FE finite_element;
    assert(finite_element.size() > 0);
    _local_operator = std::make_shared<LOP>(_grid, _config, finite_element);

    _temporal_local_operator =
      std::make_shared<TLOP>(_grid, _config, finite_element);

    MBE mbe((int)pow(3, dim));

    _spatial_grid_operator = std::make_shared<GOS>(
      *_gfs, *_constraints, *_gfs, *_constraints, *_local_operator, mbe);

    _temporal_grid_operator = std::make_shared<GOT>(*_gfs,
                                                    *_constraints,
                                                    *_gfs,
                                                    *_constraints,
                                                    *_temporal_local_operator,
                                                    mbe);

    _grid_operator =
      std::make_shared<GOI>(*_spatial_grid_operator, *_temporal_grid_operator);

    _linear_solver = std::make_shared<LS>(5000, false);

    _nonlinear_solver =
      std::make_shared<NLS>(*_grid_operator, *_x, *_linear_solver);
    _nonlinear_solver->setReassembleThreshold(0.0);
    _nonlinear_solver->setVerbosityLevel(2);
    _nonlinear_solver->setReduction(1e-8);
    _nonlinear_solver->setMinLinearReduction(1e-10);
    _nonlinear_solver->setMaxIterations(25);
    _nonlinear_solver->setLineSearchMaxIterations(10);

    using AlexMethod = Dune::PDELab::Alexander2Parameter<double>;
    _time_stepping_method = std::make_shared<AlexMethod>();
    _one_step_method = std::make_shared<OSM>(
      *_time_stepping_method, *_grid_operator, *_nonlinear_solver);
    _one_step_method->setVerbosityLevel(2);

    using GridFunction = ExpressionToGridFunctionAdapter<GridView, RF>;

    std::array<std::shared_ptr<GridFunction>, max_subdomains> functions;

    for (int i = 0; i < _size; ++i) {
      const std::string compartement = compartements[i];
      auto& model_config = _config.sub(compartement);

      int sub_domain_id =
        _config.sub("compartements").template get<int>(compartement);
      SubDomainGridView sub_grid_view =
        _grid->subDomain(sub_domain_id).leafGridView();

      auto& intial_config = _config.sub(compartement + ".initial");
      functions[i] =
        std::make_shared<ExpressionToGridFunctionAdapter<GridView, RF>>(
          _grid_view, intial_config);
    }

    using PowerGridFunction =
      PDELab::PowerGridFunction<GridFunction, max_subdomains>;
    PowerGridFunction initial(functions);
    Dune::PDELab::interpolate(initial, *_gfs, *_x);

    using Data =
      PDELab::vtk::DGFTreeCommonData<GFS, X, PDELab::vtk::DefaultPredicate>;
    std::shared_ptr<Data> data = std::make_shared<Data>(*_gfs, *_x);

    for (int i = 0; i < _size; ++i) {
      std::cout << "GFS size: " << _gfs->child(i).size() << std::endl;
      const std::string compartement = compartements[i];
      auto& model_config = _config.sub(compartement);

      int sub_domain_id =
        _config.sub("compartements").template get<int>(compartement);
      SubDomainGridView sub_grid_view =
        _grid->subDomain(sub_domain_id).leafGridView();

      std::shared_ptr<W> writer =
        std::make_shared<W>(sub_grid_view, Dune::VTK::conforming);
      std::string filename = model_config.get("name", compartement);
      struct stat st;

      if (stat(filename.c_str(), &st) != 0) {
        int stat = 0;
        stat = mkdir(filename.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
        if (stat != 0 && stat != -1)
          std::cout << "Error: Cannot create directory " << filename
                    << std::endl;
      }

      auto& intial_config = model_config.sub("initial");
      auto names = intial_config.getValueKeys();

      _sequential_writer[i] =
        std::make_shared<SW>(writer, filename, filename, "");
      // add data field for all components of the space to the VTK writer

      using LFS = PDELab::LocalFunctionSpace<GFS>;
      LFS lfs(*_gfs);

      for (int j = 0; j < names.size(); ++j) {
        using DGF = Dune::Copasi::DiscreteGridFunction<Data,
                                                       SubDomainGridView,
                                                       RF,
                                                       1,
                                                       FieldVector<RF, 1>>;
        auto dgf = std::make_shared<const DGF>(data,
                                               sub_grid_view,
                                               j * _dof_per_component,
                                               (j + 1) * _dof_per_component);
        _sequential_writer[i]->addVertexData(dgf, names[j]);
      }
      _sequential_writer[i]->write(current_time(), Dune::VTK::appendedraw);
      _sequential_writer[i]->clear();
    }
  }

  void suggest_timestep(double dt) override
  {
    DUNE_THROW(NotImplemented, "not implemented");
  };

  void step() override
  {
    // for (int i = 0; i < _size; ++i)
    //   _models[i]->step();

    double dt = 0.001;

    // do time step
    auto x_new = std::make_shared<X>(*_x);
    _one_step_method->apply(current_time(), dt, *_x, *x_new);

    // accept time step
    _x = x_new;
    current_time() += dt;

    using Data =
      PDELab::vtk::DGFTreeCommonData<GFS, X, PDELab::vtk::DefaultPredicate>;
    std::shared_ptr<Data> data = std::make_shared<Data>(*_gfs, *_x);

    const auto& compartements = _config.sub("compartements").getValueKeys();
    for (int i = 0; i < _size; ++i) {
      const std::string compartement = compartements[i];
      auto& model_config = _config.sub(compartement);
      auto& intial_config = model_config.sub("initial");
      auto names = intial_config.getValueKeys();

      int sub_domain_id =
        _config.sub("compartements").template get<int>(compartement);
      SubDomainGridView sub_grid_view =
        _grid->subDomain(sub_domain_id).leafGridView();

      using LFS = PDELab::LocalFunctionSpace<GFS>;
      LFS lfs(*_gfs);

      for (int j = 0; j < names.size(); ++j) {
        using DGF = Dune::Copasi::DiscreteGridFunction<Data,
                                                       SubDomainGridView,
                                                       RF,
                                                       1,
                                                       FieldVector<RF, 1>>;
        auto dgf = std::make_shared<const DGF>(data,
                                               sub_grid_view,
                                               j * _dof_per_component,
                                               (j + 1) * _dof_per_component);
        _sequential_writer[i]->addVertexData(dgf, names[j]);
      }
      _sequential_writer[i]->write(current_time(), Dune::VTK::appendedraw);
      _sequential_writer[i]->clear();
    }
  };

private:
  ParameterTree _config;
  using ModelBase::_logger;

  std::array<std::shared_ptr<SW>, max_subdomains> _sequential_writer;

  std::shared_ptr<Grid> _grid;
  // ModelStorage _models;

  std::shared_ptr<GFS> _gfs;
  std::unique_ptr<CC> _constraints;
  std::shared_ptr<LOP> _local_operator;
  std::shared_ptr<TLOP> _temporal_local_operator;
  std::shared_ptr<GOS> _spatial_grid_operator;
  std::shared_ptr<GOT> _temporal_grid_operator;
  std::shared_ptr<GOI> _grid_operator;
  std::shared_ptr<LS> _linear_solver;
  std::shared_ptr<NLS> _nonlinear_solver;
  std::shared_ptr<TSP> _time_stepping_method;
  std::shared_ptr<OSM> _one_step_method;

  std::shared_ptr<X> _x;

  std::size_t _size;
  std::size_t _dof_per_component;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH