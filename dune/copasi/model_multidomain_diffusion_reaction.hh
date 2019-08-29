#ifndef DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH

#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/local_operator_multidomain.hh>
#include <dune/copasi/model_base.hh>
#include <dune/copasi/model_diffusion_reaction.cc>
#include <dune/copasi/model_diffusion_reaction.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/backend/istl/novlpistlsolverbackend.hh>
#include <dune/pdelab/gridfunctionspace/dynamicpowergridfunctionspace.hh>

#include <dune/common/parametertree.hh>

#include <array>
#include <memory>

namespace Dune::Copasi {

template<class Grid,
         int FEMorder = 1,
         class OrderingTag = PDELab::EntityBlockedOrderingTag>
class ModelMultiDomainDiffusionReaction : public ModelBase
{
  // Check grid template
  static_assert(Concept::isMultiDomainGrid<Grid>(),
                "Provided grid type is not a multidomain grid");

  static_assert(Concept::isSubDomainGrid<typename Grid::SubDomainGrid>());

  //! World dimension
  static constexpr int dim = 2;
  //! Polynomial order
  static constexpr int order = FEMorder;

  using SubDomainGridView = typename Grid::SubDomainGrid::LeafGridView;
  using SubModel =
    ModelDiffusionReaction<Grid, SubDomainGridView, order, OrderingTag>;

  //! Grid view
  using GridView = typename Grid::LeafGridView;
  using GV = GridView;

  //! Domain field
  using DF = typename Grid::ctype;

  //! Range field
  using RF = double;

  //! Finite element
  using FE = Dune::PkLocalFiniteElement<DF, RF, dim, order>;

  //! Base finite element map
  using BaseFEM = PDELab::PkLocalFiniteElementMap<GV, DF, RF, order>;

  //! Finite element map
  using FEM = MultiDomainLocalFiniteElementMap<BaseFEM, SubDomainGridView>;

  //! Constraints builder
  using CON = PDELab::ConformingDirichletConstraints;

  //! Leaf vector backend
  using LVBE = PDELab::ISTL::VectorBackend<>;

  //! Leaf grid function space
  using LGFS = typename SubModel::LGFS;

  //! SubDomain grid function space
  using SDGFS = typename SubModel::GFS;

  using VBE = typename LGFS::Traits::Backend;
  using GFS = PDELab::DynamicPowerGridFunctionSpace<SDGFS, VBE>;

  //! Coefficient vector
  using X = PDELab::Backend::Vector<GFS, RF>;

  //! Constraints container
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;

  //! Local operator
  using LOP = LocalOperatorMultiDomainDiffusionReaction<Grid, FE>;

  //! Temporal local operator
  using TLOP = TemporalLocalOperatorMultiDomainDiffusionReaction<Grid, FE>;

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
  using LS = Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GOI>;

  //! Nonlinear solver
  using NLS = Dune::PDELab::Newton<GOI, LS, X>;

  //! Time stepping parameter
  using TSP = Dune::PDELab::TimeSteppingParameterInterface<double>;

  //! One step method
  using OSM = Dune::PDELab::OneStepMethod<RF, GOI, NLS, X, X>;

  //! Writer
  using W = Dune::VTKWriter<SubDomainGridView>;

  //! Sequential writer
  using SW = Dune::VTKSequenceWriter<SubDomainGridView>;

  //! Model state structure
  using State = Dune::Copasi::ModelState<Grid, GFS, X>;

  //! Constant model state structure
  using ConstState = Dune::Copasi::ConstModelState<Grid, GFS, X>;

public:
  ModelMultiDomainDiffusionReaction(std::shared_ptr<Grid> grid,
                                    const Dune::ParameterTree& config);

  ~ModelMultiDomainDiffusionReaction();

  void setup(ModelSetupPolicy setup_policy = ModelSetupPolicy::All);

  void suggest_timestep(double dt) override;

  void step() override;

protected:
  void setup_grid_function_spaces();
  // void setup_coefficient_vector();
  // void setup_constraints();
  // void setup_local_operators();
  // void setup_grid_operators();
  // void setup_solvers();
  // void setup_vtk_writer();

private:
  using ModelBase::_logger;
  Logging::Logger _solver_logger;

  ParameterTree _config;
  GV _grid_view;

  std::vector<std::shared_ptr<SW>> _sequential_writer;

  std::vector<State> _states;
  std::shared_ptr<Grid> _grid;

  std::unique_ptr<CC> _constraints;
  std::vector<std::shared_ptr<LOP>> _local_operators;
  std::vector<std::shared_ptr<TLOP>> _temporal_local_operators;
  std::vector<std::shared_ptr<GOS>> _spatial_grid_operators;
  std::vector<std::shared_ptr<GOT>> _temporal_grid_operators;
  std::vector<std::shared_ptr<GOI>> _grid_operators;
  std::vector<std::shared_ptr<LS>> _linear_solvers;
  std::vector<std::shared_ptr<NLS>> _nonlinear_solvers;
  std::vector<std::shared_ptr<TSP>> _time_stepping_methods;
  std::vector<std::shared_ptr<OSM>> _one_step_methods;

  std::size_t _domains;
  std::size_t _dof_per_component;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH