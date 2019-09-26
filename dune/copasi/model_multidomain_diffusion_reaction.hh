#ifndef DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH

#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/enum.hh>
#include <dune/copasi/local_operator_multidomain.hh>
#include <dune/copasi/model_base.hh>
#include <dune/copasi/model_diffusion_reaction.cc>
#include <dune/copasi/model_diffusion_reaction.hh>
#include <dune/copasi/multidomain_entity_transformation.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/backend/istl/novlpistlsolverbackend.hh>
#include <dune/pdelab/gridfunctionspace/dynamicpowergridfunctionspace.hh>

#include <dune/common/parametertree.hh>

#include <array>
#include <memory>

namespace Dune::Copasi {

template<class G,
         int FEMorder = 1,
         class OT = PDELab::EntityBlockedOrderingTag,
         JacobianMethod JM = JacobianMethod::Analytical>
struct ModelMultiDomainDiffusionReactionTraits
{
  using Grid = G;
  using OrderingTag = OT;
  static constexpr int order = FEMorder;
  static constexpr JacobianMethod jacobian_method = JM;
};

template<class Traits>
class ModelMultiDomainDiffusionReaction : public ModelBase
{
  using Grid = typename Traits::Grid;

  using OT = typename Traits::OrderingTag;

  static constexpr JacobianMethod JM = Traits::jacobian_method;

  // Check grid template
  static_assert(Concept::isMultiDomainGrid<Grid>(),
                "Provided grid type is not a multidomain grid");

  static_assert(Concept::isSubDomainGrid<typename Grid::SubDomainGrid>());

  //! World dimension
  static constexpr int dim = 2;
  //! Polynomial order
  static constexpr int order = Traits::order;

  using SubDomainGridView = typename Grid::SubDomainGrid::LeafGridView;

  using SubModelTraits =
    ModelDiffusionReactionTraits<Grid, SubDomainGridView, order, OT>;
  using SubModel = ModelDiffusionReaction<SubModelTraits>;

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

  //! Model state structure
  using State = Dune::Copasi::ModelState<Grid, GFS, X>;

  //! Constant model state structure
  using ConstState = Dune::Copasi::ConstModelState<Grid, GFS, X>;

  //! Coefficients mapper
  using CM = Dune::Copasi::MultiDomainModelCoefficientMapper<ConstState>;

  //! Local operator
  using LOP = LocalOperatorMultiDomainDiffusionReaction<Grid, FE, CM, JM>;

  //! Temporal local operator
  using TLOP = TemporalLocalOperatorMultiDomainDiffusionReaction<Grid, FE, JM>;

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

  //! Entity transformation between grids
  using EntityTransformation =
    Dune::Copasi::MultiDomainEntityTransformation<Grid>;

  using DataHandler =
    PDELab::vtk::DGFTreeCommonData<GFS,
                                   X,
                                   PDELab::vtk::DefaultPredicate,
                                   SubDomainGridView,
                                   EntityTransformation>;

  using ComponentLFS =
    typename PDELab::LocalFunctionSpace<GFS>::ChildType::ChildType;

  using ComponentGridFunction = PDELab::vtk::
    DGFTreeLeafFunction<ComponentLFS, DataHandler, SubDomainGridView>;

public:
  ModelMultiDomainDiffusionReaction(std::shared_ptr<Grid> grid,
                                    const Dune::ParameterTree& config);

  ~ModelMultiDomainDiffusionReaction();

  /**
   * @brief      Get the model state
   *
   * @return     Model state
   */
  std::map<std::size_t, State> states() { return _states; }

  /**
   * @brief      Get the model state
   *
   * @return     Model state
   */
  std::map<std::size_t, ConstState> const_states() const
  {
    return const_states(_states);
  }

  std::map<std::size_t, ConstState> const_states(
    const std::map<std::size_t, State>& states) const
  {
    std::map<std::size_t, ConstState> const_states(states.begin(),
                                                   states.end());
    return const_states;
  }

  /**
   * @brief      Get the model state
   *
   * @return     Model state
   */
  std::map<std::size_t, ConstState> states() const { return const_states(); }

  void setup(ModelSetupPolicy setup_policy = ModelSetupPolicy::All);

  void suggest_timestep(double dt) override;

  void step() override;

  auto get_grid_function(const std::map<std::size_t, State>&,
                         std::size_t,
                         std::size_t) const;

protected:
  void setup_grid_function_spaces();
  void setup_coefficient_vectors();
  void setup_constraints();
  auto setup_local_operator(std::size_t) const;
  void setup_local_operators();
  void setup_grid_operators();
  void setup_solvers();
  void setup_vtk_writer();
  void write_states() const;

  void update_data_handler();
  auto get_data_handler(std::map<std::size_t, State>) const;

private:
  using ModelBase::_logger;
  Logging::Logger _solver_logger;

  ParameterTree _config;
  GV _grid_view;

  std::vector<std::shared_ptr<SW>> _sequential_writer;

  std::map<std::size_t, State> _states;
  std::shared_ptr<Grid> _grid;

  std::vector<std::map<std::size_t, std::shared_ptr<DataHandler>>> _data;

  std::map<std::size_t, std::unique_ptr<CC>> _constraints;
  std::map<std::size_t, std::shared_ptr<LOP>> _local_operators;
  std::map<std::size_t, std::shared_ptr<TLOP>> _temporal_local_operators;
  std::map<std::size_t, std::shared_ptr<GOS>> _spatial_grid_operators;
  std::map<std::size_t, std::shared_ptr<GOT>> _temporal_grid_operators;
  std::map<std::size_t, std::shared_ptr<GOI>> _grid_operators;
  std::map<std::size_t, std::shared_ptr<LS>> _linear_solvers;
  std::map<std::size_t, std::shared_ptr<NLS>> _nonlinear_solvers;
  std::map<std::size_t, std::shared_ptr<TSP>> _time_stepping_methods;
  std::map<std::size_t, std::shared_ptr<OSM>> _one_step_methods;

  std::size_t _domains;
  std::size_t _dof_per_component;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH