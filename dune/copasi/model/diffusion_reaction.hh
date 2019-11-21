#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH

#include <dune/copasi/common/coefficient_mapper.hh>
#include <dune/copasi/common/enum.hh>
#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/finite_element_map/dynamic_power.hh>
#include <dune/copasi/finite_element_map/multidomain.hh>
#include <dune/copasi/model/base.hh>
#include <dune/copasi/model/local_operator_base.hh>
#include <dune/copasi/model/local_operator_FV.hh>
#include <dune/copasi/model/local_operator_CG.hh>
#include <dune/copasi/model/state.hh>
#include <dune/copasi/finite_element/p0.hh>
#include <dune/copasi/finite_element/pk.hh>
#include <dune/copasi/finite_element_map/p0.hh>
#include <dune/copasi/finite_element_map/pk.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/p0.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/function/discretegridviewfunction.hh>
#include <dune/pdelab/gridfunctionspace/dynamicpowergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <memory>

namespace Dune::Copasi {

/**
 * @brief      Traits for diffusion reaction models with Pk elements
 *
 * @tparam     G         Grid
 * @tparam     GV        Grid View
 * @tparam     FEMorder  Order of the finite element method
 * @tparam     OT        PDELab ordering tag
 * @tparam     JM        Jacobian method
 */
template<class G,
         class GV = typename G::Traits::LeafGridView,
         class OT = PDELab::EntityBlockedOrderingTag,
         JacobianMethod JM = JacobianMethod::Analytical>
struct ModelP0DiffusionReactionTraits
{
  using Grid = G;
  using GridView = GV;
  using BaseFEM =
    PDELab::P0LocalFiniteElementMap<typename Grid::ctype,
                                    double,2>;

  static constexpr bool is_sub_model = not std::is_same_v<typename Grid::Traits::LeafGridView,GridView>;

  //! Finite element map
  using FEM = std::conditional_t<
                  is_sub_model,
                  MultiDomainLocalFiniteElementMap<BaseFEM,GridView>,
                  BaseFEM
                >;

  using OrderingTag = OT;
  static constexpr JacobianMethod jacobian_method = JM;

  //! Local operator
  template <class CoefficientMapper>
  using LocalOperator = LocalOperatorDiffusionReactionFV<
    GridView,
    typename BaseFEM::Traits::FiniteElement::Traits::LocalBasisType::Traits,
    CoefficientMapper,
    jacobian_method>;

  //! Temporal local operator
  template <class CoefficientMapper>
  using TemporalLocalOperator = TemporalLocalOperatorDiffusionReactionCG<
    GridView,
    typename BaseFEM::Traits::FiniteElement::Traits::LocalBasisType::Traits,
    jacobian_method>;
};

template<class G,
         class GV = typename G::Traits::LeafGridView,
         int FEMorder = 1,
         class OT = PDELab::EntityBlockedOrderingTag,
         JacobianMethod JM = JacobianMethod::Analytical>
struct ModelPkDiffusionReactionTraits
{
  using Grid = G;
  using GridView = GV;
  using BaseFEM =
    PDELab::PkLocalFiniteElementMap<typename G::LeafGridView, double,double,FEMorder>;

  static constexpr bool is_sub_model = not std::is_same_v<typename Grid::Traits::LeafGridView,GridView>;

  //! Finite element map
  using FEM = std::conditional_t<
                  is_sub_model,
                  MultiDomainLocalFiniteElementMap<BaseFEM,GridView>,
                  BaseFEM
                >;

  using OrderingTag = OT;
  static constexpr JacobianMethod jacobian_method = JM;

  //! Local operator
  template <class CoefficientMapper>
  using LocalOperator = LocalOperatorDiffusionReactionCG<
    GridView,
    typename BaseFEM::Traits::FiniteElement::Traits::LocalBasisType::Traits,
    CoefficientMapper,
    jacobian_method>;

  //! Temporal local operator
  template <class CoefficientMapper>
  using TemporalLocalOperator = TemporalLocalOperatorDiffusionReactionCG<
    GridView,
    typename BaseFEM::Traits::FiniteElement::Traits::LocalBasisType::Traits,
    jacobian_method>;
};

template<class G, class GV, class OT, JacobianMethod JM>
struct ModelPkDiffusionReactionTraits<G,GV,0,OT,JM> : public ModelP0DiffusionReactionTraits<G,GV,OT,JM> {};

/**
 * @brief      Class for diffusion-reaction models.
 *
 * @tparam     Traits  Class that define static policies on the model
 */
template<class Traits>
class ModelDiffusionReaction : public ModelBase
{
public:
  using Grid = typename Traits::Grid;

  using OT = typename Traits::OrderingTag;

  static constexpr JacobianMethod JM = Traits::jacobian_method;

  // Check templates
  static_assert(Concept::isGrid<Grid>(), "Provided an invalid grid");

  //! Grid view
  using GV = typename Traits::GridView;

  //! Host grid view
  using HGV = typename Traits::Grid::LeafGridView;

  //! Domain field
  using DF = typename Grid::ctype;

  //! Finite element map
  using FEM = typename Traits::FEM;

  //! Range field
  using RF = typename FEM::Traits::FiniteElement::Traits::
    LocalBasisType::Traits::RangeFieldType;

  //! Constraints builder
  using CON = PDELab::P0ParallelConstraints;

  //! Entity set
  using ES = Dune::PDELab::OverlappingEntitySet<HGV>;

  //! Leaf vector backend
  using LVBE = PDELab::ISTL::VectorBackend<>;

  //! Leaf grid function space
  using LGFS = PDELab::GridFunctionSpace<ES, FEM, CON, LVBE>;

  //! Vector backend
  using VBE = LVBE;

  //! Grid function space
  using GFS = PDELab::DynamicPowerGridFunctionSpace<LGFS, VBE, OT>;

  //! Coefficient vector
  using X = PDELab::Backend::Vector<GFS, RF>;

public:
  //! Model state structure
  using State = Dune::Copasi::ModelState<Grid, GFS, X>;

  //! Constant model state structure
  using ConstState = Dune::Copasi::ConstModelState<Grid, GFS, X>;

private:
  //! Constraints container
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;

  //! Coefficients mapper
  using CM = Dune::Copasi::ModelCoefficientMapper<ConstState>;

  //! Local operator
  using LOP = typename Traits::template LocalOperator<CM>;

  //! Temporal local operator
  using TLOP = typename Traits::template TemporalLocalOperator<CM>;

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
  using W = Dune::VTKWriter<GV>;

  //! Sequential writer
  using SW = Dune::VTKSequenceWriter<GV>;

  using DataHandler =
    PDELab::vtk::DGFTreeCommonData<GFS,
                                   X,
                                   PDELab::vtk::DefaultPredicate,
                                   GV>;
  using ComponentLFS =
    typename PDELab::LocalFunctionSpace<GFS>::ChildType;

  using ComponentGridFunction = PDELab::vtk::
    DGFTreeLeafFunction<ComponentLFS, DataHandler, GV>;

public:
  /**
   * @brief      Constructs the model
   *
   * @param[in]  grid          The grid
   * @param[in]  config        The configuration file
   * @param[in]  grid_view     The grid view to operate with
   * @param[in]  setup_policy  Policy to setup model
   */
  ModelDiffusionReaction(std::shared_ptr<Grid> grid,
                         const Dune::ParameterTree& config,
                         GV grid_view,
                         ModelSetupPolicy setup_policy = ModelSetupPolicy::All);

  /**
   * @brief      Constructs the model
   * @details    This constructor only is available if the grid view
   *             is the leaf grid view of the templated grid
   *
   * @param[in]  grid          The grid
   * @param[in]  config        The configuration file
   * @param[in]  setup_policy  Policy to setup model
   */
  template<class T = int,
           class = std::enable_if_t<
             std::is_same_v<GV, typename Grid::Traits::LeafGridView>,
             T>>
  ModelDiffusionReaction(std::shared_ptr<Grid> grid,
                         const Dune::ParameterTree& config,
                         ModelSetupPolicy setup_policy = ModelSetupPolicy::All)
    : ModelDiffusionReaction(grid, config, grid->leafGridView(), setup_policy)
  {}

  /**
   * @brief      Destroys the model
   */
  ~ModelDiffusionReaction();

  /**
   * @brief      Suggest a time step to the model.
   *
   * @param[in]  dt    Suggestion for the internal time step of the model.
   */
  void suggest_timestep(double dt);

  /**
   * @brief      Performs one steps in direction to end_time().
   * @details    The time-step should never result on a bigger step than the one
   *             suggested in suggest_timestep().
   */
  void step();

  /**
   * @brief      Get mutable model states
   *
   * @return     Model states
   */
  std::map<std::size_t, State> states()
  {
    for (auto& [op, state] : _states)
      state.time = current_time();
    return _states;
  }

  /**
   * @brief      Get constat model states
   *
   * @return     Constant model states
   */
  std::map<std::size_t, ConstState> const_states() const
  {
    std::map<std::size_t, ConstState> const_states(_states.begin(),
                                                   _states.end());
    return const_states;
  }

  /**
   * @brief      Get constat model states
   *
   * @return     Constant model states
   */
  std::map<std::size_t, ConstState> states() const { return const_states(); }

  /**
   * @brief      Sets the initial state of the model
   *
   * @param[in]  model_config  A parameter tree with 'initial' and optionally
   * 'data' subsections
   */
  template<class GFGridView>
  static auto get_muparser_initial(const ParameterTree& model_config,
                                   const GFGridView& gf_grid_view, bool compile = true);

  /**
   * @brief      Sets the initial state of the model
   * @details    The input vector should have the same size as the number of
   * variables in the model. Additionally, they will be indepreted aphabetically
   * accodingly to the name set to othe input sections (e.g. 'model.diffusion'
   * section).
   *
   * @tparam     GF       A valid PDELab grid functions (see
   * @Concepts::PDELabGridFunction)
   * @param[in]  initial  Vector of grid functions for each variable
   */
  template<class GF>
  void set_initial(const std::vector<std::shared_ptr<GF>>& initial);

  auto get_data_handler(std::map<std::size_t, State> states) const;

  void update_data_handler();

  auto get_grid_function(const std::map<std::size_t, State>& states,
                         std::size_t comp) const;

  auto get_grid_function(std::size_t comp) const;

  auto get_grid_functions(const std::map<std::size_t, State>& states) const;

  auto get_grid_functions() const;

protected:
  auto setup_component_grid_function_space(std::string) const;
  auto setup_domain_grid_function_space(std::vector<std::string>) const;
  void setup_grid_function_space();
  void setup_coefficient_vectors();
  void setup_constraints();
  auto setup_local_operator(std::size_t) const;
  void setup_local_operators();
  void setup_grid_operators();
  void setup_solvers();
  void setup_vtk_writer();
  void write_states() const;

  /**
   * @brief      Setup for next time step
   */
  void setup(ModelSetupPolicy setup_policy = ModelSetupPolicy::All);

  using ModelBase::_logger;
  Logging::Logger _solver_logger;

private:
  std::size_t _components;
  ParameterTree _config;
  GV _grid_view;
  std::map<std::size_t, State> _states;
  std::multimap<std::size_t, std::string> _operator_splitting;

  std::shared_ptr<Grid> _grid;

  std::map<std::size_t, std::shared_ptr<DataHandler>> _data;

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

  std::shared_ptr<W> _writer;
  std::shared_ptr<SW> _sequential_writer;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH