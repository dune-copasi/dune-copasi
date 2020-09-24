#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/concepts/grid.hh>

#include <dune/copasi/finite_element/p0.hh>
#include <dune/copasi/finite_element/pk.hh>
#include <dune/copasi/finite_element_map/p0.hh>
#include <dune/copasi/finite_element_map/pk.hh>
#include <dune/copasi/finite_element_map/subdomain.hh>
#include <dune/copasi/finite_element_map/variadic.hh>
#include <dune/copasi/finite_element_map/virtual.hh>
#include <dune/copasi/local_operator/diffusion_reaction/base.hh>
#include <dune/copasi/local_operator/diffusion_reaction/continuous_galerkin.hh>
#include <dune/copasi/local_operator/diffusion_reaction/finite_volume.hh>
#include <dune/copasi/local_operator/diffusion_reaction/multidomain.hh>
#include <dune/copasi/local_operator/variadic.hh>
#include <dune/copasi/model/base.hh>
#include <dune/copasi/model/state.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/p0.hh>
#include <dune/pdelab/function/discretegridviewfunction.hh>
#include <dune/pdelab/gridfunctionspace/dynamicpowergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>

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
  using FEMP0 =
    PDELab::P0LocalFiniteElementMap<typename Grid::ctype, double, G::dimension>;

  static constexpr bool is_sub_model =
    not std::is_same_v<typename Grid::Traits::LeafGridView, GridView>;

  //! Finite element map
  using FEM =
    std::conditional_t<is_sub_model,
                       SubDomainLocalFiniteElementMap<FEMP0, GridView>,
                       FEMP0>;

  using OrderingTag = OT;
  static constexpr JacobianMethod jacobian_method = JM;

  //! Local operator
  using LocalOperator = LocalOperatorDiffusionReactionFV<
    GridView,
    typename FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits,
    jacobian_method>;

  //! Temporal local operator
  using TemporalLocalOperator = TemporalLocalOperatorDiffusionReactionCG<
    GridView,
    typename FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits,
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
  using BaseFEM = PDELab::
    PkLocalFiniteElementMap<typename G::LeafGridView, double, double, FEMorder>;

  static constexpr bool is_sub_model =
    not std::is_same_v<typename Grid::Traits::LeafGridView, GridView>;

  //! Finite element map
  using FEM =
    std::conditional_t<is_sub_model,
                       SubDomainLocalFiniteElementMap<BaseFEM, GridView>,
                       BaseFEM>;

  using OrderingTag = OT;
  static constexpr JacobianMethod jacobian_method = JM;

  //! Local operator
  using LocalOperator = LocalOperatorDiffusionReactionCG<
    GridView,
    typename FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits,
    jacobian_method>;

  //! Temporal local operator
  using TemporalLocalOperator = TemporalLocalOperatorDiffusionReactionCG<
    GridView,
    typename FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits,
    jacobian_method>;
};

template<class G, class GV, class OT, JacobianMethod JM>
struct ModelPkDiffusionReactionTraits<G, GV, 0, OT, JM>
  : public ModelP0DiffusionReactionTraits<G, GV, OT, JM>
{};

template<class G,
         class GV = typename G::Traits::LeafGridView,
         int PkOrder = 1,
         class OT = PDELab::EntityBlockedOrderingTag,
         JacobianMethod JM = JacobianMethod::Analytical>
struct ModelP0PkDiffusionReactionTraits
{
  using Grid = G;
  using GridView = GV;
  using FEMP0 = PDELab::
    P0LocalFiniteElementMap<typename Grid::ctype, double, Grid::dimension>;
  template<int order>
  using FEMPk = PDELab::
    PkLocalFiniteElementMap<typename G::LeafGridView, double, double, order>;

  static constexpr bool is_sub_model =
    not std::is_same_v<typename Grid::Traits::LeafGridView, GridView>;

  using FEMEntity = typename Grid::LeafGridView::template Codim<0>::Entity;
  //! Finite element map
  using FEM = std::conditional_t<
    is_sub_model,
    VariadicLocalFiniteElementMap<
      FEMEntity,
      SubDomainLocalFiniteElementMap<FEMP0, GridView>,
      SubDomainLocalFiniteElementMap<FEMPk<PkOrder>, GridView>>,
    VariadicLocalFiniteElementMap<FEMEntity, FEMP0, FEMPk<PkOrder>>>;

  using OrderingTag = OT;
  static constexpr JacobianMethod jacobian_method = JM;

  //! Local operator
  using LocalOperatorCG = LocalOperatorDiffusionReactionCG<
    GridView,
    typename FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits,
    jacobian_method>;

  using LocalOperatorFV = LocalOperatorDiffusionReactionFV<
    GridView,
    typename FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits,
    jacobian_method>;

  struct TestFunctor
  {
    template<class T>
    std::size_t operator()(const T& fe_v) const
    {
      return fe_v.type().isCube() ? 0 : 1;
    }
    template<class T0, class T1>
    std::size_t operator()(const T0& fe_u, const T1& fe_v) const
    {
      return fe_v.type().isCube() ? 0 : 1;
    }
  };

  using LocalOperatorVariadic =
    VariadicLocalOperator<TestFunctor, LocalOperatorFV, LocalOperatorCG>;

  using LocalOperator = LocalOperatorVariadic;

  //! Temporal local operator
  using TemporalLocalOperator = TemporalLocalOperatorDiffusionReactionCG<
    GridView,
    typename FEM::Traits::FiniteElement::Traits::LocalBasisType::Traits,
    jacobian_method>;
};

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
  using RF = typename FEM::Traits::FiniteElement::Traits::LocalBasisType::
    Traits::RangeFieldType;

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

  //! Constraints container
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;

  //! Local operator
  using LOP = typename Traits::LocalOperator;

  //! Temporal local operator
  using TLOP = typename Traits::TemporalLocalOperator;

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

  //! Jacobian operator
  using JO = Dune::PDELab::ISTLBackend_NOVLP_BCGS_SSORk<GOI>;

  using DataHandler = PDELab::vtk::
    DGFTreeCommonData<const GFS, const X, PDELab::vtk::DefaultPredicate, GV>;

  using ComponentLFS = typename PDELab::LocalFunctionSpace<GFS>::ChildType;

  using ComponentGridFunction =
    PDELab::vtk::DGFTreeLeafFunction<ComponentLFS, DataHandler, GV>;

public:

  //! Model state structure
  using State = Dune::Copasi::ModelState<Grid, GFS, X>;

  //! Constant model state structure
  using ConstState = Dune::Copasi::ConstModelState<Grid, GFS, X>;

  //! Grid operator type
  using GridOperator = GOI;

  //! Jacobian operator type
  using JacobianOperator = JO;

  /**
   * @brief      Constructs the model
   * @todo       Make a seccion describing the requirements of the confi file
   * @details    The model will be constructed according to the stages included in the
   *             setup policy
   *
   * @param[in]  grid          The grid
   * @param[in]  config        The configuration file
   * @param[in]  grid_view     The grid view to operate with
   * @param[in]  setup_policy  Policy to setup model
   */
  ModelDiffusionReaction(std::shared_ptr<Grid> grid,
                         const Dune::ParameterTree& config,
                         GV grid_view,
                         BitFlags<ModelSetup::Stages> setup_policy =
                           BitFlags<ModelSetup::Stages>::all_flags());

  /**
   * @brief      Constructs the model
   * @warning    This constructor only is available if the grid view
   *             is the leaf grid view of the templated grid
   * @todo       Make a seccion describing the requirements of the confi file
   * @details    The model will be constructed according to the stages included in the
   *             setup policy
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
                         BitFlags<ModelSetup::Stages> setup_policy =
                           BitFlags<ModelSetup::Stages>::all_flags())
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
   * @brief      Get mutable model state
   *
   * @return     Model state
   */
  State state()
  {
    return _state;
  }

  /**
   * @brief      Get constat model state
   *
   * @return     Constant model state
   */
  ConstState const_state() const
  {
    return ConstState{_state};
  }

  /**
   * @brief      Get constat model state
   *
   * @return     Constant model state
   */
  ConstState state() const { return const_state(); }

  /**
   * @brief Sets model internal state
   *
   * @param state  Valid model state
   */
  void set_state(const State& state)
  {
    if (not state)
      DUNE_THROW(InvalidStateException,"State must be valid");
    _state = state;
    if (_state.grid_function_space != state.grid_function_space)
    {
      if (_constraints)
        setup_constraints();
      if (_local_operator)
        setup_local_operator();
      if (_grid_operator)
        setup_grid_operator();
      if (_jacobian_operator)
        setup_jacobian_operator();
    }
  }

  /**
   * @brief      Sets the initial state of the model
   *
   * @param[in]  model_config  A parameter tree with 'initial' and optionally
   * 'data' subsections
   */
  template<class GFGridView>
  static auto get_muparser_initial(const ParameterTree& model_config,
                                   const GFGridView& gf_grid_view,
                                   bool compile = true);

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

  /**
   * @brief      Gets a grid function for a given component, and a state.
   * @details    The resulting grid function is persistent w.r.t the grid.
   *             This means that the grid function will be valid and will
   * contain exaclty the same data even if the model is modified in any form.
   * The only exception to this is when the grid is modified.
   *
   * @param[in]  state  The model state
   * @param[in]  comp    The component
   *
   * @return     The grid function.
   */
  static std::shared_ptr<ComponentGridFunction> get_grid_function(
    const ConstState& state,
    std::size_t comp);

  /**
   * @brief      Gets a grid function for a given component at the current state
   * of the model.
   * @details    The resulting grid function is persistent w.r.t the grid.
   *             This means that the grid function will be valid and will
   * contain exaclty the same data even if the model is modified in any form.
   * The only exception to this is when the grid is modified.
   *
   * @param[in]  comp    The component
   *
   * @return     The grid function.
   */
  std::shared_ptr<ComponentGridFunction> get_grid_function(
    std::size_t comp) const;

  /**
   * @brief      Gets a grid function for each component.
   * @details    The resulting grid functions are persistent w.r.t the grid.
   *             This means that the grid functions will be valid and will
   * contain exaclty the same data even if the model is modified in any form.
   * The only exception to this is when the grid is modified.
   *
   * @param[in]  state  The model state
   *
   * @return     The grid functions.
   */
  static std::vector<std::shared_ptr<ComponentGridFunction>> get_grid_functions(
    const ConstState& state);

  /**
   * @brief      Gets a grid function for each component at the current state of
   * the model.
   * @details    The resulting grid functions are persistent w.r.t the grid.
   *             This means that the grid functions will be valid and will
   * contain exaclty the same data even if the model is modified in any form.
   * The only exception to this is when the grid is modified.
   *
   * @param[in]  state  The model state
   *
   * @return     The grid functions.
   */
  std::vector<std::shared_ptr<ComponentGridFunction>> get_grid_functions()
    const;


  //! warning this is not completely const correct. grid operators may modify local operators and thus the model on certain calls
  std::shared_ptr<GridOperator> get_grid_operator() const
  {
    return _grid_operator;
  }

  std::shared_ptr<JacobianOperator> get_jacobian_operator() const
  {
    return _jacobian_operator;
  }

protected:

protected:
  /**
   * @brief      Setup model internals
   * @details    Sets up the model with specific options passdes on the policy
   *
   * @param[in]  setup_policy  The setup policy
   */
  void setup(BitFlags<ModelSetup::Stages> setup_policy);

private:

  auto setup_component_grid_function_space(const std::string&) const;
  void setup_grid_function_space();
  void setup_coefficient_vector();
  void setup_initial_condition();
  void setup_constraints();
  void setup_local_operator();
  void setup_grid_operator();
  void setup_jacobian_operator();
  void setup_vtk_writer();

  static auto get_data_handler(const ConstState& state);

  using ModelBase::_logger;

private:
  ParameterTree _config;
  std::string _compartment_name;
  GV _grid_view;
  State _state;
  std::shared_ptr<Grid> _grid;

  std::unique_ptr<CC> _constraints;
  std::shared_ptr<LOP> _local_operator;
  std::shared_ptr<TLOP> _temporal_local_operator;
  std::shared_ptr<GOS> _spatial_grid_operator;
  std::shared_ptr<GOT> _temporal_grid_operator;
  std::shared_ptr<GOI> _grid_operator;
  std::shared_ptr<JO> _jacobian_operator;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH
