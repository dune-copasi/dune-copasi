#ifndef DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/grid/multidomain_entity_transformation.hh>
#include <dune/copasi/local_operator/diffusion_reaction/multidomain.hh>
#include <dune/copasi/model/base.hh>
#include <dune/copasi/model/diffusion_reaction.cc>
#include <dune/copasi/model/diffusion_reaction.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/backend/istl/novlpistlsolverbackend.hh>
#include <dune/pdelab/gridfunctionspace/dynamicpowergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>

#include <dune/common/parametertree.hh>

#include <memory>

namespace Dune::Copasi {

/**
 * @brief      Traits for diffusion reaction models in multigrid domains
 *
 * @tparam     G         Grid
 * @tparam     FEMorder  Order of the finite element method
 * @tparam     OT        PDELab ordering tag
 * @tparam     JM        Jacobian method
 */
template<class G,
         int FEMorder = 1,
         class OT = PDELab::EntityBlockedOrderingTag,
         JacobianMethod JM = JacobianMethod::Analytical>
struct ModelMultiDomainPkDiffusionReactionTraits
{
  using Grid = G;

  // Check grid template
  static_assert(Concept::isMultiDomainGrid<Grid>(),
                "Provided grid type is not a multidomain grid");

  static_assert(Concept::isSubDomainGrid<typename Grid::SubDomainGrid>());

  using OrderingTag = OT;

  static constexpr JacobianMethod jacobian_method = JM;

  using SubModelTraits =
    ModelPkDiffusionReactionTraits<Grid,
                                   typename Grid::SubDomainGrid::LeafGridView,
                                   FEMorder,
                                   OrderingTag,
                                   jacobian_method>;
};

template<class G,
         int FEMorder = 1,
         class OT = PDELab::EntityBlockedOrderingTag,
         JacobianMethod JM = JacobianMethod::Analytical>
struct ModelMultiDomainP0PkDiffusionReactionTraits
{
  using Grid = G;

  // Check grid template
  static_assert(Concept::isMultiDomainGrid<Grid>(),
                "Provided grid type is not a multidomain grid");

  static_assert(Concept::isSubDomainGrid<typename Grid::SubDomainGrid>());

  using OrderingTag = OT;

  static constexpr JacobianMethod jacobian_method = JM;

  using SubModelTraits =
    ModelP0PkDiffusionReactionTraits<Grid,
                                     typename Grid::SubDomainGrid::LeafGridView,
                                     FEMorder,
                                     OrderingTag,
                                     jacobian_method>;
};

/**
 * @brief      Class for diffusion-reaction models in multigrid domains.
 *
 * @tparam     Traits  Class that define static policies on the model
 */
template<class Traits_>
class ModelMultiDomainDiffusionReaction : public ModelBase
{
public:
  using Traits = Traits_;

private:
  using Grid = typename Traits::Grid;

  using OT = typename Traits::OrderingTag;

  static constexpr JacobianMethod JM = Traits::jacobian_method;

  using SubDomainGridView = typename Grid::SubDomainGrid::LeafGridView;

  using SubModel = ModelDiffusionReaction<typename Traits::SubModelTraits>;

  //! Grid view
  using GridView = typename Grid::LeafGridView;
  using GV = GridView;

  //! Domain field
  using DF = typename Grid::ctype;

  //! Range field
  using RF = double;

  //! Constraints builder
  using CON = PDELab::ConformingDirichletConstraints;

  //! Leaf vector backend
  using LVBE = PDELab::ISTL::VectorBackend<>;

  //! Leaf grid function space
  using LGFS = typename SubModel::LGFS;

  //! SubDomain grid function space
  using SDGFS = typename SubModel::GFS;

  //! Vector backend
  using VBE = typename LGFS::Traits::Backend;
  using GFS = PDELab::DynamicPowerGridFunctionSpace<SDGFS, VBE>;

  //! Coefficient vector
  using X = PDELab::Backend::Vector<GFS, RF>;

  //! Constraints container
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;

  //! Local operator
  using LOP = LocalOperatorMultiDomainDiffusionReaction<
    Grid,
    typename Traits::SubModelTraits::LocalOperator,
    JM>;

  //! Temporal local operator
  using TLOP = TemporalLocalOperatorMultiDomainDiffusionReaction<
    Grid,
    typename Traits::SubModelTraits::TemporalLocalOperator,
    JM>;

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

  //! Entity transformation between grids
  using EntityTransformation =
    Dune::Copasi::MultiDomainEntityTransformation<Grid>;

  using DataHandler =
    PDELab::vtk::DGFTreeCommonData<const GFS,
                                   const X,
                                   PDELab::vtk::DefaultPredicate,
                                   SubDomainGridView,
                                   EntityTransformation>;

  using ComponentLFS =
    typename PDELab::LocalFunctionSpace<GFS>::ChildType::ChildType;

  using ComponentGridFunction = PDELab::vtk::
    DGFTreeLeafFunction<ComponentLFS, DataHandler, SubDomainGridView>;

public:

  //! Model state structure
  using State = Dune::Copasi::ModelState<Grid, GFS, X>;

  //! Constant model state structure
  using ConstState = Dune::Copasi::ConstModelState<Grid, GFS, X>;

  //! Grid operator type
  using StationaryGridOperator = GOS;

  //! Grid operator type
  using InstationaryGridOperator = GOI;

  /**
   * @brief      Constructs a new instance.
   * @todo       Make a seccion describing the requirements of the confi file
   * @details    The model will be constructed according to the stages included in the
   *             setup policy
   *
   * @param[in]  grid    The grid
   * @param[in]  config  The configuration
   */
  ModelMultiDomainDiffusionReaction(std::shared_ptr<Grid> grid,
                                    const Dune::ParameterTree& config,
                                    BitFlags<ModelSetup::Stages> setup_policy =
                                      BitFlags<ModelSetup::Stages>::all_flags());

  /**
   * @brief      Destroys the object.
   */
  ~ModelMultiDomainDiffusionReaction();

  /**
   * @brief      Get mutable model state
   *
   * @return     Model state
   */
  State state() { return _state; }

  /**
   * @brief      Get constat model state
   *
   * @return     Constant model state
   */
  ConstState const_state() const { return _state; }

  /**
   * @brief      Get the model state
   *
   * @return     Model state
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
   * @details    The input vector of vectors should have the same size as the
   * number of domains variables in the model, and each vector for each
   * subdomain has to have the same size as the number of variables in the
   * compartment. Additionally, variables will be indepreted aphabetically
   * accodingly to the name set to othe input sections (e.g.
   * 'model.<compartment>.diffusion' section).
   *
   * @tparam     GF       A valid PDELab grid functions (see
   * @Concepts::PDELabGridFunction)
   * @param[in]  initial  Vector of vecotrs of grid functions, one for each
   * variable
   */
  template<class GF>
  void set_initial(
    const std::vector<std::vector<std::shared_ptr<GF>>>& initial);

  /**
   * @brief      Gets a grid function for a given component, a sub domain, and a
   * state.
   * @details    The resulting grid function is persistent w.r.t the grid.
   *             This means that the grid function will be valid and will
   * contain exaclty the same data even if the model is modified in any form.
   * The only exception to this is when the grid is modified.
   *
   * @param[in]  state  The model state
   * @param[in]  domain  The domain
   * @param[in]  comp    The component
   *
   * @return     The grid function.
   */
  std::shared_ptr<ComponentGridFunction> get_grid_function(
    const ConstState& state,
    std::size_t domain,
    std::size_t comp) const;

  /**
   * @brief      Gets a grid function for a given component, and a sub domain at
   * the current state of the model.
   * @details    The resulting grid function is persistent w.r.t the grid.
   *             This means that the grid function will be valid and will
   * contain exaclty the same data even if the model is modified in any form.
   * The only exception to this is when the grid is modified.
   *
   * @param[in]  domain  The domain
   * @param[in]  comp    The component
   *
   * @return     The grid function.
   */
  std::shared_ptr<ComponentGridFunction> get_grid_function(
    std::size_t domain,
    std::size_t comp) const;

  /**
   * @brief      Gets a grid function for each component, and each sub domain.
   * @details    The resulting grid functions are persistent w.r.t the grid.
   *             This means that the grid functions will be valid and will
   * contain exaclty the same data even if the model is modified in any form.
   * The only exception to this is when the grid is modified.
   *
   * @param[in]  state  The model state
   *
   * @return     The grid functions.
   */
  std::vector<std::vector<std::shared_ptr<ComponentGridFunction>>>
  get_grid_functions(const ConstState& state) const;

  /**
   * @brief      Gets a grid function for each component, and each sub domain at
   * the current state of the model.
   * @details    The resulting grid functions are persistent w.r.t the grid.
   *             This means that the grid functions will be valid and will
   * contain exaclty the same data even if the model is modified in any form.
   * The only exception to this is when the grid is modified.
   *
   * @param[in]  state  The model state
   *
   * @return     The grid functions.
   */
  std::vector<std::vector<std::shared_ptr<ComponentGridFunction>>>
  get_grid_functions() const;

  //! warning this is not completely const correct. grid operators may modify local operators and thus the model on certain calls
  std::shared_ptr<StationaryGridOperator> get_stationary_grid_operator() const
  {
    return _spatial_grid_operator;
  }

  //! warning this is not completely const correct. grid operators may modify local operators and thus the model on certain calls
  std::shared_ptr<InstationaryGridOperator> get_instationary_grid_operator() const
  {
    return _grid_operator;
  }

protected:
  /**
   * @brief      Setup model internals
   * @details    Sets up the model with specific options passdes on the policy
   *
   * @param[in]  setup_policy  The setup policy
   */
  void setup(BitFlags<ModelSetup::Stages> setup_policy);

private:
  void setup_grid_function_space();
  void setup_coefficient_vector();
  void setup_initial_condition();
  void setup_constraints();
  auto setup_local_operator(std::size_t) const;
  void setup_local_operator();
  void setup_grid_operator();
  void setup_solvers();
  void setup_vtk_writer();

  auto get_data_handler(const ConstState&) const;

  using ModelBase::_logger;

private:

  ParameterTree _config;
  GV _grid_view;
  State _state;
  std::shared_ptr<Grid> _grid;

  std::unique_ptr<CC> _constraints;
  std::shared_ptr<LOP> _local_operator;
  std::shared_ptr<TLOP> _temporal_local_operator;
  std::shared_ptr<GOS> _spatial_grid_operator;
  std::shared_ptr<GOT> _temporal_grid_operator;
  std::shared_ptr<GOI> _grid_operator;

  std::size_t _domains;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH
