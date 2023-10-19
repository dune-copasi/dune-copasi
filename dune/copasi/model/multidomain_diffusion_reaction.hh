#ifndef DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH

#ifndef DUNE_COPASI_HAVE_MEMBRANE_SPACE
#define DUNE_COPASI_HAVE_MEMBRANE_SPACE 0
#endif

#include <dune/copasi/common/enum.hh>
#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/grid/multidomain_entity_transformation.hh>
#include <dune/copasi/local_operator/diffusion_reaction/multidomain.hh>
#include <dune/copasi/model/diffusion_reaction.cc>
#include <dune/copasi/model/diffusion_reaction.hh>

#include <dune/copasi/finite_element_map/intersection.hh>

// #include <dune/pdelab/backend/istl.hh>
// #include <dune/pdelab/backend/istl/novlpistlsolverbackend.hh>
// #include <dune/pdelab/gridfunctionspace/dynamicpowergridfunctionspace.hh>
// #include <dune/pdelab/gridfunctionspace/vtk.hh>


#include <dune/common/parametertree.hh>

#include <memory>

namespace Dune::Copasi {


template<Dune::Concept::Grid G, std::size_t Order = 1, bool MultiThread = 1>
struct ModelMultiDomainDiffusionReactionPkTraits : public ModelDiffusionPkReactionTraits<G, typename G::SubDomainGrid::LeafGridView, Order, MultiThread>
{
  using EntitySet = std::conditional_t<MultiThread, Assembler::EntitySetThreadPool<typename G::LeafGridView>, typename G::LeafGridView>;
  using MultiCompartmentMergingStrategy = Assembler::Lexicographic<false>;

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  using SkeletonFiniteElementMap = PDELab::fem::PkLocalFiniteElementMapBase<EntitySet,double,double,Order,EntitySet::dimension-1>;
  using MembraneSpeciesFiniteElementMap = IntersectionLocalFiniteElementMap<SkeletonFiniteElementMap, EntitySet>;
  using MultiCompartmentMembraneMergingStrategy = Assembler::Lexicographic<false>; // top level blocking
#endif

};


/**
 * @brief      Class for diffusion-reaction models in multigrid domains.
 *
 * @tparam     Traits  Class that define static policies on the model
 */
template<class Traits_>
class ModelMultiDomainDiffusionReaction
{
public:
  using Traits = Traits_;

private:
  using Grid = typename Traits::Grid;

  //! Type for the set of grid entities that define a compartment (e.g. a grid view)
  using CompartmentEntitySet = typename Traits::CompartmentEntitySet;

  using SubModel = ModelDiffusionReaction<Traits>;

  //! Grid view
  using EntitySet = typename Traits::EntitySet;

  //! SubDomain grid function space
  using CompartmentSpace = typename SubModel::CompartmentSpace;

  using CompartmentMergingStrategy = typename SubModel::CompartmentMergingStrategy;

  // Choose a merging strategy for the degrees of freedom on the species between different compartment
  using MultiCompartmentMergingStrategy = typename Traits::MultiCompartmentMergingStrategy;
  using MultiCompartmentSpace = Assembler::VectorDiscreteFunctionSpace<MultiCompartmentMergingStrategy, CompartmentSpace>;  // TODO: Use assembler GV to update persistent regions

  //! Membranes!
#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  static constexpr auto compartments_path = Assembler::multiIndex(Indices::_0);
  using SkeletonFiniteElementMap = typename Traits::SkeletonFiniteElementMap;
  using MembraneSpeciesFiniteElementMap = typename Traits::MembraneSpeciesFiniteElementMap;
  using MembraneSpeciesMergingStrategy = Assembler::EntityGrouping<EntitySet, false>;
  using MembraneSpeciesSpace = Assembler::LeafDiscreteFunctionSpace<MembraneSpeciesMergingStrategy, MembraneSpeciesFiniteElementMap, PDELab::NoConstraints>;

  using MembraneMergingStrategy = Assembler::EntityGrouping<EntitySet, CompartmentMergingStrategy::Blocked>;
  using MembraneSpace = Assembler::VectorDiscreteFunctionSpace<MembraneMergingStrategy, MembraneSpeciesSpace>;

  // membranes and compartments should always have same blocking 
  using MultiMembraneMergingStrategy = MultiCompartmentMergingStrategy;
  using MultiMembraneSpace = Assembler::VectorDiscreteFunctionSpace<MultiMembraneMergingStrategy, MembraneSpace>;

  using MultiCompartmentMembraneMergingStrategy = typename Traits::MultiCompartmentMembraneMergingStrategy;
  using MultiCompartmentMembraneSpace = Assembler::TupleDiscreteFunctionSpace<MultiCompartmentMembraneMergingStrategy, MultiCompartmentSpace, MultiMembraneSpace>;

  // Create the final space with the options above
  using Space = Assembler::DiscreteFunctionSpace<MultiCompartmentMembraneSpace, EntitySet>;
#else
  static constexpr auto compartments_path = Assembler::multiIndex();
  using Space = Assembler::DiscreteFunctionSpace<MultiCompartmentSpace, EntitySet>;
#endif

 // Choose a quantity type to represent the species (e.g. double)
  using SpeciesQuantity = double; // TODO use units

  // Choose a backend whose containers have the same value: SpeciesQuantity
  // (e.g. a flat space container may give Dune::BlockVector<SpeciesQuantity> or std::vector<SpeciesQuantity>)
  using CoefficientsBackend = Assembler::ISTLUniformBackend<SpeciesQuantity>;
  // Extract container according to the backend and the space blocking strategy
  // (e.g. Dune::BlockVector<Dune::FieldVector<SpeciesQuantity,2>> or std::vector<std::array<SpeciesQuantity,2>>)
  using Coefficients = typename Space::template Container<CoefficientsBackend>;

  // Choose a quantity type to represent the residual (e.g. double)
  using ResidualQuantity = double; // TODO use units
  using ResidualBackend = Assembler::ISTLUniformBackend<ResidualQuantity>;
  using Residual = typename Space::template Container<ResidualBackend>;

  using Time = double;

  using StepOperator = Assembler::OneStepOperator<Coefficients, Time, Time>;


//   //! Local operator
//   using LOP = LocalOperatorMultiDomainDiffusionReaction<
//     Grid,
//     typename Traits::SubModelTraits::LocalOperator,
//     Traits::jacobian_method_interface>;

//   //! Temporal local operator
//   using TLOP = TemporalLocalOperatorMultiDomainDiffusionReaction<
//     Grid,
//     typename Traits::SubModelTraits::TemporalLocalOperator,
//     Traits::jacobian_method_interface>;

//   //! Matrix backend
//   using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;

//   //! Spatial grid operator
//   using GOS =
//     Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, RF, RF, RF, CC, CC>;

//   //! Temporal grid operator
//   using GOT =
//     Dune::PDELab::GridOperator<GFS, GFS, TLOP, MBE, RF, RF, RF, CC, CC>;

//   //! Instationary grid operator
//   using GOI = Dune::PDELab::OneStepGridOperator<GOS, GOT>;

//   //! Entity transformation between grids
//   using EntityTransformation =
//     Dune::Copasi::MultiDomainEntityTransformation<Grid>;

//   using DataHandler =
//     PDELab::vtk::DGFTreeCommonData<const GFS,
//                                    const X,
//                                    PDELab::vtk::DefaultPredicate,
//                                    SubDomainGridView,
//                                    EntityTransformation>;

//   using ComponentLFS =
//     typename PDELab::LocalFunctionSpace<GFS>::ChildType::ChildType;

//   using ComponentGridFunction = PDELab::vtk::
//     DGFTreeLeafFunction<ComponentLFS, DataHandler, SubDomainGridView>;

public:

  //! Model state structure
  using State = Dune::Copasi::ModelState<Grid, Space, Coefficients>;

  //! Constant model state structure
  using ConstState = typename State::Const;

//   //! Grid operator type
//   using StationaryGridOperator = GOS;

//   //! Grid operator type
//   using InstationaryGridOperator = GOI;

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

//   /**
//    * @brief      Destroys the object.
//    */
//   ~ModelMultiDomainDiffusionReaction();

  /**
   * @brief      Get mutable model state
   *
   * @return     Model state
   */
  State& state() { return _state; }

  /**
   * @brief      Get the model state
   *
   * @return     Model state
   */
  ConstState state() const { return _state; }

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
  void interpolate(State& state,
    const std::map<std::array<std::string, 2>, GF>& initial);

//   /**
//    * @brief      Gets a grid function for a given component, a sub domain, and a
//    * state.
//    * @details    The resulting grid function is persistent w.r.t the grid.
//    *             This means that the grid function will be valid and will
//    * contain exaclty the same data even if the model is modified in any form.
//    * The only exception to this is when the grid is modified.
//    *
//    * @param[in]  state  The model state
//    * @param[in]  domain  The domain
//    * @param[in]  comp    The component
//    *
//    * @return     The grid function.
//    */
//   std::shared_ptr<ComponentGridFunction> get_grid_function(
//     const ConstState& state,
//     std::size_t domain,
//     std::size_t comp) const;

//   /**
//    * @brief      Gets a grid function for a given component, and a sub domain at
//    * the current state of the model.
//    * @details    The resulting grid function is persistent w.r.t the grid.
//    *             This means that the grid function will be valid and will
//    * contain exaclty the same data even if the model is modified in any form.
//    * The only exception to this is when the grid is modified.
//    *
//    * @param[in]  domain  The domain
//    * @param[in]  comp    The component
//    *
//    * @return     The grid function.
//    */
//   std::shared_ptr<ComponentGridFunction> get_grid_function(
//     std::size_t domain,
//     std::size_t comp) const;

//   /**
//    * @brief      Gets a grid function for each component, and each sub domain.
//    * @details    The resulting grid functions are persistent w.r.t the grid.
//    *             This means that the grid functions will be valid and will
//    * contain exaclty the same data even if the model is modified in any form.
//    * The only exception to this is when the grid is modified.
//    *
//    * @param[in]  state  The model state
//    *
//    * @return     The grid functions.
//    */
//   std::vector<std::vector<std::shared_ptr<ComponentGridFunction>>>
//   get_grid_functions(const ConstState& state) const;

//   /**
//    * @brief      Gets a grid function for each component, and each sub domain at
//    * the current state of the model.
//    * @details    The resulting grid functions are persistent w.r.t the grid.
//    *             This means that the grid functions will be valid and will
//    * contain exaclty the same data even if the model is modified in any form.
//    * The only exception to this is when the grid is modified.
//    *
//    * @param[in]  state  The model state
//    *
//    * @return     The grid functions.
//    */
//   std::vector<std::vector<std::shared_ptr<ComponentGridFunction>>>
//   get_grid_functions() const;

  std::unique_ptr<StepOperator> make_step_operator(const ConstState& state, const ParameterTree&) const
  {
    return setup_step_operator(state);
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
  void setup_grid_function_space(State& state);
  void setup_coefficient_vector(State& state);
  void setup_initial_condition(State& state);
//   void setup_constraints();
//   auto setup_local_operator(std::size_t) const;
//   void setup_local_operator();
  std::unique_ptr<StepOperator> setup_step_operator(const ConstState&) const;
//   void setup_solvers();
  void setup_vtk_writer(State& state);

//   auto get_data_handler(const ConstState&) const;

//   using ModelBase::_logger;

private:

  ParameterTree _config;
  Logging::Logger _logger;
  State _state;

  std::shared_ptr<StepOperator> _step_operator;

#if DUNE_COPASI_HAVE_MEMBRANE_SPACE
  using SubDomainIndex = typename Grid::SubDomainIndex;
  std::vector<std::array<SubDomainIndex,2>> _membrane_map;
#endif
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MULTIDOMAIN_DIFFUSION_REACTION_HH
