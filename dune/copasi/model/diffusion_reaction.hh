#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH

#include <dune/copasi/finite_element/p0.hh>
#include <dune/copasi/finite_element/pk.hh>
#include <dune/copasi/finite_element_map/p0.hh>
#include <dune/copasi/finite_element_map/pk.hh>
#include <dune/copasi/finite_element_map/subdomain.hh>
#include <dune/copasi/finite_element_map/variadic.hh>
#include <dune/copasi/finite_element_map/virtual.hh>
#include <dune/copasi/local_operator/diffusion_reaction/continuous_galerkin.hh>
// #include <dune/copasi/local_operator/diffusion_reaction/finite_volume.hh>
#include <dune/copasi/model/state.hh>

#include <dune/assembler/discrete_function_space/backend_istl.hh>
#include <dune/assembler/discrete_function_space/discrete_function_space.hh>
#include <dune/assembler/discrete_function_space/merging_strategy.hh>
#include <dune/assembler/discrete_function_space/discrete_function.hh>
#include <dune/assembler/runge_kutta/shu_osher_tableau.hh>
#include <dune/assembler/runge_kutta/operator.hh>
#include <dune/assembler/common/entityset_threadpool.hh>
#include <dune/assembler/common/entityset_tbb.hh>

#include <dune/pdelab/constraints/noconstraints.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>

#include <dune/grid/concepts/grid.hh>
#include <dune/grid/concepts/gridview.hh>

#if HAVE_UNITS
#include <units/isq/si/amount_of_substance.h>
#include <units/isq/si/time.h>
#include <units/quantity_io.h>
#include <units/math.h>

namespace Dune {
  template<class T>
  struct IsNumber;

  // allow dune to treat units as numbers
  template<units::Dimension D, units::UnitOf<D> U, units::Representation Rep>
  struct IsNumber<units::quantity<D,U,Rep>> : std::true_type {};


  namespace FloatCmp {

    template<class T> struct EpsilonType;

    template<units::Dimension D, units::UnitOf<D> U, units::Representation Rep>
    struct EpsilonType<units::quantity<D,U,Rep>> {
      using Type = Rep;
    };
  }
}

#endif

namespace Dune::Copasi {

template<Dune::Concept::Grid G, Dune::Concept::EntitySet<0> ES, std::size_t Order = 1, bool MultiThread = true>
struct ModelDiffusionPkReactionTraits {
  using Grid = G;
  using CompartmentEntitySet = std::conditional_t<MultiThread, Assembler::EntitySetThreadPool<ES>, ES>;
  static_assert(Order != 0);
  using SpeciesFiniteElementMap = PDELab::PkLocalFiniteElementMap<CompartmentEntitySet, double, double, Order>;

  using CompartmentMergingStrategy = Assembler::EntityGrouping<CompartmentEntitySet, false>;
};

/**
 * @brief      Class for diffusion-reaction models.
 *
 * @tparam     Traits  Class that define static policies on the model
 */
template<class Traits>
class ModelDiffusionReaction
{
public:
  // A dune grid type
  using Grid = typename Traits::Grid;

  //! Type for the set of grid entities that define a compartment (e.g. a grid view)
  using CompartmentEntitySet = typename Traits::CompartmentEntitySet;

  // Choose a mapping that gives an finite element for each entity in the grid
  using SpeciesFiniteElementMap = typename Traits::SpeciesFiniteElementMap;
  // Choose a merging strategy for the degrees of freedom on each finite element
  // Degrees of freedom will be grouped by entity affinity and won't be blocked
  using SpeciesMergingStrategy = Assembler::EntityGrouping<CompartmentEntitySet, false>;
  // Constrcut space for individual biochemical species
  using SpeciesSpace = Assembler::LeafDiscreteFunctionSpace<SpeciesMergingStrategy, SpeciesFiniteElementMap, PDELab::NoConstraints>;

  // Choose a merging strategy for the degrees of freedom on the species within a compartment
  // Degrees of freedom will be grouped by entity affinity and will be blocked
  // Each block will contain degrees of freedom for all the species within a entity
  using CompartmentMergingStrategy = typename Traits::CompartmentMergingStrategy;
  // Compose several species that live within a compartment
  using CompartmentSpace = Assembler::VectorDiscreteFunctionSpace<CompartmentMergingStrategy, SpeciesSpace>;
  // Choose a quantity type to represent the species (e.g. double)
  using SpeciesQuantityRepresentation = double;
#if HAVE_UNITS
  using AmountOfSubstanceQuantity = units::isq::si::amount_of_substance<units::isq::si::mole, SpeciesQuantityRepresentation>;
  using VolumeQuantity = typename CompartmentEntitySet::template Codim<0>::Entity::Geometry::Volume;
  using SpeciesQuantity = decltype(AmountOfSubstanceQuantity{}/VolumeQuantity{});
#else
  using SpeciesQuantity = SpeciesQuantityRepresentation;
#endif

  // Create the final space with the options above
  using Space = Assembler::DiscreteFunctionSpace<CompartmentSpace, CompartmentEntitySet>;

  // Choose a backend whose containers have the same value: SpeciesQuantity
  // (e.g. a flat space container may give Dune::BlockVector<SpeciesQuantity> or std::vector<SpeciesQuantity>)
  using CoefficientsBackend = Assembler::ISTLUniformBackend<SpeciesQuantity>;
  // Extract container according to the backend and the space blocking strategy
  // (e.g. Dune::BlockVector<Dune::FieldVector<SpeciesQuantity,2>> or std::vector<std::array<SpeciesQuantity,2>>)
  using Coefficients = typename Space::template Container<CoefficientsBackend>;

  using SpeciesGridFunction = Assembler::ScalarDiscreteFunction<Space, Coefficients>;

  using TimeRepresentation = double;
#if HAVE_UNITS
  using TimeQuantity = units::isq::si::time<units::isq::si::second, TimeRepresentation>;
#else
  using TimeQuantity = TimeRepresentation;
#endif

  // Choose a quantity type to represent the residual (e.g. double)
  using ResidualQuantity = SpeciesQuantity; // RK: dt in stiffnes numerator
  // using ResidualQuantity = decltype(SpeciesQuantity{}/TimeQuantity{}); // RK: dt in mass denominator
  using ResidualBackend = Assembler::ISTLUniformBackend<ResidualQuantity>;

  using Residual = typename Space::template Container<ResidualBackend>;

  using StepOperator = Assembler::OneStepOperator<Coefficients, TimeQuantity, TimeQuantity>;

  //! Model state structure
  using State = Dune::Copasi::ModelState<Grid, Space, Coefficients, TimeQuantity>;

  //! Constant model state structure
  using ConstState = typename State::Const;

  /**
   * @brief      Constructs the model
   * @todo       Make a seccion describing the requirements of the confi file
   * @details    The model will be constructed according to the stages included in the
   *             setup policy
   *
   * @param[in]  grid          The grid
   * @param[in]  config        The configuration file
   * @param[in]  compartment_entity_set     The grid view to operate with
   * @param[in]  setup_policy  Policy to setup model
   */
  ModelDiffusionReaction(std::shared_ptr<Grid> grid,
                         const Dune::ParameterTree& config,
                         const CompartmentEntitySet& compartment_entity_set,
                         BitFlags<ModelSetup::Stages> setup_policy =
                           BitFlags<ModelSetup::Stages>::all_flags());

  /**
   * @brief      Destroys the model
   */
  ~ModelDiffusionReaction();

  /**
   * @brief      Get mutable model state
   *
   * @return     Model state
   */
  State& state()
  {
    return _state;
  }

  /**
   * @brief      Get constat model state
   *
   * @return     Constant model state
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
  void interpolate(State&, const std::map<std::string, GF>& initial);

  static SpeciesGridFunction as_function(const ConstState& state, Assembler::Concept::MultiIndex auto species_path);

  SpeciesGridFunction as_function(Assembler::Concept::MultiIndex auto species_path) const;

  std::unique_ptr<StepOperator> make_step_operator(const ConstState& state, const ParameterTree&) const
  {
    return setup_step_operator(state);
  }

// protected:

  auto make_compartment_function_space(const CompartmentEntitySet&) const;
protected:
  /**
   * @brief      Setup model internals
   * @details    Sets up the model with specific options passdes on the policy
   *
   * @param[in]  setup_policy  The setup policy
   */
  void setup(BitFlags<ModelSetup::Stages> setup_policy, const CompartmentEntitySet&);

private:

  auto make_component_grid_function_space(const CompartmentEntitySet&, const std::string&) const;
  void setup_grid_function_space(State&, const CompartmentEntitySet&);
  void setup_coefficient_vector(State&);
  void setup_initial_condition(State&);
  void setup_vtk_writer(State&);
  std::unique_ptr<StepOperator> setup_step_operator(const ConstState&) const;

private:
  ParameterTree _config;
  Logging::Logger _logger;
  State _state;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH
