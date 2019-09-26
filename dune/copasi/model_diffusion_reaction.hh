#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH

#include <dune/copasi/coefficient_mapper.hh>
#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/enum.hh>
#include <dune/copasi/grid_function_writer.hh>
#include <dune/copasi/local_operator.hh>
#include <dune/copasi/model_base.hh>
#include <dune/copasi/model_state.hh>
#include <dune/copasi/multidomain_local_finite_element_map.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/function/discretegridviewfunction.hh>
#include <dune/pdelab/gridfunctionspace/dynamicpowergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/newton/newton.hh>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <memory>

namespace Dune::Copasi {

template<class G,
         class GV = typename G::Traits::LeafGridView,
         int FEMorder = 1,
         class OT = PDELab::EntityBlockedOrderingTag,
         JacobianMethod JM = JacobianMethod::Analytical>
struct ModelDiffusionReactionTraits
{
  using Grid = G;
  using GridView = GV;
  using OrderingTag = OT;
  static constexpr JacobianMethod jacobian_method = JM;
  static constexpr int order = FEMorder;
};

/**
 * @brief      Class for diffusion-reaction models.
 *
 * @tparam     components  Number of components
 * @tparam     Param       Parameterization class
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

  //! World dimension
  static constexpr int dim = 2;

  //! Polynomial order
  static constexpr int order = Traits::order;

  //! Grid view
  using GV = typename Traits::GridView;

  //! Host grid view
  using HGV = typename Traits::Grid::LeafGridView;

  //! Domain field
  using DF = typename Grid::ctype;

  //! Range field
  using RF = double;

  //! Finite element
  using FE = Dune::PkLocalFiniteElement<DF, RF, dim, order>;

  //! Base finite element map
  using BaseFEM = PDELab::PkLocalFiniteElementMap<HGV, DF, RF, order>;

  //! Finite element map
  using FEM = MultiDomainLocalFiniteElementMap<BaseFEM, GV>;

  //! Constraints builder
  using CON = PDELab::ConformingDirichletConstraints;

  //! Entity set
  using ES = Dune::PDELab::NonOverlappingEntitySet<HGV>;

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
  using LOP = LocalOperatorDiffusionReaction<GV, FE, CM, JM>;

  //! Temporal local operator
  using TLOP = TemporalLocalOperatorDiffusionReaction<GV, FE, JM>;

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

public:
  /**
   * @brief      Constructs the model
   *
   * @param[in]  grid    The grid
   * @param[in]  config  The configuration file
   */
  ModelDiffusionReaction(std::shared_ptr<Grid> grid,
                         const Dune::ParameterTree& config,
                         GV grid_view,
                         ModelSetupPolicy setup_policy = ModelSetupPolicy::All);

  /**
   * @brief      Constructs the model
   *
   * @param[in]  grid    The grid
   * @param[in]  config  The configuration file
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
    std::map<std::size_t, ConstState> const_states(_states.begin(),
                                                   _states.end());
    return const_states;
  }

  /**
   * @brief      Get the model state
   *
   * @return     Model state
   */
  std::map<std::size_t, ConstState> states() const { return const_states(); }

  /**
   * @brief      Sets the state of the model
   *
   * @param[in]  input_state  The state to set in the model
   *
   * @tparam     T            Type of valid input states. Valid states are:
   *                          * Arithmetic values: Set all components with the
   *                            same value everywhere in the domain
   *                          * Field vector: Set each component with the values
   *                            in the field vector everywhere in the domain
   *                          * PDELab callable: A lambda function which returns
   *                            a field vector with the components state for
   *                            every position in the domain
   *                          * PDELab grid function: A function following the
   *                            PDELab grid function interface
   */
  template<class T>
  void set_state(const T& input_state);

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