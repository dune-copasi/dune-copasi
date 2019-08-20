#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH

#include <dune/copasi/concepts/grid.hh>
#include <dune/copasi/grid_function_writer.hh>
#include <dune/copasi/local_operator.hh>
#include <dune/copasi/model_base.hh>
#include <dune/copasi/model_state.hh>
#include <dune/copasi/multidomain_local_finite_element_map.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/function/discretegridviewfunction.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/newton/newton.hh>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <memory>

namespace Dune::Copasi {

/**
 * @brief      Class for diffusion-reaction models.
 *
 * @tparam     components  Number of components
 * @tparam     Param       Parameterization class
 */
template<class Grid, class GridView>
class ModelDiffusionReaction : public ModelBase
{
  // Check templates
  static_assert(Concept::isGrid<Grid>(), "Provided and invalid grid");

  //! World dimension
  static constexpr int dim = 2;

  //! Polynomial order
  static constexpr int order = 1;

  //! Grid view
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
  using FEM = MultiDomainLocalFiniteElementMap<BaseFEM, GridView>;

  //! Constraints builder
  using CON = PDELab::ConformingDirichletConstraints;

  //! Leaf vector backend
  using LVBE = PDELab::ISTL::VectorBackend<>;

  //! Leaf grid function space
  using LGFS = PDELab::GridFunctionSpace<GV, FEM, CON, LVBE>;

  //! Vector backend
  using VBE = LVBE;

  //! Ordering tag
  using OrderingTag = PDELab::LexicographicOrderingTag;

  //! Grid function space
  using GFS = LGFS;

  //! Coefficient vector
  using X = PDELab::Backend::Vector<GFS, RF>;

  //! Constraints container
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;

  //! Local operator
  using LOP = LocalOperatorDiffusionReaction<GV, FE>;

  //! Temporal local operator
  using TLOP = TemporalLocalOperatorDiffusionReaction<GV, FE>;

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
  using W = Dune::VTKWriter<GV>;

  //! Sequential writer
  using SW = Dune::Copasi::GridFunctionVTKSequenceWriter<GV>;

  // //! Discrete grid function
  // using DGF = Dune::Copasi::DiscreteGridFunction<GFS, X>;

public:
  //! Model state structure
  using State = Dune::Copasi::ModelState<Grid, GFS, X>;

  //! Constant model state structure
  using ConstState = Dune::Copasi::ModelState<const Grid, const GFS, const X>;

  /**
   * @brief      Constructs the model
   *
   * @param[in]  grid    The grid
   * @param[in]  config  The configuration file
   */
  ModelDiffusionReaction(std::shared_ptr<Grid> grid,
                         GV grid_view,
                         const Dune::ParameterTree& config);

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
  State state() { return _state; };

  /**
   * @brief      Get the model state
   *
   * @return     Model state
   */
  ConstState state() const { return _state; };
  ;

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
  /**
   * @brief      Setup operator for next time step
   */
  void operator_setup();

  using ModelBase::_logger;

private:
  std::size_t _components;
  std::size_t _dof_per_component;
  ParameterTree _config;
  GV _grid_view;
  State _state;
  std::shared_ptr<FEM> _finite_element_map;
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

  std::shared_ptr<W> _writer;
  std::shared_ptr<SW> _sequential_writer;

  std::shared_ptr<X>& _x;     //! reference to coefficients pointer
  std::shared_ptr<GFS>& _gfs; //! reference to grid function space pointer
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH