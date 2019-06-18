#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH

#include <dune/copasi/model_base.hh>

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/backend/istl.hh>

#include <dune/grid/uggrid.hh>

#include <memory>

namespace Dune::Copasi {

/**
 * @brief      Class for diffusion-reaction models.
 *
 * @tparam     components  Number of components
 */
template<int components>
class ModelDiffusionReaction : public ModelBase {

  //! World dimension
  static constexpr int dim = 2;

  //! Grid type
  using Grid = Dune::UGGrid<dim>;

  //! Grid view
  using GV = typename Grid::LeafGridView;

  //! Domain field
  using DF = typename Grid::ctype;

  //! Range field
  using RF = double;

  //! Finite element map
  using FEM = PDELab::QkLocalFiniteElementMap<GV,DF,RF,1>;

  //! Constraints builder
  using CON = PDELab::ConformingDirichletConstraints;

  //! Root vector backend
  using RVBE = PDELab::ISTL::VectorBackend<>;

  //! Root grid function space
  using RGFS = PDELab::GridFunctionSpace<GV,FEM,CON,RVBE>;

  //! Vector backend
  using VBE = RVBE;

  //! Ordering tag
  using OrderingTag = PDELab::LexicographicOrderingTag;

  //! Grid function space
  using GFS = PDELab::PowerGridFunctionSpace<RGFS,components,VBE,OrderingTag>;

  //! Coefficient vector
  using X = PDELab::Backend::Vector<GFS,RF>;

  //! Constraints container
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;

public:

  ModelDiffusionReaction(std::shared_ptr<Grid> grid, 
                         const Dune::ParameterTree& config);

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

protected:

  /**
   * @brief      Setup operator for next time step
   */
  void operator_setup();

private:
  std::shared_ptr<Grid> _grid;
  std::shared_ptr<GFS> _gfs;

  GV _grid_view;
};

} // Dune::Copasi namespace

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_HH