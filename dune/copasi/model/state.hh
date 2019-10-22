#ifndef DUNE_COPASI_MODEL_STATE_HH
#define DUNE_COPASI_MODEL_STATE_HH

#include <dune/common/exceptions.hh>

#include <memory>

namespace Dune::Copasi {

/// forward declaration of model state
template<class G, class GFS, class X>
struct ModelState;

/**
 * @brief      Constant model state container
 *
 * @tparam     G     Grid type
 * @tparam     GFS   Grid function space type
 * @tparam     X     Coefficients vector
 */
template<class G, class GFS, class X>
struct ConstModelState
{
  //! Grid
  using Grid = G;

  //! Grid function space
  using GridFunctionSpace = GFS;

  //! Coefficients vector
  using Coefficients = X;

  std::shared_ptr<const Grid> grid;
  std::shared_ptr<const GridFunctionSpace> grid_function_space;
  std::shared_ptr<const Coefficients> coefficients;
  const double time;

  // ConstModelState() = default;

  ConstModelState(const ModelState<G, GFS, X>& model_state)
    : grid(model_state.grid)
    , grid_function_space(model_state.grid_function_space)
    , coefficients(model_state.coefficients)
    , time(model_state.time)
  {}

  /**
   * @brief      Write the model state into a file
   *
   * @param[in]  file_name  The file name
   */
  void write(const std::string& file_name) const
  {
    DUNE_THROW(Dune::NotImplemented, "ModelState writer is not implemented");
  }
};

/**
 * @brief      Model state container
 *
 * @tparam     G     Grid type
 * @tparam     GFS   Grid function space type
 * @tparam     X     Coefficients vector
 */
template<class G, class GFS, class X>
struct ModelState
{
  //! Grid
  using Grid = G;

  //! Grid function space
  using GridFunctionSpace = GFS;

  //! Coefficients vector
  using Coefficients = X;

  std::shared_ptr<Grid> grid;
  std::shared_ptr<GridFunctionSpace> grid_function_space;
  std::shared_ptr<Coefficients> coefficients;
  double time;

  // ModelState() = default;

  /**
   * @brief      Write the model state into a file
   *
   * @param[in]  file_name  The file name
   */
  void write(const std::string& file_name) const
  {
    ConstModelState<G, GFS, X> const_state(*this);
    const_state.write(file_name);
  }

  /**
   * @brief      Read a model state from a file
   *
   * @param[in]  file_name  The file name
   */
  void read(const std::string& file_name)
  {
    DUNE_THROW(Dune::NotImplemented, "ModelState reader is not implemented");
  }

  /**
   * @brief      implicit conversion to constant state
   */
  operator ConstModelState<G, GFS, X>()
  {
    return ConstModelState<G, GFS, X>{
      grid, grid_function_space, coefficients, time
    };
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_STATE_HH