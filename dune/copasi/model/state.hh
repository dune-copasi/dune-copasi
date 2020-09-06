#ifndef DUNE_COPASI_MODEL_STATE_HH
#define DUNE_COPASI_MODEL_STATE_HH

#include <dune/common/exceptions.hh>

#include <dune/copasi/common/filesystem.hh>
#include <limits>
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
  using Grid = const G;

  //! Grid function space
  using GridFunctionSpace = const GFS;

  //! Coefficients vector
  using Coefficients = const X;

  std::shared_ptr<const Grid> grid;
  std::shared_ptr<const GridFunctionSpace> grid_function_space;
  std::shared_ptr<const Coefficients> coefficients;
  const double time = std::numeric_limits<double>::quiet_NaN();
  std::function<
    void(const ConstModelState&, const fs::path&, bool)>
    writer;
  std::function<void(ModelState<G, GFS, X>&, const fs::path&)>
    reader;

  /**
   * @brief      Constructor from non-const state
   *
   * @param[in]  model_state     Non-const state
   */
  ConstModelState(const ModelState<G, GFS, X>& model_state)
    : grid(model_state.grid)
    , grid_function_space(model_state.grid_function_space)
    , coefficients(model_state.coefficients)
    , time(model_state.time)
    , writer(model_state.writer)
    , reader(model_state.reader)
  {}

  /**
   * @brief      Write the model state into a file
   *
   * @param[in]  path     The path to write file
   * @param[in]  append   True if write should apppend file to older state write
   */
  void write(const fs::path& path, bool append) const
  {
    if (not writer)
      DUNE_THROW(Dune::InvalidStateException, "ModelState writer is not setup");
    writer(*this, path, append);
  }

  /**
   * @brief      implicit conversion to bool
   * @details    A vaild state contains a grid, a grid function space and
   *             coefficients
   */
  operator bool() const
  {
    return grid and grid_function_space and coefficients;
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
  double time = std::numeric_limits<double>::quiet_NaN();
  std::function<
    void(const ConstModelState<G, GFS, X>&, const fs::path&, bool)>
    writer;
  std::function<void(ModelState&, const fs::path&)> reader;

  /**
   * @brief      Write the model state into a file
   *
   * @param[in]  path     The path to write file
   * @param[in]  append   True if write should apppend file to older state write
   */
  void write(const fs::path& path, bool append) const
  {
    if (not writer)
      DUNE_THROW(Dune::InvalidStateException, "ModelState writer is not setup");
    writer(*this, path, append);
  }

  /**
   * @brief      Read a model state from a file
   *
   * @param[in]  path  The file path
   */
  void read(const fs::path& path)
  {
    if (not reader)
      DUNE_THROW(Dune::InvalidStateException, "ModelState reader is not setup");
    reader(*this, path);
  }

  /**
   * @brief      implicit conversion to constant state
   */
  operator ConstModelState<G, GFS, X>()
  {
    return ConstModelState<G, GFS, X>{
      grid, grid_function_space, coefficients, time, writer
    };
  }

  /**
   * @brief      implicit conversion to bool
   * @details    A vaild state contains a grid, a grid function space and
   *             coefficients
   */
  operator bool() const
  {
    return grid and grid_function_space and coefficients;
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_STATE_HH
