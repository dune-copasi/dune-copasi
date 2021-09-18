#ifndef DUNE_COPASI_MODEL_STATE_HH
#define DUNE_COPASI_MODEL_STATE_HH

#include <dune/copasi/common/filesystem.hh>

#include <dune/pdelab/backend/interface.hh>

#include <dune/common/exceptions.hh>

#include <limits>
#include <memory>

namespace Dune::Copasi {

/**
 * @brief      Model state container
 *
 * @tparam     G     Grid type
 * @tparam     GFS   Grid function space type
 * @tparam     RF    Range field
 */
template<class X>
struct PDELabState
{

  //! Grid function space
  using GridFunctionSpace = typename X::GridFunctionSpace;

  //! Grid
  using Grid = typename GridFunctionSpace::Traits::GridView::Grid;

  //! Coefficients vector
  using Coefficients = X;

private:
  std::shared_ptr<Grid> _grid;
  std::shared_ptr<GridFunctionSpace> _grid_function_space;
  std::shared_ptr<Coefficients> _coefficients;

public:
  double time = std::numeric_limits<double>::quiet_NaN();
  std::function<void(const PDELabState&, const fs::path&, bool)> writer;
  std::function<void(PDELabState&, const fs::path&)> reader;

  PDELabState(const std::shared_ptr<Grid>& grid,
              const std::shared_ptr<GridFunctionSpace>& grid_function_space,
              const std::shared_ptr<Coefficients>& coefficients)
    : _grid(grid)
    , _grid_function_space(grid_function_space)
    , _coefficients(coefficients)
  {}

  ~PDELabState() = default;

  /**
   * @brief      Write the model state into a file
   *
   * @param[in]  path     The path to write file
   * @param[in]  append   True if write should apppend file to older state
   *                      write
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

  const Grid& grid() const { return *_grid; }

  Grid& grid() { return *_grid; }

  std::shared_ptr<const Grid> grid_storage() const { return _grid; }

  std::shared_ptr<Grid>& grid_storage() { return _grid; }

  const GridFunctionSpace& space() const
  {
    assert(_grid_function_space);
    return *_grid_function_space;
  }

  GridFunctionSpace& space()
  {
    assert(_grid_function_space);
    return *_grid_function_space;
  }

  std::shared_ptr<const GridFunctionSpace> space_storage() const
  {
    return _grid_function_space;
  }

  std::shared_ptr<GridFunctionSpace>& space_storage()
  {
    return _grid_function_space;
  }

  const Coefficients& coefficients() const
  {
    assert(_coefficients);
    return *_coefficients;
  }

  Coefficients& coefficients()
  {
    assert(_coefficients);
    return *_coefficients;
  }

  std::shared_ptr<const Coefficients> coefficients_storage() const
  {
    return _coefficients;
  }

  std::shared_ptr<Coefficients>& coefficients_storage()
  {
    return _coefficients;
  }

  /**
   * @brief      implicit conversion to bool
   * @details    A vaild state contains a grid, a grid function space and
   *             coefficients
   */
  operator bool() const
  {
    return _grid and _grid_function_space and _coefficients;
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_STATE_HH
