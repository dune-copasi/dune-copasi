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
 * @tparam     S     Space type
 * @tparam     X     Coefficients type
 */
template<class G, class S, class C, class TP = double>
struct ModelState
{
  //! Grid
  using Grid = G;

  //! Grid function space
  using Space = S;

  //! Coefficients vector
  using Coefficients = C;

  //! Point in time
  using TimePoint = TP;

  //! Non mutable version of this object
  using Const = ModelState<const G, const S,const C, const TP>;

private:
  std::shared_ptr<Grid> _grid;
  std::optional<Space> _space;
  std::optional<Coefficients> _coefficients;
  TimePoint _time_point;

public:
  std::function<void(const ModelState&, const fs::path&, bool)> writer;
  std::function<void(ModelState&, const fs::path&)> reader;

  ModelState(const std::shared_ptr<Grid>& grid,
             const std::optional<Space>& space,
             const std::optional<Coefficients>& coefficients,
             TimePoint time_point = std::make_shared<TimePoint>(std::numeric_limits<TimePoint>::quiteNaN()))
    : _grid(grid)
    , _space(space)
    , _coefficients(coefficients)
    , _time_point(time_point)
  {}

  ModelState() = default;
  ~ModelState() = default;

   operator Const() const {
    return {_grid, _space, _coefficients, _time_point};
   }

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

  const Space& space() const
  {
    assert(_space);
    return *_space;
  }

  Space& space()
  {
    assert(_space);
    return *_space;
  }

  TimePoint& time_point()
  {
    return _time_point;
  }

  TimePoint time_point() const
  {
    return _time_point;
  }

  void set_time_point(TimePoint time_point)
  {
    _time_point = time_point;
  }

  void set_space(const Space& space){
    _space.emplace(space);
  }

  void set_space(Space&& space){
    _space.emplace(std::move(space));
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

  void set_coefficients(const Coefficients& coefficients) {
    _coefficients.emplace(coefficients);
  }

  void set_coefficients(Coefficients&& coefficients) {
    _coefficients.emplace(std::move(coefficients));
  }

  /**
   * @brief      implicit conversion to bool
   * @details    A vaild state contains a grid, a grid function space and
   *             coefficients
   */
  operator bool() const
  {
    return _grid and _space and _coefficients;
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_STATE_HH
