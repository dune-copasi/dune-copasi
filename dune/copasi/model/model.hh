#ifndef DUNE_COPASI_MODEL_MODEL_HH
#define DUNE_COPASI_MODEL_MODEL_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/operator/operator.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/grid/concepts/grid.hh>
#include <dune/grid/concepts/gridview.hh>

#include <dune/common/parametertree.hh>

#include <any>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>

namespace Dune::Copasi {

template<Dune::Concept::Grid Grid_,
         Dune::Concept::GridView GridView_ = typename Grid_::LeafGridView,
         class RangeQuatinty_ = double,
         class TimeQuantity_ = double>
struct Model
{

  using Grid = Grid_;
  using GridView = GridView_;
  using TimeQuantity = TimeQuantity_;
  using RangeQuatinty = RangeQuatinty_;

  struct State
  {
    std::shared_ptr<const Grid> grid;
    TimeQuantity time{ 0. };
    std::any coefficients = nullptr;
    std::any basis = nullptr;
  };

  using GridFunction = Dune::Functions::GridViewFunction<
    RangeQuatinty(typename GridView::template Codim<0>::Geometry::GlobalCoordinate),
    GridView>;

  Model() = default;
  Model(const Model&) = delete;
  Model(Model&&) = delete;

  Model& operator=(const Model&) = delete;
  Model& operator=(Model&&) = delete;

  virtual ~Model() = default;

  [[nodiscard]] virtual std::unique_ptr<State> make_state(const std::shared_ptr<const Grid>&,
                                                          const ParameterTree&) const = 0;

  virtual void interpolate(State&, const std::unordered_map<std::string, GridFunction>&) const = 0;

  [[nodiscard]] virtual std::unordered_map<std::string, GridFunction> make_initial(
    const Grid&,
    const ParameterTree&) const = 0;

  // copies the state into the grid function
  [[nodiscard]] GridFunction make_compartment_function(const State& state,
                                                       std::string_view name) const
  {
    return make_compartment_function(std::make_shared<State>(state), name);
  }

  // shares the state with the grid function
  [[nodiscard]] virtual GridFunction make_compartment_function(const std::shared_ptr<const State>&,
                                                               std::string_view) const = 0;

  virtual std::nullptr_t make_membrane_function(const State&, std::string_view) const
  {
    throw format_exception(NotImplemented{}, "Model has no membrane functions");
  }

  virtual std::nullptr_t make_membrane_function(const std::shared_ptr<const State>&,
                                                std::string_view) const
  {
    throw format_exception(NotImplemented{}, "Model has no membrane functions");
  }

  virtual void write_vtk(const State&, const std::filesystem::path&, bool) const
  {
    throw format_exception(NotImplemented{}, "Model write has not been implemented");
  }

  [[nodiscard]] virtual std::unique_ptr<OneStep<State>> make_step_operator(
    const State&,
    const ParameterTree&) const = 0;

  virtual std::map<std::string, double> reduce(const State&, const ParameterTree&) const
  {
    throw format_exception(NotImplemented{}, "Model reduce has not been implemented");
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MODEL_HH
