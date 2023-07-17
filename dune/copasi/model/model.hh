#ifndef DUNE_COPASI_MODEL_MODEL_HH
#define DUNE_COPASI_MODEL_MODEL_HH

#include <dune-copasi-config.h>
#include <dune/copasi/common/filesystem.hh>

#include <dune/pdelab/operator/operator.hh>

#include <dune/functions/gridfunctions/gridviewfunction.hh>

#include <dune/grid/concepts/grid.hh>
#include <dune/grid/concepts/gridview.hh>

#include <dune/common/parametertree.hh>

#include <any>
#include <memory>
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
    DUNE_THROW(NotImplemented, "\tModel has no membrane functions");
  }

  virtual std::nullptr_t make_membrane_function(const std::shared_ptr<const State>&,
                                                std::string_view) const
  {
    DUNE_THROW(NotImplemented, "\tModel has no membrane functions");
  }

  virtual void write(const State&, const fs::path&, bool) const
  {
    DUNE_THROW(NotImplemented, "\tModel write has not been implemented");
  }

  [[nodiscard]] virtual std::unique_ptr<PDELab::OneStep<State>> make_step_operator(
    const State&,
    const ParameterTree&) const = 0;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_MODEL_HH
