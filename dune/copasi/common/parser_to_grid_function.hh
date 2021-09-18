#ifndef DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH
#define DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH

#include <dune/copasi/parser/factory.hh>
#include <dune/copasi/common/tiff_grayscale.hh>

#include <dune/pdelab/common/function.hh>

#include <dune/logging/logging.hh>

#include <dune/common/dynvector.hh>
#include <dune/common/parametertree.hh>

#include <string>
#include <type_traits>
#include <memory>

namespace Dune::Copasi {

/**
 * @brief      Converts an interface to match an expression to a PDELab grid
 *             function.
 * @details    The resulting grid view is only for scalar expressions
 *
 * @tparam     GV    Grid View
 * @tparam     RF    Range Field
 */
template<class GV, class RF>
class ParserToGridFunctionAdapter
  : public PDELab::GridFunctionBase<
      PDELab::GridFunctionTraits<GV, RF, 1, FieldVector<RF, 1>>,
      ParserToGridFunctionAdapter<GV, RF>>
{
  using RangeField = RF;
public:
  using Traits = PDELab::GridFunctionTraits<GV, RangeField, 1, FieldVector<RangeField, 1>>;

  /**
   * @brief      Constructs a new instance.
   *
   * @param[in]  grid_view          The grid view
   * @param[in]  equation           The math expression
   * @param[in]  do_compile_parser  Bool to compile parser at object
   * construction
   * @param[in]  other_variables    Extra varialbes names to be available in the
   * expression
   */
  ParserToGridFunctionAdapter(const GV& grid_view,
                              const std::shared_ptr<Parser>& parser)
    : _logger(Logging::Logging::componentLogger({}, "model"))
    , _gv{grid_view}
    , _parser{parser}
    , _time{0}
    , _global_pos(0)
  {
    assert(not _parser->compiled());
    _parser->define_variable("t", &_time);
    _parser->define_variable("x", &_global_pos[0]);
    _parser->define_constant("dim", GV::dimension);
    if (GV::dimension >= 2) _parser->define_variable("y", &_global_pos[1]);
    if (GV::dimension >= 3) _parser->define_variable("z", &_global_pos[2]);
  }

public:
  /**
   * @brief      Gets a reference to the grid view
   *
   * @return     The grid view.
   */
  inline const GV& getGridView() const { return _gv; }

  /**
   * @brief      Evaluates extended function on a element
   *
   * @param[in]  e     Entity to operate with
   * @param[in]  x     Local coordinates in the entity
   * @param      y     Resulting value
   *
   * @tparam     E     Entity
   * @tparam     D     Domain
   * @tparam     R     Range
   */
  template<class E, class D, class R>
  void evaluate(const E& e, const D& pos_local, R& result) const
  {
    assert(_parser->compiled());
    _global_pos = e.geometry().global(pos_local);
    result = _parser->eval();
  }

  /**
   * @brief      Get parser
   *
   * @return     Reference to internal parser
   */
  Parser& parser()
  {
    return *_parser;
  }

  /**
   * @brief      Sets the time.
   *
   * @param[in]  t     The new time
   */
  void setTime(const double& t) { _time = t; }

private:

  Logging::Logger _logger;
  GV _gv;
  std::shared_ptr<Parser> _parser;
  RF _time;
  mutable FieldVector<RF,GV::dimension> _global_pos;
};

/**
 * @brief      Gets the muparser expressions.
 *
 * @param[in]  expressions_config  The expressions configuration
 * @param[in]  gf_grid_view        The grid view for the grid function
 * @param[in]  compile             True to compile expression at construction
 *
 * @tparam     GFGridView          Grid view type
 * @tparam     RF                  Range field type
 *
 * @return     Vector with muparser expressions pointers
 */
template<class GFGridView, class RF = double>
auto make_multicomponent_grid_function(
  const ParameterTree& expressions_config,
  const GFGridView& gf_grid_view,
  bool compile = true)
{
  auto vars = expressions_config.getValueKeys();
  sort(begin(vars),end(vars));

  using GridFunction = ParserToGridFunctionAdapter<GFGridView, RF>;
  std::vector<std::shared_ptr<GridFunction>> functions;

  for (std::size_t i = 0; i < vars.size(); ++i) {
    std::shared_ptr parser = make_parser(expressions_config[vars[i]]);
    auto gf = std::make_shared<GridFunction>(gf_grid_view, parser);
    if (compile) parser->compile();
    functions.emplace_back(std::move(gf));
  }

  return PDELab::DynamicPowerGridFunction<GridFunction>{functions};
}

template<class GridFunction>
auto
add_tiff_to_grid_function(GridFunction& grid_function,
                          const ParameterTree& tiff_data,
                          bool compile = true)
{
  // get TIFF data if available and make it available to grid function parses
  for (const auto& data_key : tiff_data.getValueKeys()) {
    auto tiff_func_ptr =
      std::make_shared<TIFFGrayscale<unsigned short>>(tiff_data[data_key]);
    Dune::TypeTree::forEachLeafNode(
      grid_function, [&](auto& leaf_gf, auto tree_path) {

        using LeafGridFunction = std::decay_t<decltype(leaf_gf)>;
        using RangeField = typename LeafGridFunction::Traits::RangeFieldType;

        leaf_gf.parser().define_function(
          tiff_data[data_key], [=](const RangeField& x, const RangeField& y) {
            return std::invoke(*tiff_func_ptr, x, y);
          });
        if (compile)
          leaf_gf.parser().compile();
      });
  }

  if (compile)
    Dune::TypeTree::forEachLeafNode(
      grid_function, [&](auto& leaf_gf, auto tree_path) {
          leaf_gf.parser().compile();
      });

  return grid_function;
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_GRID_FUNCTION_EXPRESSION_ADAPTER_HH
