#ifndef DUNE_COPASI_COMMON_GRID_FUNCTION_EXPRESSION_ADAPTER_HH
#define DUNE_COPASI_COMMON_GRID_FUNCTION_EXPRESSION_ADAPTER_HH

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
 * @tparam     GV     Grid View
 * @tparam     Time   Time Field
 */
template<class GV, class Time>
class ParserToGridFunctionAdapter
  // : public PDELab::GridFunctionBase<
  //     PDELab::GridFunctionTraits<GV, RF, 1, FieldVector<RF, 1>>,
  //     ParserToGridFunctionAdapter<GV, RF>>
{
  using Position = typename GV::template Codim<0>::Geometry::GlobalCoordinate;
public:
  // using Traits = PDELab::GridFunctionTraits<GV, RangeField, 1, FieldVector<RangeField, 1>>;

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
  {
    assert(not _parser->compiled());
    _parser->define_constant("dim", GV::dimension);
#if HAVE_UNITS
    // _parser->define_variable("t", &(_time.number()));
    // _parser->define_variable("x", &(_global_pos[0].number()));
    // if (GV::dimension >= 2) _parser->define_variable("y", &(_global_pos[1].number()));
    // if (GV::dimension >= 3) _parser->define_variable("z", &(_global_pos[2].number()));
#else
    _parser->define_variable("t", &_time);
    _parser->define_variable("x", &_global_pos[0]);
    if (GV::dimension >= 2) _parser->define_variable("y", &_global_pos[1]);
    if (GV::dimension >= 3) _parser->define_variable("z", &_global_pos[2]);
#endif
  }

  ParserToGridFunctionAdapter(const ParserToGridFunctionAdapter&) = delete;
  ParserToGridFunctionAdapter(ParserToGridFunctionAdapter&&) = delete;

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
  template<class Range>
  void evaluate(const auto& e, const auto& pos_local, Range& result) const
  {
    assert(_parser->compiled());
#if HAVE_UNITS
    DUNE_THROW(NotImplemented, "");
#else
    _global_pos = e.geometry().global(pos_local);
#endif
    result = Range{_parser->eval()};
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
  void setTime(Time t) { _time = t; }

private:

  Logging::Logger _logger;
  GV _gv;
  std::shared_ptr<Parser> _parser;
  Time _time;
  mutable Position _global_pos;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_COMMON_GRID_FUNCTION_EXPRESSION_ADAPTER_HH
