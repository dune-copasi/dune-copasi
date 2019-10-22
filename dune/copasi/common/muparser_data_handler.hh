#ifndef DUNE_MUPARSER_DATA_HANDLER_HH
#define DUNE_MUPARSER_DATA_HANDLER_HH

#include <dune/common/hybridutilities.hh>
#include <dune/common/parametertree.hh>

#include <muParser.h>

#include <memory>
#include <string>
#include <vector>

namespace Dune::Copasi {

/**
 * @brief      MuParser data handler.
 * @details    MuParser only allows to define functions if they are defined
 *             statically. Hence, this is a trick to make happy muparser about
 *             them. We basically store all possible functions in a static
 *             vector, allow mu parser to call them in a static function. The
 *             trick is that we assig the content of each function at run time.
 *             This is mainly for tiff files, but it can be extended for other
 *             cases.
 * @todo       Add static check to check that T is a callable of two arguments
 *
 * @tparam     T     The type to manage statically. It has to have the
 *                   parenthesis operator for two arguments.
 */
template<class T>
struct MuParserDataHandler
{

  /**
   * @brief      Destroys the object.
   */
  ~MuParserDataHandler()
  {
    MuParserDataHandler<T>::_names.clear();
    MuParserDataHandler<T>::_functions.clear();
  }

  /**
   * @brief      The static functions for MuParser
   *
   * @param[in]  x     The x position in global coordinates
   * @param[in]  y     The x position in global coordinates
   *
   * @tparam     i     Index to make each function to be unique
   *
   * @return     The return value of the wrapped function
   */
  template<int i>
  static double function_wrapper(double x, double y)
  {
    assert(MuParserDataHandler<T>::_functions[i]);
    return (*MuParserDataHandler<T>::_functions[i])(x, y);
  }

  /**
   * @brief      Adds tiff functions.
   * @todo       Disable this function if T is not a tiff grayscale
   *
   * @param[in]  data_config  The data configuration
   */
  void add_tiff_functions(const ParameterTree& data_config)
  {
    const auto& keys = data_config.getValueKeys();
    for (std::size_t i = 0; i < keys.size(); i++) {
      _names.push_back(keys[i]);
      std::string filename = data_config[keys[i]];
      _functions.push_back(std::make_shared<T>(filename));
    }
  }

  /**
   * @brief      Sets the functions to muParser.
   *
   * @param      parser         The muParser parser
   *
   * @tparam     max_functions  Maximum set of static functions to set to
   *                            muParser
   */
  template<std::size_t max_functions = 20>
  void set_functions(mu::Parser& parser)
  {
    auto indices =
      Dune::range(std::integral_constant<std::size_t, max_functions>{});
    Dune::Hybrid::forEach(indices, [&](auto i) {
      if (i < _functions.size())
        parser.DefineFun(_names[i], function_wrapper<i>);
    });
  }

  /// vector of strigs descriving the name for each function
  static std::vector<std::string> _names;

  /// vector of callables to be used for each static function
  static std::vector<std::shared_ptr<T>> _functions;
};

/// initialization fo the static vectore of names
template<class T>
std::vector<std::string> MuParserDataHandler<T>::_names = {};

/// initialization of vector of callables
template<class T>
std::vector<std::shared_ptr<T>> MuParserDataHandler<T>::_functions = {};

} // namespace Dune::Copasi

#endif // DUNE_MUPARSER_DATA_HANDLER_HH