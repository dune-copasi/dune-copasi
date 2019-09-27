#ifndef DUNE_MUPARSER_DATA_HANDLER_HH
#define DUNE_MUPARSER_DATA_HANDLER_HH

#include <dune/common/hybridutilities.hh>
#include <dune/common/parametertree.hh>

#include <muParser.h>

#include <memory>
#include <string>
#include <vector>

namespace Dune::Copasi {

// This is a trick to make happy muparser about static functions

template<class T>
struct MuParserDataHandler
{
  MuParserDataHandler() {}

  ~MuParserDataHandler()
  {
    MuParserDataHandler<T>::_names.clear();
    MuParserDataHandler<T>::_functions.clear();
  }

  template<int i>
  static double function_wrapper(double x, double y)
  {
    assert(MuParserDataHandler<T>::_functions[i]);
    return (*MuParserDataHandler<T>::_functions[i])(x, y);
  }

  void add_tiff_functions(const ParameterTree& data_config)
  {
    const auto& keys = data_config.getValueKeys();
    for (std::size_t i = 0; i < keys.size(); i++) {
      _names.push_back(keys[i]);
      std::string filename = data_config[keys[i]];
      _functions.push_back(std::make_shared<T>(filename));
    }
  }

  template<std::size_t max_data = 100>
  void set_functions(mu::Parser& parser)
  {
    auto indices = Dune::range(std::integral_constant<std::size_t, max_data>{});
    Dune::Hybrid::forEach(indices, [&](auto i) {
      if (i < _functions.size())
        parser.DefineFun(_names[i], function_wrapper<i>);
    });
  }

  static std::vector<std::string> _names;
  static std::vector<std::shared_ptr<T>> _functions;
};

template<class T>
std::vector<std::string> MuParserDataHandler<T>::_names = {};

template<class T>
std::vector<std::shared_ptr<T>> MuParserDataHandler<T>::_functions = {};

} // namespace Dune::Copasi

#endif // DUNE_MUPARSER_DATA_HANDLER_HH