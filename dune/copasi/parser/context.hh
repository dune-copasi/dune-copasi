#ifndef DUNE_COPASI_PARSER_CONTEXT_HH
#define DUNE_COPASI_PARSER_CONTEXT_HH

#include <dune/copasi/parser/factory.hh>

#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>

#include <unordered_map>
#include <vector>

namespace Dune::Copasi {

class ParserContext
{

public:
  explicit ParserContext(const ParameterTree& config = {});

  void add_context(Parser& parser) const;

  typename Parser::Function0D make_function(ParserType, const std::array<std::string,0>&, std::string_view) const;
  typename Parser::Function1D make_function(ParserType, const std::array<std::string,1>&, std::string_view) const;
  typename Parser::Function2D make_function(ParserType, const std::array<std::string,2>&, std::string_view) const;
  typename Parser::Function3D make_function(ParserType, const std::array<std::string,3>&, std::string_view) const;
  typename Parser::Function4D make_function(ParserType, const std::array<std::string,4>&, std::string_view) const;

  // parse functions of the form `arg0, arg1, ...: expr`. Returns a vector of arguments and a string
  // containing the expression
  static std::tuple<std::vector<std::string_view>, std::string_view>
  parse_function_expression(std::string_view fnc_expr);

private:

  void add_independent_context(Parser& parser) const;

  ParserType _default_context_parser;
  std::unordered_map<std::string, double> _constants;
  std::unordered_map<std::string, std::function<typename Parser::Function1D()>> _func_maker_1;
  std::unordered_map<std::string, std::function<typename Parser::Function2D()>> _func_maker_2;
  std::unordered_map<std::string, std::function<typename Parser::Function3D()>> _func_maker_3;
  std::unordered_map<std::string, std::function<typename Parser::Function4D()>> _func_maker_4;
  std::unordered_map<std::string, ParameterTree> _functions_expr;
  std::unordered_map<std::string, std::function<double(std::size_t)>> _cell_data;
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_PARSER_CONTEXT_HH
