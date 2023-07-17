#include <dune/copasi/parser/factory.hh>
#include <dune/copasi/parser/parser.hh>

#if __has_include(<muParser.h>)
#include <dune/copasi/parser/mu.hh>
#endif

#if __has_include(<exprtk.hpp>)
#include <dune/copasi/parser/exprtk.hh>
#endif

#if __has_include(<symengine_config.h>)
#include <dune/copasi/parser/symengine.hh>
#endif

#include <dune/common/exceptions.hh>

#include <cassert>
#include <memory>
#include <string_view>

namespace Dune::Copasi {

const ParserType default_parser =
#if __has_include(<symengine_config.h>)
  ParserType::SymEngineSBML;
#elif __has_include(<muParser.h>)
  ParserType::MuParser;
#elif __has_include(<exprtk.hpp>)
  ParserType::ExprTk;
#else
#warning "No known parser was found"
  ParserType::None;
#endif

auto
make_parser(ParserType parser_type) -> std::unique_ptr<Parser>
{
  if (parser_type == ParserType::MuParser) {
#if __has_include(<muParser.h>)
    return std::make_unique<MuParser>();
#else
    DUNE_THROW(Exception, "\tMuParser is not available");
#endif
  }

  if (parser_type == ParserType::SymEngine) {
#if __has_include(<symengine_config.h>)
    return std::make_unique<SymEngineParser>(SymEngineParser::Type::Native);
#else
    DUNE_THROW(Exception, "\tExprTkParser is not available");
#endif
  }

  if (parser_type == ParserType::SymEngineSBML) {
#if __has_include(<symengine_config.h>)
    return std::make_unique<SymEngineParser>(SymEngineParser::Type::SBML);
#else
    DUNE_THROW(Exception, "\tExprTkParser is not available");
#endif
  }

  if (parser_type == ParserType::ExprTk) {
#if __has_include(<exprtk.hpp>)
    return std::make_unique<ExprTkParser>();
#else
    DUNE_THROW(Exception, "\tExprTkParser is not available");
#endif
  }

  assert(parser_type == ParserType::None);
  DUNE_THROW(Exception, "\tNo parser type selected");
}

std::unique_ptr<Parser>
make_parser(std::string_view parser_type)
{
  return make_parser(string2parser.at(parser_type));
}

} // namespace Dune::Copasi
