#ifndef DUNE_COPASI_PARSER_FACTORY_HH
#define DUNE_COPASI_PARSER_FACTORY_HH

#include <dune/copasi/parser/parser.hh>

#if HAVE_MUPARSER
#include <dune/copasi/parser/mu.hh>
#endif // HAVE_MUPARSER

#if HAVE_EXPRTK
#include <dune/copasi/parser/exprtk.hh>
#endif // HAVE_EXPRTK

#if HAVE_SYMENGINE
#include <dune/copasi/parser/symengine.hh>
#endif // HAVE_SYMENGINE

#include <dune/common/exceptions.hh>

#include <memory>

namespace Dune::Copasi {

enum class ParserType
{
  None,
  MuParser,
  SymEngine,
  SymEngineSBML,
  ExprTk
};

static ParserType default_parser =
#if HAVE_MUPARSER
  ParserType::MuParser;
#elif HAVE_SYMENGINE
  ParserType::SymEngineSBML;
#elif HAVE_EXPRTK
  ParserType::ExprTk;
#else
#warning "No known parser was found"
  ParserType::None;
#endif

inline std::unique_ptr<Parser>
make_parser(ParserType parser_type = default_parser)
{
  if (parser_type == ParserType::MuParser) {
#if HAVE_MUPARSER
    return std::make_unique<MuParser>();
#else
    DUNE_THROW(Exception, "MuParser is not available");
#endif // HAVE_MUPARSER
  }

  if (parser_type == ParserType::SymEngine) {
#if HAVE_SYMENGINE
    return std::make_unique<SymEngineParser>(SymEngineParser::Type::Native);
#else
    DUNE_THROW(Exception, "ExprTkParser is not available");
#endif // HAVE_SYMENGINE
  }

  if (parser_type == ParserType::SymEngineSBML) {
#if HAVE_SYMENGINE
    return std::make_unique<SymEngineParser>(SymEngineParser::Type::SBML);
#else
    DUNE_THROW(Exception, "ExprTkParser is not available");
#endif // HAVE_SYMENGINE
  }

  if (parser_type == ParserType::ExprTk) {
#if HAVE_EXPRTK
    return std::make_unique<ExprTkParser>();
#else
    DUNE_THROW(Exception, "ExprTkParser is not available");
#endif // HAVE_EXPRTK
  }

  assert(parser_type == ParserType::None);
  DUNE_THROW(Exception, "No parser type selected");
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_PARSER_FACTORY_HH
