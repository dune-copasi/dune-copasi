#ifndef DUNE_COPASI_PARSER_FACTORY_HH
#define DUNE_COPASI_PARSER_FACTORY_HH

#include <dune/copasi/parser/parser.hh>

#include <memory>
#include <unordered_map>

namespace Dune::Copasi {

enum class ParserType
{
  None,
  SymEngine,
  SymEngineSBML,
  MuParser,
  ExprTk
};

extern const ParserType default_parser;

inline const static std::unordered_map<ParserType, std::string_view> parser2string = {
  { ParserType::None, "None" },
  { ParserType::SymEngine, "SymEngine" },
  { ParserType::SymEngineSBML, "SymEngineSBML" },
  { ParserType::MuParser, "MuParser" },
  { ParserType::ExprTk, "ExprTk" }
};

inline const static std::unordered_map<std::string_view, ParserType> string2parser = {
  { "None", ParserType::None },
  { "SymEngine", ParserType::SymEngine },
  { "SymEngineSBML", ParserType::SymEngineSBML },
  { "MuParser", ParserType::MuParser },
  { "ExprTk", ParserType::ExprTk }
};

inline const static std::string default_parser_str =
  std::string{ parser2string.at(default_parser) };

[[nodiscard]] std::unique_ptr<Parser>
make_parser(ParserType parser_type = default_parser);

[[nodiscard]] std::unique_ptr<Parser>
make_parser(std::string_view parser_type);

} // namespace Dune::Copasi

#endif // DUNE_COPASI_PARSER_FACTORY_HH
