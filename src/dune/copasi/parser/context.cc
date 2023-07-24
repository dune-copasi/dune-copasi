
#include <dune/copasi/parser/context.hh>

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/tiff_grayscale.hh>
#include <dune/copasi/parser/factory.hh>
#include <dune/copasi/parser/parser.hh>

#include <dune/randomfield/randomfield.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parametertree.hh>

#include <spdlog/spdlog.h>

#include <fmt/core.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cctype>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <ranges>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Dune::Copasi {

template<unsigned int dimension>
class RandomFieldTraits
{
public:
  enum
  {
    dim = dimension
  };
  using RangeField = double;
  using Scalar = FieldVector<double, 1>;
  using DomainField = double;
  using Domain = FieldVector<double, dim>;

  int levels;
  std::vector<DomainField> maxExt;
  std::vector<unsigned int> minCells, maxCells;

  explicit RandomFieldTraits(const ParameterTree& config)
  {
    levels = config.get<int>("grid.refinement_level", 1);
    maxExt = config.get<std::vector<DomainField>>("grid.extensions");
    maxCells = config.get<std::vector<unsigned int>>("grid.cells");

    if (maxExt.size() != maxCells.size()) {
      throw format_exception(IOError{},
                             "'cell[{}]' and 'extension[{}]' vectors differ in size",
                             maxCells.size(),
                             maxExt.size());
    }
    minCells = maxCells;
    for (int i = 0; i < levels - 1; i++) {
      for (unsigned int j = 0; j < maxCells.size(); j++) {
        if (minCells[j] % 2 != 0) {
          throw format_exception(
            Exception{},
            "cannot create enough levels for hierarchical grid, check number of cells");
        }
        minCells[j] /= 2;
      }
    }
  }

  [[nodiscard]] FieldVector<DomainField, dim> L() const
  {
    FieldVector<DomainField, dim> Lvector;

    for (unsigned int i = 0; i < dim; i++) {
      Lvector[i] = maxExt[i];
    }

    return Lvector;
  }

  [[nodiscard]] std::array<int, dim> N() const
  {
    std::array<int, dim> Nvector{};

    for (unsigned int i = 0; i < dim; i++) {
      Nvector[i] = minCells[i];
    }

    return Nvector;
  }

  [[nodiscard]] std::bitset<dim> B() const { return std::bitset<dim>(false); }
};

ParserContext::ParserContext(const ParameterTree& config)
{
  for (const auto& sub : config.getSubKeys()) {
    const auto& sub_config = config.sub(sub, true);
    const std::string type = sub_config["type"];
    if (type == "constant") {
      _constants[sub] = sub_config.template get<double>("value");
    } else if (type == "function") {
      auto parser_type = string2parser.at(sub_config.get("parser_type", default_parser_str));
      _functions_expr[sub] = std::make_pair(parser_type, sub_config["expression"]);
    } else if (type == "tiff") {
      auto tiff = std::make_shared<TIFFGrayscale>(sub_config["file"]);
      _functions_2[sub] = [tiff](const auto& pos_x, const auto& pos_y) {
        return std::invoke(*tiff, pos_x, pos_y);
      };
    } else if (type == "interpolation") {
      auto domain = sub_config.template get<std::vector<double>>("domain");
      auto range = sub_config.template get<std::vector<double>>("range");
      if (not std::is_sorted(domain.begin(), domain.end())) {
        throw format_exception(IOError{}, "The interpolation domain must be sorted");
      }
      if (domain.size() < 2 or domain.size() != range.size()) {
        throw format_exception(
          IOError{},
          "Interpolation range and domain must have at least two points and be the same size");
      }
      _functions_1[sub] = [_domain = std::move(domain),
                           _range = std::move(range)](const auto& pos) {
        auto domain_it = lower_bound(_domain.begin(), _domain.end(), pos);
        if (domain_it == _domain.begin()) {
          return _range.front();
        }
        if (domain_it == _domain.end()) {
          return _range.back();
        }
        auto dist = std::distance(_domain.begin(), domain_it);
        return std::lerp(_range[dist - 1],
                         _range[dist],
                         (pos - _domain[dist - 1]) / (_domain[dist] - _domain[dist - 1]));
      };
    } else if (type == "random_field") {
      spdlog::info("Generating '{}' random field", sub);
      auto rmf_dim = sub_config.template get<std::vector<double>>("grid.extensions").size();
      auto seed = sub_config.template get<unsigned int>("seed", 0);
      auto file = sub_config.get("write_to", "");
      Hybrid::switchCases(std::index_sequence<1, 2, 3>{}, rmf_dim, [&](auto dim) {
        auto field = std::make_shared<RandomField::RandomField<RandomFieldTraits<dim>>>(sub_config);
        if (seed == 0) {
          field->generate();
        } else {
          field->generate(seed);
        }
        const RandomFieldTraits<dim> traits(sub_config);
        const YaspGrid<dim> yaspGrid(traits.L(), traits.N(), traits.B(), 1);
        if (not file.empty()) {
          // TODO(sospinar): change for dune-vtk which writes structured grids much faster
          spdlog::info("Writing random field");
          field->writeToVTK(file, yaspGrid.leafGridView());
        }

        if constexpr (dim == 1) {
          _functions_1[sub] = [_field = std::move(field)](auto pos_x) {
            FieldVector<double, 1> res;
            _field->evaluate(FieldVector<double, 1>{ pos_x }, res);
            return res[0];
          };
        }
        if constexpr (dim == 2) {
          _functions_2[sub] = [_field = std::move(field)](auto pos_x, auto pos_y) {
            FieldVector<double, 1> res;
            _field->evaluate(FieldVector<double, 2>{ pos_x, pos_y }, res);
            return res[0];
          };
        }
        if constexpr (dim == 3) {
          _functions_3[sub] = [_field = std::move(field)](auto pos_x, auto pos_y, auto pos_z) {
            FieldVector<double, 1> res;
            _field->evaluate(FieldVector<double, 3>{ pos_x, pos_y, pos_z }, res);
            return res[0];
          };
        }
      });
    } else if (type == "cell_data") {
      throw format_exception(NotImplemented{}, "cell data is not implemented");
    } else {
      throw format_exception(IOError{}, "Unknown type {}", type);
    }
  }
}

void
ParserContext::add_context(Parser& parser) const
{
  add_independent_context(parser);

  for (const auto& [name, pair] : _functions_expr) {
    const auto& [parser_type, function_expression] = pair;
    auto [args, expr] = parse_function_expression(function_expression);

    std::shared_ptr sub_parser = make_parser(parser_type);

    // setup storage of arguments and define variables in the parser
    std::shared_ptr<double[]> arg_vals(new double[args.size()]);
    for (std::size_t i = 0; i != args.size(); ++i) {
      sub_parser->define_variable(std::string{ args[i] }, &arg_vals[i]);
    }

    sub_parser->set_expression(std::string{ expr });
    add_independent_context(*sub_parser);
    sub_parser->compile();
    if (args.empty()) {
      parser.define_function(
        name, [_sub_parser = std::move(sub_parser)]() { return std::invoke(*_sub_parser); });
    } else if (args.size() == 1) {
      parser.define_function(
        name, [_sub_parser = std::move(sub_parser), _arg_vals = std::move(arg_vals)](auto arg0) {
          _arg_vals[0] = arg0;
          return std::invoke(*_sub_parser);
        });
    } else if (args.size() == 2) {
      parser.define_function(name,
                             [_sub_parser = std::move(sub_parser),
                              _arg_vals = std::move(arg_vals)](auto arg0, auto arg1) {
                               _arg_vals[0] = arg0;
                               _arg_vals[1] = arg1;
                               return std::invoke(*_sub_parser);
                             });
    } else if (args.size() == 3) {
      parser.define_function(name,
                             [_sub_parser = std::move(sub_parser),
                              _arg_vals = std::move(arg_vals)](auto arg0, auto arg1, auto arg2) {
                               _arg_vals[0] = arg0;
                               _arg_vals[1] = arg1;
                               _arg_vals[2] = arg2;
                               return std::invoke(*_sub_parser);
                             });
    } else {
      throw NotImplemented{};
    }
  }
}

// parse functions of the form `arg0, arg1, ...: expr`. returns a vector of arguments and a string
// containing the expression
std::tuple<std::vector<std::string_view>, std::string_view>
ParserContext::parse_function_expression(std::string_view fnc_expr)
{
  // colon splits the argument between arguments and exprssion
  auto colon_pos = fnc_expr.find(':');
  if (colon_pos == std::string_view::npos) {
    throw format_exception(IOError{}, "Function arguments must precede a colon ':'");
  }

  // get expression part
  const std::string_view expr = fnc_expr.substr(colon_pos + 1);

  // find each argument individually
  std::vector<std::string_view> args;
  const auto begin_delimiter = std::numeric_limits<std::size_t>::max();
  std::size_t delimiter = begin_delimiter;

  while (std::max<std::size_t>(delimiter, 0) != colon_pos) {
    // find the end of this argument
    auto arg_count = fnc_expr.substr(delimiter + 1).find_first_of(",:");
    // assert(arg_count != std::string_view::npos && delimiter+arg_count+1 <= colon_pos);

    // extract argument
    std::string_view arg = fnc_expr.substr(delimiter + 1, arg_count);
    // trim whitespace in argument
    if (auto ltrim = arg.find_first_not_of(" \t\r\n"); ltrim != std::string_view::npos) {
      arg = arg.substr(ltrim);
    }
    if (auto rtrim = arg.find_last_not_of(" \t\r\n"); rtrim != std::string_view::npos) {
      arg = arg.substr(0, rtrim + 1);
    }

    // check empty arguments
    if (arg.empty()) {
      if ((delimiter == begin_delimiter) and
          (delimiter + arg_count + 1 == colon_pos)) { // fine: zero arguments
        break;
      }
      // wrong: argument after a comma is empty!
      throw format_exception(IOError{}, "Function arguments shall be named: {}", fnc_expr);
    }

    // check that argument contains valid identifiers
    if (std::isalpha(arg[0]) == 0) {
      throw format_exception(
        IOError{},
        "Function argument '{}' must start with an alphabetic character.\nExpression: {}",
        arg,
        fnc_expr);
    }
    if (not std::ranges::all_of(arg, [](auto lchar) { return std::isalnum(lchar); })) {
      throw format_exception(
        IOError{},
        "Function argument '{}' must contain only alphanumeric character.\nExpression: {}",
        arg,
        fnc_expr);
    }

    // store argument view in vector
    args.emplace_back(arg);

    // advance next argument after the end of the current one
    delimiter += arg_count + 1;
  }

  return make_tuple(std::move(args), expr);
}

void
ParserContext::add_independent_context(Parser& parser) const
{
  for (const auto& [name, value] : _constants) {
    parser.define_constant(name, value);
  }

  for (const auto& [name, func] : _functions_1) {
    parser.define_function(name, func);
  }

  for (const auto& [name, func] : _functions_2) {
    parser.define_function(name, func);
  }

  for (const auto& [name, func] : _functions_3) {
    parser.define_function(name, func);
  }
}

} // namespace Dune::Copasi
