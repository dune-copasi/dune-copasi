
#include <dune/copasi/parser/context.hh>

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/tiff_grayscale.hh>
#include <dune/copasi/parser/factory.hh>
#include <dune/copasi/parser/parser.hh>

#include <dune/pdelab/common/trace.hh>

#if __has_include(<parafields/randomfield.hh>)
#include <parafields/randomfield.hh>
#include <dune/grid/yaspgrid.hh>
#endif

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/parametertree.hh>

#include <spdlog/spdlog.h>

#include <fmt/core.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
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
};

ParserContext::ParserContext(const ParameterTree& config)
{
  _default_context_parser = string2parser.at(config.get("parser_type", default_parser_str));
  for (const auto& sub : config.getSubKeys()) {
    const auto& sub_config = config.sub(sub, true);
    const std::string type = sub_config["type"];
    if (type == "constant") {
      _constants[sub] = sub_config.template get<double>("value");
    } else if (type == "function") {
      _functions_expr[sub] = sub_config;
    } else if (type == "tiff") {
      auto tiff = std::make_shared<TIFFGrayscale>(sub_config["path"]);
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
#if __has_include(<parafields/randomfield.hh>)
      spdlog::info("Generating '{}' random field", sub);
      auto rmf_dim = sub_config.template get<std::vector<double>>("grid.extensions").size();
      auto seed = sub_config.template get<unsigned int>("seed", 0);
      auto file = sub_config.get("writer.vtk.path", "");
      Hybrid::switchCases(std::index_sequence<1, 2, 3>{}, rmf_dim, [&](auto dim) {
        auto field = std::make_shared<parafields::RandomField<RandomFieldTraits<dim>>>(sub_config);
        {
          TRACE_EVENT("dune", "Parafield::Generate");
          if (seed == 0) {
            field->generate();
          } else {
            field->generate(seed);
          }
        }

        auto upper_right = sub_config.get("grid.extensions", FieldVector<double, dim>(1.));
        if (not file.empty()) {
          TRACE_EVENT("dune", "Parafield::Write");
          // make grid for writing it into vtk file
          std::array<int, dim> cells{};
          std::fill_n(std::begin(cells), dim, 1);
          cells = sub_config.get("grid.cells", cells);
          auto levels = config.get<int>("grid.refinement_level", 1);
          YaspGrid<dim> yasp_grid(upper_right, cells);
          yasp_grid.globalRefine(levels);
          // TODO(sospinar): change for dune-vtk which writes structured grids much faster
          spdlog::info("Writing random field");
          field->writeToVTK(file, yasp_grid.leafGridView());
        }

        if constexpr (dim == 1) {
          _functions_1[sub] = [upper_right, _field = std::move(field)](auto pos_x) {
            pos_x = std::clamp(0., pos_x, upper_right[0]);
            FieldVector<double, 1> res;
            _field->evaluate(FieldVector<double, 1>{ pos_x }, res);
            return res[0];
          };
        }
        if constexpr (dim == 2) {
          _functions_2[sub] = [upper_right, _field = std::move(field)](auto pos_x, auto pos_y) {
            pos_x = std::clamp(0., pos_x, upper_right[0]);
            pos_y = std::clamp(0., pos_y, upper_right[1]);
            FieldVector<double, 1> res;
            _field->evaluate(FieldVector<double, 2>{ pos_x, pos_y }, res);
            return res[0];
          };
        }
        if constexpr (dim == 3) {
          _functions_3[sub] = [upper_right, _field = std::move(field)](auto pos_x, auto pos_y, auto pos_z) {
            pos_x = std::clamp(0., pos_x, upper_right[0]);
            pos_y = std::clamp(0., pos_y, upper_right[1]);
            pos_z = std::clamp(0., pos_z, upper_right[2]);
            FieldVector<double, 1> res;
            _field->evaluate(FieldVector<double, 3>{ pos_x, pos_y, pos_z }, res);
            return res[0];
          };
        }
      });
#else
      throw format_exception(NotImplemented{}, "Library was compiled without random fields support");
#endif
    } else if (type == "cell_data") {
      throw format_exception(NotImplemented{}, "cell data is not implemented");
    } else {
      throw format_exception(IOError{}, "Unknown type {}", type);
    }
  }
}

template<std::size_t arg_size>
auto make_callable(ParserType parser_type, const std::array<std::string,arg_size>& args, std::string_view expression, auto add_independent_context) {
  std::shared_ptr sub_parser = make_parser(parser_type);

  // setup storage of arguments and define variables in the parser
  std::shared_ptr<double[]> arg_vals(new double[args.size()]);
  for (std::size_t i = 0; i != arg_size; ++i) {
    sub_parser->define_variable(std::string{ args[i] }, &arg_vals[i]);
  }

  sub_parser->set_expression(std::string{ expression });
  add_independent_context(*sub_parser);
  sub_parser->compile();
  return [_sub_parser = std::move(sub_parser), _arg_vals = std::move(arg_vals)](auto... args) {
    auto arg_tuple = std::make_tuple(args...);
    Hybrid::forEach(std::make_index_sequence<sizeof...(args)>{}, [&](auto i){
      _arg_vals[i] = std::get<i>(arg_tuple);
    });
    return std::invoke(*_sub_parser);
  };
}


typename Parser::Function0D ParserContext::make_function(ParserType parser_type, const std::array<std::string,0>& args, std::string_view expression) const {
  return make_callable(parser_type, args, expression, [&](auto& parser){add_independent_context(parser);});
}

typename Parser::Function1D ParserContext::make_function(ParserType parser_type, const std::array<std::string,1>& args, std::string_view expression) const {
  return make_callable(parser_type, args, expression, [&](auto& parser){add_independent_context(parser);});
}

typename Parser::Function2D ParserContext::make_function(ParserType parser_type, const std::array<std::string,2>& args, std::string_view expression) const {
  return make_callable(parser_type, args, expression, [&](auto& parser){add_independent_context(parser);});
}

typename Parser::Function3D ParserContext::make_function(ParserType parser_type, const std::array<std::string,3>& args, std::string_view expression) const {
  return make_callable(parser_type, args, expression, [&](auto& parser){add_independent_context(parser);});
}

typename Parser::Function4D ParserContext::make_function(ParserType parser_type, const std::array<std::string,4>& args, std::string_view expression) const {
  return make_callable(parser_type, args, expression, [&](auto& parser){add_independent_context(parser);});
}

void
ParserContext::add_context(Parser& parser) const
{
  add_independent_context(parser);

  // add functions to the parser
  for (const auto& name_config_pair : _functions_expr) {
    const std::string &name{name_config_pair.first};
    const ParameterTree& config{name_config_pair.second};
    auto parser_type = _default_context_parser;
    if (config.hasKey("parser_type"))
      parser_type = string2parser.at(config["parser_type"]);
    std::vector<std::string_view> args;
    std::string_view expr;
    std::tie(args, expr) = parse_function_expression(config["expression"]);
    Hybrid::switchCases(std::make_index_sequence<5>{}, args.size(), [&](auto dim){
      std::array<std::string,dim> array_args;
      std::move(args.begin(), args.end(), array_args.begin());
      auto fnc = make_function(parser_type, array_args, expr);
      if (config.get("interpolate", false)) {
        // replace actual function with an interpolation
        if constexpr (dim == 1) {
          const auto intervals_u = config.get("interpolation.intervals", 1000u);
          const auto intervals_f = static_cast<double>(intervals_u);
          if (std::trunc(intervals_f/1e5) != 0.0)
            throw format_exception(IOError{}, "Number of interpolation intervals is too big!");
          if (intervals_u < 1u)
            throw format_exception(IOError{}, "At least one interval is required in function {}", name);
          auto domain = config.get("interpolation.domain." + array_args[0],  std::array<double,2>{0., 1.});
          if (domain[0] >= domain[1])
            throw format_exception(IOError{}, "Domain arguments ({}, {}) of function {} are not ordered", domain[0], domain[1], name);
          auto ooo = config.get("interpolation.out_of_bounds", "error");
          std::vector<double> grid_eval(intervals_u+1, 0.);
          double width = (domain[1] - domain[0])/intervals_f;
          // evaluate function over the grid
          for(std::size_t i = 0; i != grid_eval.size(); ++i)
            grid_eval[i] = fnc(domain[0] + static_cast<double>(i)*width);
          auto ctn = intervals_f/(domain[1] - domain[0]);
          // store interpolation instead of function
          auto fnc_raw = [domain, ctn, intervals_u, _grid_eval = std::move(grid_eval)](auto pos) {
            // convert position to one of the intervals
            double interval = (pos - domain[0])*ctn;
            double t = std::modf(interval, &interval);
            auto interval_st = static_cast<unsigned int>(interval);
            auto i = std::max(0u, interval_st);
            auto j = std::min(interval_st, intervals_u+1u);
            return std::lerp(_grid_eval[i], _grid_eval[j], t);
          };
          if (ooo == "clamp") {
            fnc = [domain, _fnc_raw = std::move(fnc_raw)](auto pos) {
              pos = std::clamp(domain[0], pos, domain[1]);
              return _fnc_raw(pos);
            };
          } else if (ooo == "error") {
            fnc = [domain, name, _fnc_raw = std::move(fnc_raw)](auto pos) {
              if ((pos < domain[0]) or (pos > domain[1]))
                throw format_exception(MathError{}, "Interpolation of function {} is out of bounds: {} ∉ [{}, {}]", name, pos, domain[0], domain[1]);
              return _fnc_raw(pos);
            };
          } else {
            throw format_exception(IOError{}, "Not known '{}.interpolation.out_of_bounds = {}'", name, ooo);
          }
        } else {
          throw format_exception(NotImplemented{}, "Cannot interpolate '{}' function with {} arguments", name, dim());
        }
      }
      parser.define_function(name, std::move(fnc));
    }, [&](){
      throw format_exception(NotImplemented{}, "Cannot parse '{}' function with {} arguments", name, args.size());
    });
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
    if (not std::ranges::all_of(arg, [](auto lchar) { return std::isalnum(lchar) or lchar == '_'; })) {
      throw format_exception(
        IOError{},
        "Function argument '{}' must contain only alphanumeric character or underscores.\nExpression: {}",
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
