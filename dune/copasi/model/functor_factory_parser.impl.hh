#ifndef DUNE_COPASI_MODEL_LOCAL_VALUES_PARSED_FUNCTOR_FACTORY_IMPL_HH
#define DUNE_COPASI_MODEL_LOCAL_VALUES_PARSED_FUNCTOR_FACTORY_IMPL_HH

#include <dune/copasi/model/functor_factory_parser.hh>
#include <dune/copasi/model/local_domain.hh>
#include <dune/copasi/common/axis_names.hh>

#include <regex>

#ifdef DUNE_COPASI_PRECOMPILED_MODE
#warning "Including this file in pre-compiled mode may defeat the purpose of pre-compilation"
#endif

namespace Dune::Copasi {

namespace Impl {
static inline const std::regex scientific_notation_regex("/[+\\-]?(?:0|[1-9]\\d*)(?:\\.\\d+)?(?:[eE][+\\-]?\\d+)?/");
static inline const std::regex float_regex("-?([0-9]+)?([\\.]?)([0-9]+)?");
static inline const std::regex zero_regex("-?([0]+)?([\\.]?)([0]+)?");
} // namespace Impl

template<std::size_t dim>
class LocalEquations;

template<std::size_t dim>
auto
FunctorFactoryParser<dim>::make_scalar(std::string_view /*prefix*/,
                                       const ParameterTree& config,
                                       const LocalDomain<dim>& local_values,
                                       int codim) const -> ScalarFunctor
{
  return parse_scalar_expression(config, local_values, codim);
}

template<std::size_t dim>
auto
FunctorFactoryParser<dim>::make_vector(std::string_view /*prefix*/,
                                       const ParameterTree& config,
                                       const LocalDomain<dim>& local_values,
                                       int codim) const -> VectorFunctor
{
  // create one parser for each entry of the vector
  std::array<ScalarFunctor, dim> vector_parser;
  const std::size_t max_axis = std::min<std::size_t>({axis_names.size(), dim});
  bool is_active = false;
  for (std::size_t i = 0; i != max_axis; ++i) {
    is_active |= bool(vector_parser[i] = parse_scalar_expression(
                        config.sub(axis_names[i], true), local_values, codim));
  }
  if (not is_active) {
    return nullptr;
  }

  // move parsers into a lambda that evaluates vector components
  return { [_vector_parser = std::move(vector_parser)]() noexcept {
    Vector vec;
    for (std::size_t i = 0; i != dim; ++i) {
      vec[i] = (_vector_parser[i] ? _vector_parser[i]()[0] : 0.);
    }
    return vec;
  } };
}

template<std::size_t dim>
auto
FunctorFactoryParser<dim>::make_tensor_apply(std::string_view prefix,
                                             const ParameterTree& config,
                                             const LocalDomain<dim>& local_values,
                                             int codim) const
  -> TensorApplyFunctor
{
  // diffusion apply parser
  std::string type = config.get("type", "scalar");
  if (type == "scalar") {
    if (auto parser = parse_scalar_expression(config, local_values, codim)) {
      return [parser = std::move(parser)](Vector vec)
               noexcept { return parser()[0] * vec; };
    }
    return nullptr;
  } else if (type == "tensor") {
    // create one parser for each entry of the tensor
    std::array<std::array<ScalarFunctor, dim>, dim> tensor_parser;
    const std::size_t max_axis = std::min<std::size_t>({axis_names.size(), dim});
    bool is_active = false;
    for (std::size_t i = 0; i != max_axis; ++i) {
      for (std::size_t j = 0; j != max_axis; ++j) {
        if (config.hasSub(axis_names[i] + axis_names[j]))
          is_active |= bool(tensor_parser[i][j] = parse_scalar_expression(
                              config.sub(axis_names[i] + axis_names[j]),
                              local_values,
                              codim));
      }
    }
    if (not is_active) {
      return nullptr;
    }

    // move parsers into a lambda that evaluates matrix-vector product
    return { [_tensor_parser = std::move(tensor_parser)](Vector in) noexcept {
      Vector out;
      for (std::size_t i = 0; i != dim; ++i) {
        out[i] = 0.;
        for (std::size_t j = 0; j != dim; ++j) {
          out[i] += in[j] * (_tensor_parser[i][j] ? _tensor_parser[i][j]()[0] : 0.);
        }
      }
      return out;
    } };
  } else {
    spdlog::error("not known type 'scalar_value.{}.type = {}'", prefix, type);
    std::terminate();
  }
}

template<std::size_t dim>
auto
FunctorFactoryParser<dim>::parse_scalar_expression(const ParameterTree& config,
                                                   const LocalDomain<dim>& local_values,
                                                   int codim) const
  -> ScalarFunctor
{
  auto expression = config.get("expression", std::string{});
  if (expression.empty() or std::regex_match(expression, Impl::zero_regex)) {
    return nullptr;
  }
  if (std::regex_match(expression, Impl::float_regex)
      or std::regex_match(expression, Impl::scientific_notation_regex)) {
    double value = std::stod(expression);
    return [_value = value]() noexcept { return Scalar{ _value }; };
  } else {
    auto parser_type =
      string2parser.at(config.get("parser_type", std::string{ parser2string.at(_parser_type) }));
    auto parser_ptr = make_parser(parser_type);
    parser_ptr->define_variable("time", &(local_values.time));
    parser_ptr->define_variable("entity_volume", &(local_values.entity_volume));
    parser_ptr->define_variable("in_volume", &(local_values.in_volume));
    parser_ptr->define_variable("in_boundary", &(local_values.in_boundary));
    parser_ptr->define_variable("in_skeleton", &(local_values.in_skeleton));
    parser_ptr->define_constant("no_value", std::numeric_limits<double>::max());
    for (std::size_t i = 0; i != axis_names.size(); ++i) {
      auto pos_arg = fmt::format("position_{}", axis_names[i]);
      auto norm_arg = fmt::format("normal_{}", axis_names[i]);
      if (i < dim) {
        parser_ptr->define_variable(pos_arg, &(local_values.position)[i]);
        if (codim == 1)
          parser_ptr->define_variable(norm_arg, &(local_values.normal)[i]);
      } else {
        parser_ptr->define_constant(pos_arg, 0.);
        if (codim == 1)
          parser_ptr->define_constant(norm_arg, 0.);
      }
    }

    // bind cell keys with cell key values
    for (std::size_t j = 0; j < local_values.cell_keys.size(); ++j) {
      parser_ptr->define_variable(local_values.cell_keys[j], &local_values.cell_values[j]);
    }

    local_values.forEachValue(Codim<0>{}, [&](auto name, const auto& value, const auto& gradient){
      parser_ptr->define_variable(std::string{name}, &(value[0]));
      for (std::size_t i = 0; i != axis_names.size(); ++i) {
        auto grad_arg = fmt::format("grad_{}_{}", name, axis_names[i]);
        if (i < dim) {
          parser_ptr->define_variable(grad_arg, &(gradient)[i]);
        } else {
          parser_ptr->define_constant(grad_arg, 0.);
        }
      }
    });

    parser_ptr->set_expression(expression);
    if (_parser_context)
      _parser_context->add_context(*parser_ptr);
    parser_ptr->compile();
    // try to run the parser once, if compilation is wrong,
    // this will throw outside of the functor
    [[maybe_unused]] auto dummy = std::invoke(*parser_ptr);
    return [_parser_ptr = std::move(parser_ptr)]() noexcept {
      return Scalar{ std::invoke(*_parser_ptr) };
    };
  }
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_LOCAL_VALUES_PARSED_FUNCTOR_FACTORY_IMPL_HH
