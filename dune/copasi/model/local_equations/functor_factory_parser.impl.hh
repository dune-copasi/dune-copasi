#ifndef DUNE_COPASI_MODEL_LOCAL_VALUES_PARSED_FUNCTOR_FACTORY_IMPL_HH
#define DUNE_COPASI_MODEL_LOCAL_VALUES_PARSED_FUNCTOR_FACTORY_IMPL_HH

#include <dune-copasi-config.h>

#include <dune/copasi/model/local_equations/functor_factory_parser.hh>
#include <dune/copasi/model/local_equations/local_equations.hh>

#include <regex>

#ifdef DUNE_COPASI_PRECOMPILED_MODE
#warning "Including this file in pre-compiled mode may defeat the purpose of pre-compilation"
#endif

namespace Dune::Copasi {

namespace Impl {
// TODO: match scientific notation
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
                                       bool is_membrane_expression) const -> ScalarFunctor
{
  return parse_scalar_expression(config, local_values, is_membrane_expression);
}

template<std::size_t dim>
auto
FunctorFactoryParser<dim>::make_vector(std::string_view /*prefix*/,
                                       const ParameterTree& config,
                                       const LocalDomain<dim>& local_values,
                                       bool is_membrane_expression) const -> VectorFunctor
{
  // create one parser for each entry of the vector
  std::array<ScalarFunctor, dim> vector_parser;
  bool is_active = false;
  std::vector<std::string> dim_name = { "x", "y", "z" };
  for (std::size_t i = 0; i != dim; ++i) {
    is_active |= bool(vector_parser[i] = parse_scalar_expression(
                        config.sub(dim_name.at(i), true), local_values, is_membrane_expression));
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
                                             bool is_membrane_expression) const
  -> TensorApplyFunctor
{
  // diffusion apply parser
  std::vector<std::string> dim_name = { "x", "y", "z" };
  std::string type = config.get("type", "scalar");
  if (type == "scalar") {
    if (auto parser = parse_scalar_expression(config, local_values, is_membrane_expression)) {
      return [parser = std::move(parser)](Vector vec)
               DUNE_COPASI_FUNCTOR_NOEXCEPT { return parser()[0] * vec; };
    }
    return nullptr;
  } else if (type == "tensor") {
    // create one parser for each entry of the tensor
    std::array<std::array<ScalarFunctor, dim>, dim> tensor_parser;
    bool is_active = false;
    for (std::size_t i = 0; i != dim; ++i) {
      for (std::size_t j = 0; j != dim; ++j) {
        is_active |= bool(tensor_parser[i][j] = parse_scalar_expression(
                            config.sub(dim_name.at(i) + dim_name.at(j), true),
                            local_values,
                            is_membrane_expression));
      }
    }
    if (not is_active) {
      return nullptr;
    }

    // move parsers into a lambda that evaluates matrix-vector product
    return { [_tensor_parser = std::move(tensor_parser)](Vector in) DUNE_COPASI_FUNCTOR_NOEXCEPT {
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
                                                   bool is_membrane_expression) const
  -> ScalarFunctor
{
  const auto& expression = config["expression"];
  if (expression.empty() or std::regex_match(expression, Impl::zero_regex)) {
    return nullptr;
  }
  if (std::regex_match(expression, Impl::float_regex)) {
    double value = std::stod(expression);
    return [_value = value] DUNE_COPASI_FUNCTOR_NOEXCEPT { return Scalar{ _value }; };
  } else {
    auto parser_type =
      string2parser.at(config.get("parser_type", std::string{ parser2string.at(_parser_type) }));
    std::vector<std::string> dim_name = { "x", "y", "z" };
    auto parser_ptr = make_parser(parser_type);
    parser_ptr->define_variable("time", &(local_values.time));
    parser_ptr->define_variable("entity_volume", &(local_values.entity_volume));
    if (not is_membrane_expression)
      parser_ptr->define_variable("cell_index", &(local_values.cell_index));
    for (std::size_t i = 0; i != 3; ++i) {
      auto pos_arg = fmt::format("position_{}", dim_name.at(i));
      auto norm_arg = fmt::format("normal_{}", dim_name.at(i));
      if (i < dim) {
        parser_ptr->define_variable(pos_arg, &(local_values.position)[i]);
        if (is_membrane_expression)
          parser_ptr->define_variable(norm_arg, &(local_values.normal)[i]);
      } else {
        parser_ptr->define_constant(pos_arg, 0.);
        if (is_membrane_expression)
          parser_ptr->define_constant(norm_arg, 0.);
      }
    }

    LocalEquations<dim> const* d = dynamic_cast<LocalEquations<dim> const*>(&local_values);
    if (d != nullptr)
      PDELab::forEach(d->nodes(), [&](auto& compartments) {
        for (auto& compartment_fncs : compartments)
          for (auto& component_fncs : compartment_fncs) {
            parser_ptr->define_variable(component_fncs.name, &(component_fncs.value[0]));
            for (std::size_t i = 0; i != 3; ++i) {
              auto grad_arg = fmt::format("grad_{}_{}", component_fncs.name, dim_name.at(i));
              if (i < dim) {
                parser_ptr->define_variable(grad_arg, &(component_fncs.gradient)[i]);
              } else {
                parser_ptr->define_constant(grad_arg, 0.);
              }
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
    return [_parser_ptr = std::move(parser_ptr)] noexcept {
      return Scalar{ std::invoke(*_parser_ptr) };
    };
  }
}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_LOCAL_VALUES_PARSED_FUNCTOR_FACTORY_IMPL_HH
