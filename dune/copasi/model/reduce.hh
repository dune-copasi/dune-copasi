#ifndef DUNE_COPASI_MODEL_REDUCE_HH
#define DUNE_COPASI_MODEL_REDUCE_HH

#include <dune/copasi/model/local_equations/local_equations.hh>
#include <dune/copasi/model/local_equations/functor_factory_parser.hh>
#include <dune/copasi/parser/context.hh>
#include <dune/copasi/finite_element/local_basis_cache.hh>

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/common/local_container.hh>
#include <dune/pdelab/common/trace.hh>

#include <dune/grid/common/partitionset.hh>

#include <dune/common/float_cmp.hh>

#include <function2/function2.hpp>
#include <spdlog/spdlog.h>

#include <fmt/color.h>
#include <fmt/core.h>

#include <string>
#include <map>

namespace Dune::Copasi {




// Evaluates a evaluation/reduction/transformation algorigthm over the grid.
// Starting with val = 0, the this function evaluates `val = reduction(val, evaluation(), weight)` for each quadrature point of the grid. Finally, the result gets transformed by `value = transformation(value)`
// For each sub-tree in the parameter tree, the evaluation and reduce expressions are extracted to form the evaluation/reduce operation for its key.
template<PDELab::Concept::Basis Basis>
  requires Concept::LocalBasisTree<typename Basis::LocalView::Tree>
inline static std::map<std::string, double>
reduce(const Basis& basis,
                const auto& coefficients,
                auto time,
                const ParameterTree& config,
                std::shared_ptr<const FunctorFactory<Basis::EntitySet::GridView::dimension>> functor_factory = nullptr)
{
  constexpr std::size_t dim = Basis::EntitySet::GridView::dimension;
  TRACE_EVENT("dune", "Reduce");

  if (basis.entitySet().size(0) == 0) {
    return {};
  }

  std::shared_ptr<const ParserContext> parser_context;
  auto functor_factory_parser = std::dynamic_pointer_cast<const FunctorFactoryParser<dim>>(functor_factory);
  if (functor_factory_parser) {
    parser_context = functor_factory_parser->parser_context();
  }
  if (not parser_context) {
    parser_context = std::make_shared<const ParserContext>();
  }

  auto parser_type = string2parser.at(config.get("parser_type", default_parser_str));

  auto lbasis = basis.localView();
  PDELab::LocalContainerBuffer lcoeff{basis, coefficients};

  auto first_compartment_finite_elment = [](const Concept::CompartmentLocalBasisNode auto& lnode) noexcept {
    for (std::size_t i = 0; i != lnode.degree(); ++i)
      if (lnode.child(i).size() != 0)
        return lnode.child(i).finiteElement();
    std::terminate();
  };

  auto first_finite_element = overload(first_compartment_finite_elment,
  [first_compartment_finite_elment](const Concept::MultiCompartmentLocalBasisNode auto& lnode) noexcept {
    for (std::size_t i = 0; i != lnode.degree(); ++i)
      if (lnode.child(i).size() != 0)
        return first_compartment_finite_elment(lnode.child(i));
    std::terminate();
  });

  using FEM = std::decay_t<decltype(first_finite_element(lbasis.tree()))>;

  auto leqs = LocalEquations<dim>::make(lbasis);
  leqs->time = time;
  LocalBasisCache<typename FEM::Traits::LocalBasisType::Traits> fe_cache;
  std::vector<FieldVector<double, dim>> gradphi;

  std::map<std::string, double> values;
  double factor = 0.;
  std::vector<fu2::unique_function<void() const>> evals;
  // configure parsers of the evaluation/reduction passes
  for (auto key : config.getSubKeys()) {
    if (config.hasKey(key + ".evaluation.expression")) {
      auto evaluation = functor_factory->make_scalar(key + ".evaluation", config.sub(key).sub("evaluation"), *leqs);
      if (not evaluation)
        continue;
      values[key] = 0.;
      fu2::unique_function<double(double, double, double) const> reduction = [](auto init, auto val, auto weight) { return init + val*weight; };
      if (config.hasKey(key + ".reduction.expression")) {
        auto [args, expr] = parser_context->parse_function_expression(config.sub(key).sub("reduction")["expression"]);
        if (args.size() != 3)
          throw format_exception(IOError{}, "Reduction arguments must be exactly 3");
        auto sub_parser_type = string2parser.at(config.sub(key).get("reduction.parser_type", std::string{parser2string.at(parser_type)}));
        auto str_args = std::array{std::string{args[0]}, std::string{args[1]}, std::string{args[2]}};
        reduction = parser_context->make_function(sub_parser_type, str_args, expr);
      }
      // store actual reduction operation in a functor
      evals.emplace_back([_val = std::ref(values[key]), _factor = std::ref(factor), _evaluation = std::move(evaluation), _reduction = std::move(reduction)](){
        _val.get() = _reduction(_val, _evaluation(), _factor);
      });
    }
  }

  if (values.empty()) {
    return {};
  }

  spdlog::info("Evaluating reduce operators");

  for (const auto& cell : elements(basis.entitySet(), Partitions::interior)) {
    auto geo = cell.geometry();

    lbasis.bind(cell);
    lcoeff.load(lbasis, std::false_type{});
    const auto& finite_element = first_finite_element(lbasis.tree());
    fe_cache.bind(finite_element);
    const auto& rule = fe_cache.rule();

    assert(geo.affine());
    const auto& jac = geo.jacobianInverse(rule[0].position());

    for (std::size_t q = 0; q != rule.size(); ++q) {
      const auto [position, weight] = rule[q];
      leqs->position = geo.global(position);

      const auto& phi = fe_cache.evaluateFunction(q);
      const auto& jacphi = fe_cache.evaluateJacobian(q);

      gradphi.resize(finite_element.size());
      for (std::size_t dof = 0; dof != gradphi.size(); ++dof)
        gradphi[dof] = (jacphi[dof] * jac)[0];

      // evaluate values at quadrature point
      forEachLeafNode(lbasis.tree(), [&](const auto& node) {
        if (node.size() == 0)
          return;
        auto& value = leqs->get_value(node);
        auto& gradient = leqs->get_gradient(node);
        value = 0.;
        gradient = 0.;
        for (std::size_t dof = 0; dof != node.size(); ++dof) {
          value += phi[dof] * lcoeff(node, dof);
          gradient += gradphi[dof] * lcoeff(node, dof);
        }
      });

      // evaluate functions at the current quadrature point
      factor = weight * geo.integrationElement(position);
      for (const auto& eval : evals)
        eval();
    }

    lbasis.unbind();
    leqs->clear();
  }

  if (basis.entitySet().grid().comm().size() > 1){
    throw format_exception(ParallelError{}, "Transoform and Reduce is not implemented in parallel");
  }

  // report results
  std::size_t max_key_chars = 0;
  std::size_t max_val_chars = 5;
  for (auto [key, value] : values)
    max_key_chars = std::max(max_key_chars, key.size());

  spdlog::info("   ┌{0:─^{1}}┐", "", max_key_chars+max_val_chars+13);
  bool do_throw = false;
  for (auto& [key, value] : values) {

    if (config.sub(key).hasKey("transformation.expression")) {
      auto [args, expr] = parser_context->parse_function_expression(config.sub(key)["transformation.expression"]);
      if (args.size() != 1)
        throw format_exception(IOError{}, "Warning function must have exactly 1 argument");
      auto sub_parser_type = string2parser.at(config.sub(key).get("transformation.parser_type", std::string{parser2string.at(parser_type)}));
      auto transformation = parser_context->make_function(sub_parser_type, std::array{std::string{args[0]}}, expr);
      value = transformation(value);
    }

    auto sty_key = fmt::styled(key, fmt::emphasis::bold);
    auto verbose = config.sub(key).get("verbose", 1);
    std::function<void()> report_value = [&](){
      if (verbose) {
        spdlog::info("   | {0:>{2}} := {1: .{3}e} |", sty_key, value, max_key_chars, max_val_chars);
      }
    };

    if (config.sub(key).hasKey("warn.expression")) {
      auto [args, expr] = parser_context->parse_function_expression(config.sub(key)["warn.expression"]);
      if (args.size() != 1)
        throw format_exception(IOError{}, "Warning function must have exactly 1 argument");
      auto sub_parser_type = string2parser.at(config.sub(key).get("warn.parser_type", std::string{parser2string.at(parser_type)}));
      auto warn = parser_context->make_function(sub_parser_type, std::array{std::string{args[0]}}, expr);
      if (FloatCmp::ne(warn(value), 0.)) {
        report_value = [&]() {
          spdlog::warn("| {0:>{2}} := {1: .{3}e} |", sty_key, value, max_key_chars, max_val_chars);
        };
      }
    }

    if (config.sub(key).hasKey("error.expression")) {
      auto [args, expr] = parser_context->parse_function_expression(config.sub(key)["error.expression"]);
      if (args.size() != 1)
        throw format_exception(IOError{}, "Error function must have exactly 1 argument");
      auto sub_parser_type = string2parser.at(config.sub(key).get("error.parser_type", std::string{parser2string.at(parser_type)}));
      auto error = parser_context->make_function(sub_parser_type, std::array{std::string{args[0]}}, expr);
      if (FloatCmp::ne(error(value), 0.)) {
        do_throw = true;
        report_value = [&]() {
          spdlog::error("  | {0:>{2}} := {1: .{3}e} |", sty_key, value, max_key_chars, max_val_chars);
        };
      }
    }

    report_value();
  }
  spdlog::info("   └{0:─^{1}}┘", "", max_key_chars+max_val_chars+13);

  if (do_throw) {
    throw Exception{};
  }

  return values;

}

} // namespace Dune::Copasi

#endif // DUNE_COPASI_MODEL_REDUCE_HH