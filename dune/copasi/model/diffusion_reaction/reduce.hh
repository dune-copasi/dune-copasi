#ifndef DUNE_COPASI_MODEL_DIFFUSION_REACTION_REDUCE_HH
#define DUNE_COPASI_MODEL_DIFFUSION_REACTION_REDUCE_HH

#include <dune/copasi/model/diffusion_reaction/local_equations.hh>
#include <dune/copasi/model/functor_factory_parser.hh>
#include <dune/copasi/parser/context.hh>
#include <dune/copasi/finite_element/local_basis_cache.hh>
#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/fmt_style.hh>

#include <dune/pdelab/concepts/basis.hh>
#include <dune/pdelab/common/local_container.hh>
#include <dune/pdelab/common/trace.hh>
#include <dune/pdelab/common/for_each_entity.hh>

#include <dune/grid/common/partitionset.hh>

#include <dune/common/float_cmp.hh>

#include <function2/function2.hpp>
#include <spdlog/spdlog.h>

#include <fmt/core.h>

#if HAVE_TBB
#include <tbb/enumerable_thread_specific.h>
#endif // HAVE_TBB

#include <string>
#include <map>

namespace Dune::Copasi::DiffusionReaction {

class ReductionError : public Dune::Exception {};

// Evaluates a evaluation/reduction/transformation algorigthm over the grid.
// Starting with val = 0, the this function evaluates `val = reduction(val, evaluation(), weight)` for each quadrature point of the grid. Finally, the result gets transformed by `value = transformation(value)`
// For each sub-tree in the parameter tree, the evaluation and reduce expressions are extracted to form the evaluation/reduce operation for its key.
template<class ExecutionPolicy, PDELab::Concept::Basis Basis, class Coefficients>
requires (Concept::LocalBasisTree<typename Basis::LocalView::Tree> &&
           PDELab::Execution::is_execution_policy_v<ExecutionPolicy> )
static std::map<std::string, double>
reduce(const ExecutionPolicy exec,
       const Basis& basis,
       const Coefficients& coefficients,
       auto time,
       const ParameterTree& config,
       std::shared_ptr<const FunctorFactory<Basis::EntitySet::GridView::dimension>>
         functor_factory = nullptr)
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

  using FEM = std::decay_t<decltype(first_finite_element(basis.localView().tree()))>;
  using Geometry = typename Basis::EntitySet::template Codim<0>::Entity::Geometry;

  struct ThreadLocalData
  {
    typename Basis::LocalView lbasis;
    PDELab::LocalContainerBuffer<Basis, const Coefficients> lcoeff;
    std::unique_ptr<LocalEquations<dim>> leqs;
    LocalBasisCache<typename FEM::Traits::LocalBasisType::Traits> fe_cache = {};
    std::optional<typename Geometry::JacobianInverse> geojacinv_opt = std::nullopt;
    std::vector<double> values = {};
    std::vector<fu2::unique_function<double() const>> evaluations = {};
    std::vector<fu2::unique_function<double(double, double) const>> reductions = {};
  };

  auto init_thread_data = [&]() -> ThreadLocalData {
    auto lbasis = basis.localView();
    auto leqs_ptr = LocalEquations<dim>::make(lbasis);
    leqs_ptr->time = time;
    const auto& leqs = *leqs_ptr;
    ThreadLocalData data{ std::move(lbasis), { basis, &coefficients }, std::move(leqs_ptr) };

    auto sz = config.getSubKeys().size();
    data.values.reserve(sz);
    data.evaluations.reserve(sz);
    data.reductions.reserve(sz);
    for (const auto& key : config.getSubKeys()) {
      const auto& sub_config = config.sub(key, true);
      data.values.emplace_back(sub_config.get("initial.value", 0.));
      auto& evaluation = data.evaluations.emplace_back();
      if (sub_config.hasKey("evaluation.expression"))
        evaluation =
          functor_factory->make_scalar(key + ".evaluation", sub_config.sub("evaluation"), leqs);
      if (not evaluation)
        evaluation = [] { return 0.; };
      auto& reduction = data.reductions.emplace_back(std::plus<>{});
      if (sub_config.hasKey("reduction.expression")) {
        auto [args, expr] =
          parser_context->parse_function_expression(sub_config.sub("reduction")["expression"]);
        if (args.size() == 3)
          spdlog::warn(
            "Reduction expression for '{0}' uses 3 arguments whereas 2 are now required.\nThis "
            "functionality has changed: the third argument became the contextual keyword "
            "'integration_factor' that may be used to express '{0}.evaluation.expression'.",
            key);
        if (args.size() != 2)
          throw format_exception(IOError{}, "Reduction arguments must be exactly 2");
        auto sub_parser_type = string2parser.at(
          sub_config.get("reduction.parser_type", std::string{ parser2string.at(parser_type) }));
        auto str_args = std::array{ std::string{ args[0] }, std::string{ args[1] } };
        reduction = parser_context->make_function(sub_parser_type, str_args, expr);
      }
    }
    return data;
  };

#if HAVE_TBB
  tbb::enumerable_thread_specific<ThreadLocalData> thread_data{ std::move(init_thread_data) };
#else
  struct ThreadData : public std::array<ThreadLocalData,1> {
    ThreadLocalData& local() {
      return this->front();
    }
  };
  ThreadData thread_data{init_thread_data()};
#endif

  if (config.getSubKeys().empty()) {
    return {};
  }

  spdlog::info("Evaluating reduce operators");
  PDELab::forEachEntity(exec, basis.entitySetPartition(), [&thread_data](const auto& cell) {
    auto& data = thread_data.local();
    auto geo = cell.geometry();

    data.lbasis.bind(cell);
    data.lcoeff.load(data.lbasis, std::false_type{});

    std::size_t const order = 4;
    const auto& quad_rule =
      QuadratureRules<typename Geometry::ctype, Geometry::coorddimension>::rule(geo.type(), order);

    for (std::size_t q = 0; q != quad_rule.size(); ++q) {
      const auto [position, weight] = quad_rule[q];
      if (not geo.affine() or not data.geojacinv_opt)
        data.geojacinv_opt.emplace(geo.jacobianInverse(position));
      const auto& geojacinv = *(data.geojacinv_opt);
      data.leqs->position = geo.global(position);
      data.leqs->integration_factor = weight * geo.integrationElement(position);

      // evaluate values at quadrature point
      forEachLeafNode(data.lbasis.tree(), [&](const auto& node) {
        if (node.size() == 0)
          return;
        auto& value = data.leqs->get_value(node);
        auto& gradient = data.leqs->get_gradient(node);
        value = 0.;
        gradient = 0.;
        data.fe_cache.bind(node.finiteElement(), quad_rule);
        const auto& phi = data.fe_cache.evaluateFunction(q);
        const auto& jacphi = data.fe_cache.evaluateJacobian(q);
        for (std::size_t dof = 0; dof != node.size(); ++dof) {
          value += data.lcoeff(node, dof) * phi[dof];
          gradient += data.lcoeff(node, dof) * (jacphi[dof] * geojacinv)[0];
        }
      });

      // evaluate all functions at the current quadrature point
      for (std::size_t i = 0; i != data.values.size(); ++i)
        data.values[i] = data.reductions[i](data.evaluations[i](), data.values[i]);
    }

    data.lbasis.unbind();
    data.leqs->clear();
  });

  if (basis.entitySet().grid().comm().size() > 1){
    throw format_exception(ParallelError{}, "Transoform and Reduce is not implemented in parallel");
  }

  // gather values from every thread_data
  std::vector<double> values;
  for (const auto& key : config.getSubKeys())
    values.emplace_back(config.sub(key, true).get("initial.value", 0.));
  for (std::size_t i = 0; i != values.size(); ++i)
    for (const auto& data : thread_data)
      values[i] = data.reductions[i](data.values[i], values[i]);

  // report results
  std::size_t max_key_chars = 0;
  std::size_t max_val_chars = 5;
  for (const auto& key : config.getSubKeys())
    max_key_chars = std::max(max_key_chars, key.size());

  std::map<std::string, double> key_value;

  spdlog::info("   ┌{0:─^{1}}┐", "", max_key_chars + max_val_chars + 13);
  std::string error_msg;
  auto values_it = values.begin();
  for (const auto& key : config.getSubKeys()) {
    auto& value = key_value[key] = *(values_it++);

    if (config.sub(key).hasKey("transformation.expression")) {
      auto [args, expr] = parser_context->parse_function_expression(config.sub(key)["transformation.expression"]);
      if (args.size() != 1)
        throw format_exception(IOError{}, "Warning function must have exactly 1 argument");
      auto sub_parser_type = string2parser.at(config.sub(key).get("transformation.parser_type", std::string{parser2string.at(parser_type)}));
      auto transformation = parser_context->make_function(sub_parser_type, std::array{std::string{args[0]}}, expr);
      value = transformation(value);
    }

    bool processed = false;

    if (config.sub(key).hasKey("error.expression")) {
      auto [args, expr] = parser_context->parse_function_expression(config.sub(key)["error.expression"]);
      if (args.size() != 1)
        throw format_exception(IOError{}, "Error function must have exactly 1 argument");
      auto sub_parser_type = string2parser.at(config.sub(key).get("error.parser_type", std::string{parser2string.at(parser_type)}));
      auto error = parser_context->make_function(sub_parser_type, std::array{std::string{args[0]}}, expr);
      if (FloatCmp::ne(error(value), 0.)) {
        processed = true;
        spdlog::error("  | {0:>{2}} := {1: .{3}e} |", DUNE_COPASI_FMT_STYLED_BOLD(key), value, max_key_chars, max_val_chars);
        if (not error_msg.empty())
          error_msg += '\n';
        error_msg += fmt::format("Reduction on the token '{}' raised an error because the "
                                 "expression '{}' with evaluates to false with '{} := {}'",
                                 key,
                                 expr,
                                 args[0],
                                 value);
      }
    }

    if (not processed and config.sub(key).hasKey("warn.expression")) {
      auto [args, expr] = parser_context->parse_function_expression(config.sub(key)["warn.expression"]);
      if (args.size() != 1)
        throw format_exception(IOError{}, "Warning function must have exactly 1 argument");
      auto sub_parser_type = string2parser.at(config.sub(key).get("warn.parser_type", std::string{parser2string.at(parser_type)}));
      auto warn = parser_context->make_function(sub_parser_type, std::array{std::string{args[0]}}, expr);
      if (FloatCmp::ne(warn(value), 0.)) {
        processed = true;
        spdlog::warn("| {0:>{2}} := {1: .{3}e} |", DUNE_COPASI_FMT_STYLED_BOLD(key), value, max_key_chars, max_val_chars);
      }
    }

    if (not processed and not config.sub(key).get("quiet", false)) {
      spdlog::info("   | {0:>{2}} := {1: .{3}e} |", DUNE_COPASI_FMT_STYLED_BOLD(key), value, max_key_chars, max_val_chars);
    }
  }
  spdlog::info("   └{0:─^{1}}┘", "", max_key_chars+max_val_chars+13);

  if (not error_msg.empty()) {
    throw format_exception(ReductionError{}, "{}", error_msg);
  }

  return key_value;
}

} // namespace Dune::Copasi::DiffusionReaction

#endif // DUNE_COPASI_MODEL_DIFFUSION_REACTION_REDUCE_HH
