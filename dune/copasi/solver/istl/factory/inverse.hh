#ifndef DUNE_COPASI_SOLVER_ISTL_FACTORY_INVERSE_HH
#define DUNE_COPASI_SOLVER_ISTL_FACTORY_INVERSE_HH

#include <dune/copasi/solver/istl/factory/direct.hh>
#include <dune/copasi/solver/istl/factory/iterative.hh>
#include <dune/copasi/solver/istl/factory/preconditioner.hh>

#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solver.hh>

#ifdef DUNE_COPASI_DEFAULT_DIRECT_SOLVER
#define DUNE_COPASI_DEFAULT_LINEAR_INVERSE_OPERATOR DUNE_COPASI_DEFAULT_DIRECT_SOLVER
#else
#define DUNE_COPASI_DEFAULT_LINEAR_INVERSE_OPERATOR DUNE_COPASI_DEFAULT_ITERATIVE_SOLVER
#endif

namespace Dune::Copasi::ISTL {

namespace Impl {

template<Concept::LinearOperator Op, class Alloc = std::allocator<Op>>
std::shared_ptr<InverseOperator<typename Op::domain_type, typename Op::range_type>>
makeInverseOperator(const std::shared_ptr<Op>& op,
                    const ParameterTree& config,
                    const Alloc& alloc = Alloc(),
                    ADLTag = {})
{
  assert(op);
  auto type_name = config.get("type", DUNE_COPASI_DEFAULT_LINEAR_INVERSE_OPERATOR);
  if (config.get("verbosity", 1) > 0)
    spdlog::info("Creating linear solver with type '{}'", type_name);
  if (auto& registry = getIterativeSolverRegistry<Op, Alloc>(); registry.contains(type_name)) {
    auto prec = makePreconditioner(op, config.sub("preconditioner"));
    if constexpr (requires { op->getCommunication(); })
      throw format_exception(NotImplemented{}, "Parallel solvers have not been implemented!");
    using ScalarProd = ScalarProduct<typename Op::domain_type>;
    using ScalarProdAlloc =
      typename std::allocator_traits<Alloc>::template rebind_alloc<ScalarProd>;
    auto sp = std::allocate_shared<ScalarProd>(ScalarProdAlloc(alloc));
    return registry.create(type_name, op, sp, prec, config, alloc);
  } else if constexpr (Concept::AssembledLinearOperator<Op>) {
    if (auto& registry = getDirectSolverRegistry<Op, Alloc>(); registry.contains(type_name)) {
      return registry.create(type_name, op, config, alloc);
    }
  }

  // not found: diagnose error
  auto iterative_keys = getIterativeSolverRegistry<Op, Alloc>().keys();
  std::set<std::string> direct_keys;
  if constexpr (Concept::AssembledLinearOperator<Op>)
    direct_keys = getDirectSolverRegistry<Op, Alloc>().keys();

  throw format_exception(IOError{},
                         "The key '{}' is not a known inverse operator type.\nPossible iterative "
                         "solvers: {}\nPossible direct solvers: {}",
                         type_name,
                         iterative_keys,
                         direct_keys);
}

} // namespace Impl

template<Concept::LinearOperator Op, class Alloc = std::allocator<Op>>
std::shared_ptr<InverseOperator<typename Op::domain_type, typename Op::range_type>>
makeInverseOperator(const std::shared_ptr<Op>& op,
                    const ParameterTree& config,
                    const Alloc& alloc = Alloc())
{
  return Impl::makeInverseOperator(op, config, alloc);
}

} // Dune::Copasi::ISTL

#endif // DUNE_COPASI_SOLVER_ISTL_FACTORY_INVERSE_HH
