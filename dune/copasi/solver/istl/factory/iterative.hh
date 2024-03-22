#ifndef DUNE_COPASI_SOLVER_ISTL_FACTORY_ITERATIVE_HH
#define DUNE_COPASI_SOLVER_ISTL_FACTORY_ITERATIVE_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/common/parameterized_object.hh>
#include <dune/copasi/solver/istl/concepts.hh>
#include <dune/copasi/solver/istl/util.hh>

#include <dune/istl/blocklevel.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/solver.hh>
#include <dune/istl/solvers.hh>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>

#include <fmt/ranges.h>

#define DUNE_COPASI_DEFAULT_ITERATIVE_SOLVER "BiCGSTAB"

namespace Dune::Copasi::ISTL {

// Iterative solver factory
template<class X, class Y, class A>
using IterativeSolverFactorySignature =
  std::shared_ptr<IterativeSolver<X, Y>>(const std::shared_ptr<LinearOperator<X, Y>>&,
                                         const std::shared_ptr<ScalarProduct<X>>&,
                                         const std::shared_ptr<Preconditioner<X, Y>>&,
                                         const ParameterTree&,
                                         const A&);
template<class X, class Y, class A>
using IterativeSolverRegistry =
  ParameterizedObjectFactoryWrapper<IterativeSolverFactorySignature<X, Y, A>>;

namespace Impl {
template<class Solver, class Alloc>
auto
defaultIterativeSolverFactory()
{
  using X = typename Solver::domain_type;
  using Y = typename Solver::range_type;
  return [](const std::shared_ptr<LinearOperator<X, Y>>& op,
            const std::shared_ptr<ScalarProduct<X>>& sp,
            const std::shared_ptr<Preconditioner<X, Y>>& prec,
            const ParameterTree& config,
            const Alloc& alloc) -> std::shared_ptr<IterativeSolver<X, Y>> {
    // make an allocator for this solver
    using SolverAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<Solver>;
    SolverAlloc salloc(alloc);

    // read config options
    const auto& conv_cond = config.sub("convergence_condition");
    auto reduction = conv_cond.get("relative_tolerance", 1e-8);
    auto [min_it, max_it] = conv_cond.get("iteration_range", std::array<std::size_t, 2>{ 1, 500 });
    if (min_it != 1)
      spdlog::warn(
        "Minimum number of iterations '{}' cannot be set for solver '{}. It has been set to '1'.'",
        min_it,
        className<Solver>());
    auto verbosity = config.get("verbosity", 0);

    // construct iterative solver instance (and handle special constructor cases :-|)
    if constexpr (std::derived_from<Solver, RestartedGMResSolver<X, Y>>) {
      auto restart = config.get("restart", 40);
      return std::allocate_shared<Solver>(salloc, op, sp, prec, reduction, restart, max_it, verbosity);
    } else if constexpr (std::derived_from<Solver, GeneralizedPCGSolver<X>> or std::derived_from<Solver, RestartedFCGSolver<X>>) {
      auto restart = config.get("restart", 10);
      return std::allocate_shared<Solver>(salloc, op, sp, prec, reduction, max_it, verbosity, restart);
    } else {
      return std::allocate_shared<Solver>(salloc, op, sp, prec, reduction, max_it, verbosity);
    }
  };
}
}

template<Concept::LinearOperator O, class Alloc>
const IterativeSolverRegistry<typename O::domain_type, typename O::range_type, Alloc>&
getIterativeSolverRegistry()
{
  using X = typename O::domain_type;
  using Y = typename O::range_type;
  // initialize a unique registry instance
  static IterativeSolverRegistry<X, Y, Alloc> registry;
  // register possible solver factories
  static std::once_flag flag;
  std::call_once(flag, [] {
    registry.define("Loop",
                    Impl::defaultIterativeSolverFactory<LoopSolver<X>, Alloc>());
    registry.define("Gradient",
                    Impl::defaultIterativeSolverFactory<GradientSolver<X>, Alloc>());
    registry.define("CG",
                    Impl::defaultIterativeSolverFactory<CGSolver<X>, Alloc>());
    registry.define("BiCGSTAB",
                    Impl::defaultIterativeSolverFactory<BiCGSTABSolver<X>, Alloc>());
    registry.define("MINRES",
                    Impl::defaultIterativeSolverFactory<MINRESSolver<X>, Alloc>());
    registry.define("RestartedGMRes",
                    Impl::defaultIterativeSolverFactory<RestartedGMResSolver<X, Y>, Alloc>());
    registry.define("RestartedFlexibleGMRes",
                    Impl::defaultIterativeSolverFactory<RestartedFlexibleGMResSolver<X, Y>, Alloc>());
    registry.define("GeneralizedPCG",
                    Impl::defaultIterativeSolverFactory<GeneralizedPCGSolver<X>, Alloc>());
    registry.define("RestartedFCG",
                    Impl::defaultIterativeSolverFactory<RestartedFCGSolver<X>, Alloc>());
    registry.define("CompleteFCG",
                    Impl::defaultIterativeSolverFactory<CompleteFCGSolver<X>, Alloc>());
  });
  return registry;
}

template<Concept::LinearOperator Op, class Alloc = std::allocator<Op>>
std::shared_ptr<IterativeSolver<typename Op::domain_type, typename Op::range_type>>
makeIterativeSolver(
  const std::shared_ptr<Op>& op,
  const std::shared_ptr<ScalarProduct<typename Op::domain_type>>& sp,
  const std::shared_ptr<Preconditioner<typename Op::domain_type, typename Op::range_type>>& prec,
  const ParameterTree& config,
  const Alloc& alloc = Alloc())
{
  assert(op);
  assert(sp);
  assert(prec);
  auto& registry = getIterativeSolverRegistry<Op, Alloc>();
  auto type_name = config.get("type", DUNE_COPASI_DEFAULT_ITERATIVE_SOLVER);
  if (config.get("verbosity", 1) > 0)
    spdlog::info("Creating iterative solver with type '{}'", type_name);
  if (registry.contains(type_name))
    return registry.create(type_name, op, sp, prec, config, alloc);
  else
    throw format_exception(
      IOError{},
      "The key '{}' is not a known iterative solver type. Allowed types are {}",
      type_name,
      registry.keys());
}

} // Dune::Copasi::ISTL

#endif // DUNE_COPASI_SOLVER_ISTL_FACTORY_ITERATIVE_HH
