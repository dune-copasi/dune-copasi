#ifndef DUNE_COPASI_SOLVER_ISTL_FACTORY_DIRECT_HH
#define DUNE_COPASI_SOLVER_ISTL_FACTORY_DIRECT_HH

#include <dune/copasi/common/exceptions.hh>
#include <dune/copasi/solver/istl/concepts.hh>
#include <dune/copasi/solver/istl/dense_inverse.hh>

#if HAVE_SUITESPARSE_UMFPACK
#include <dune/copasi/solver/istl/umfpack.hh>
#define DUNE_COPASI_DEFAULT_DIRECT_SOLVER "UMFPack"
#endif

#include <dune/istl/solver.hh>
#include <dune/istl/superlu.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/common/parameterizedobject.hh>

#include <fmt/ranges.h>

namespace Dune::Copasi::ISTL {

// Iterative solver factory
template<class M, class X, class Y, class A>
using DirectSolverFactorySignature =
  std::shared_ptr<InverseOperator<X, Y>>(const std::shared_ptr<AssembledLinearOperator<M, X, Y>>&,
                                         const ParameterTree&,
                                         const A&);
template<class M, class X, class Y, class A>
using DirectSolverRegistry = ParameterizedObjectFactory<DirectSolverFactorySignature<M, X, Y, A>>;

namespace Impl {
template<class Op, class Alloc>
auto
denseInverseFactory()
{
  using M = typename Op::matrix_type;
  using X = typename Op::domain_type;
  using Y = typename Op::range_type;
  return [](const std::shared_ptr<AssembledLinearOperator<M, X, Y>>& op,
            const ParameterTree& config,
            const Alloc& alloc) -> std::shared_ptr<InverseOperator<X, Y>> {
    static_assert(blockLevel<M>() == 1, "Matrices for this solver can only be 1 level depth");
    std::shared_ptr<InverseOperator<X, Y>> inverse_op;
    using namespace Dune::Indices;
    const auto& mat = op->getmat();
    Hybrid::switchCases(
      Hybrid::integralRange(1_uc, 20_uc),
      mat.N(),
      [&](auto n) {
        if (mat.M() != n)
          throw format_exception(InvalidStateException{},
                                 "Matrix {}x{} is not square!",
                                 mat.N(),
                                 mat.M());
        // compute inverse out of the matrix
        using DenseMatrix = FieldMatrix<typename FieldTraits<M>::field_type, n, n>;
        DenseMatrix inverse;
        for (auto row_it = mat.begin(); row_it != mat.end(); ++row_it)
          for (auto col_it = row_it->begin(); col_it != row_it->end(); ++col_it)
            inverse[row_it.index()][col_it.index()] = *col_it;
        inverse.invert(config.get("pivoting", true));

        // make an allocator for this solver
        using Solver = DenseInverse<DenseMatrix, X, Y>;
        using SolverAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<Solver>;
        SolverAlloc ialloc(alloc);

        inverse_op = std::allocate_shared<Solver>(ialloc, std::move(inverse));
      },
      [] { throw format_exception(NotImplemented{}, "Dense inverses are only implemented up to 19 entries"); });
    return inverse_op;
  };
}

#if HAVE_SUPERLU
template<class Op, class Alloc>
auto
superLUFactory()
{
  using M = typename Op::matrix_type;
  using X = typename Op::domain_type;
  using Y = typename Op::range_type;
  return [](const std::shared_ptr<AssembledLinearOperator<M, X, Y>>& op,
            const ParameterTree& config,
            const Alloc& alloc) -> std::shared_ptr<InverseOperator<X, Y>> {
    static_assert(blockLevel<M>() == 1, "Matrices for this solver can only be 1 level depth");
    std::shared_ptr<InverseOperator<X, Y>> inverse_op;
    auto verbosity = config.get("verbosity", 0);
    auto reuse_vector = config.get("reuse_vector", 0);
    // make an allocator for this solver
    using Solver = Dune::SuperLU<M>;
    using SolverAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<Solver>;
    SolverAlloc ialloc(alloc);
    // TODO(sospinar): check if forwarding allocator to superlu is possible
    inverse_op = std::allocate_shared<Solver>(ialloc, op->getmat(), verbosity, reuse_vector);
    return inverse_op;
  };
}
#endif // HAVE_SUPERLU

} // namespace Impl

template<Concept::AssembledLinearOperator O, class Alloc>
const DirectSolverRegistry<typename O::matrix_type,
                           typename O::domain_type,
                           typename O::range_type,
                           Alloc>&
getDirectSolverRegistry()
{
  using M = typename O::matrix_type;
  using X = typename O::domain_type;
  using Y = typename O::range_type;
  // initialize a unique registry instance
  static DirectSolverRegistry<M, X, Y, Alloc> registry;
  // register possible solver factories
  static std::once_flag flag;
  std::call_once(flag, [] {
    if constexpr (blockLevel<M>() == 1)
      registry.define("DenseInverse", Impl::denseInverseFactory<O, Alloc>());
#if HAVE_SUITESPARSE_UMFPACK
    registry.define("UMFPack", UMFPackWapper<O, Alloc>::factory());
#endif // HAVE_SUITESPARSE_UMFPACK
#if HAVE_SUPERLU
    if constexpr (blockLevel<M>() == 1)
      registry.define("SuperLU", Impl::superLUFactory<O, Alloc>());
#endif // HAVE_SUPERLU
  });
  return registry;
}

template<Concept::AssembledLinearOperator Op, class Alloc = std::allocator<Op>>
std::shared_ptr<InverseOperator<typename Op::domain_type, typename Op::range_type>>
makeDirectSolver(const std::shared_ptr<Op>& op,
                 const ParameterTree& config,
                 const Alloc& alloc = Alloc())
{
  assert(op);
  auto& registry = getDirectSolverRegistry<Op, Alloc>();
#ifdef DUNE_COPASI_DEFAULT_DIRECT_SOLVER
  auto type_name = config.get("type", DUNE_COPASI_DEFAULT_DIRECT_SOLVER);
#else
  auto type_name = config.get("type", "None");
#endif
  if (registry.contains(type_name))
    return registry.create(type_name, op, config, alloc);
  else
    throw format_exception(IOError{},
                           "The key '{}' is not a known direct solver type. Allowed types are {}",
                           type_name,
                           /*registry.keys()*/ "<unknown>");
}

} // Dune::Copasi::ISTL

#endif // DUNE_COPASI_SOLVER_ISTL_FACTORY_DIRECT_HH
