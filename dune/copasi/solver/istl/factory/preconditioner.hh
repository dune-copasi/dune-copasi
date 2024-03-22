#ifndef DUNE_COPASI_SOLVER_ISTL_FACTORY_PRECONDITIONER_HH
#define DUNE_COPASI_SOLVER_ISTL_FACTORY_PRECONDITIONER_HH

#include <dune/copasi/common/parameterized_object.hh>
#include <dune/copasi/solver/istl/block_jacobi.hh>
#include <dune/copasi/solver/istl/concepts.hh>
#include <dune/copasi/solver/istl/util.hh>

#include <dune/istl/blocklevel.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/common/parameterizedobject.hh>

#include <fmt/ranges.h>

#define DUNE_COPASI_DEFAULT_PRECONDITIONER "SSOR"

namespace Dune::Copasi::ISTL {

template<Concept::LinearOperator O, class A>
using PreconditionerFactorySignature =
  std::shared_ptr<Preconditioner<typename O::domain_type, typename O::range_type>>(
    const std::shared_ptr<O>&,
    const ParameterTree&,
    const A&);
template<Concept::LinearOperator O, class A>
using PreconditionerRegistry = ParameterizedObjectFactoryWrapper<PreconditionerFactorySignature<O, A>>;

namespace Impl {

template<class Prec, class Alloc, bool matrix_free = false>
auto
defaultPreconditionerFactory(std::bool_constant<matrix_free> = {})
{
  using X = typename Prec::domain_type;
  using Y = typename Prec::range_type;
  return [&](const std::shared_ptr<LinearOperator<X, Y>>& op,
             const ParameterTree& config,
             const Alloc& alloc) -> std::shared_ptr<Preconditioner<X, Y>> {
    // make an allocator for this preconditioner
    using PrecAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<Prec>;
    PrecAlloc palloc(alloc);

    // construct preconditioner instance
    if constexpr (matrix_free) {
      return std::allocate_shared<Prec>(palloc, config);
    } else {
      using M = typename Prec::matrix_type; // cast operator to something the preconditioner expects
      auto aop = std::dynamic_pointer_cast<AssembledLinearOperator<M, X, Y>>(op);
      if (not aop)
        throw format_exception(InvalidStateException{}, "Linear operator does not hold a matrix!");
      return std::allocate_shared<Prec>(palloc, aop, config);
    }
  };
}

template<class X, class Y, class Alloc>
auto
inverseOperator2PreconditionerFactory()
{
  return [](const std::shared_ptr<LinearOperator<X, Y>>& op,
            const ParameterTree& config,
            const Alloc& alloc) -> std::shared_ptr<Preconditioner<X, Y>> {
    // wrapper: make sure inverse operator is stored somewhere
    struct Prec : public Dune::InverseOperator2Preconditioner<InverseOperator<X, Y>>
    {
      using InverseOperator2Preconditioner<InverseOperator<X, Y>>::InverseOperator2Preconditioner;
      std::shared_ptr<InverseOperator<X, Y>> inverse;
    };

    // make an allocator for this preconditioner
    using PrecAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<Prec>;
    PrecAlloc palloc(alloc);

    // make actual inverse operator and bind memory it to the wrapper
    // note that ADL tag is important because this function is not be yet defined!
    auto inverse = makeInverseOperator(op, config.sub("solver"), alloc, Impl::ADLTag{});
    auto prec = std::allocate_shared<Prec>(palloc, *inverse);
    prec->inverse = std::move(inverse);
    return prec;
  };
}

}

template<Concept::LinearOperator O, class Alloc>
const PreconditionerRegistry<O, Alloc>&
getPreconditionerRegistry()
{
  using X = typename O::domain_type;
  using Y = typename O::range_type;
  static PreconditionerRegistry<O, Alloc> registry;
  static std::once_flag flag;
  std::call_once(flag, [] {
    registry.define("Richardson", Impl::defaultPreconditionerFactory<Richardson<X, Y>, Alloc>(std::true_type()));
    registry.define("InverseOperator", Impl::inverseOperator2PreconditionerFactory<X, Y, Alloc>());
    // preconditioner that need a matrix
    if constexpr (Concept::AssembledLinearOperator<O>) {
      using M = typename O::matrix_type;
      registry.define("Jacobi", Impl::defaultPreconditionerFactory<SeqJac<M, X, Y, blockLevel<M>()>, Alloc>());
      registry.define("SSOR", Impl::defaultPreconditionerFactory<SeqSSOR<M, X, Y, blockLevel<M>()>, Alloc>());
      registry.define("SOR", Impl::defaultPreconditionerFactory<SeqSOR<M, X, Y, blockLevel<M>()>, Alloc>());
      registry.define("GaussSeidel", Impl::defaultPreconditionerFactory<SeqGS<M, X, Y, blockLevel<M>()>, Alloc>());
      // preconditioner that only have valid instantiation with matrices of depth level 1
      if constexpr (blockLevel<M>() == 1) {
        registry.define("DILU", Impl::defaultPreconditionerFactory<SeqDILU<M, X, Y>, Alloc>());
        registry.define("ILU", Impl::defaultPreconditionerFactory<SeqILU<M, X, Y>, Alloc>());
        registry.define("ILDL", Impl::defaultPreconditionerFactory<SeqILDL<M, X, Y>, Alloc>());
      }
      if constexpr (blockLevel<M>() >= 1) {
        registry.define("BlockJacobi", BlockJacobi<O, Alloc>::factory());
      }
    }
  });
  return registry;
}

template<Concept::LinearOperator Op, class Alloc = std::allocator<Op>>
std::shared_ptr<Preconditioner<typename Op::domain_type, typename Op::range_type>>
makePreconditioner(const std::shared_ptr<Op>& op,
                   const ParameterTree& config,
                   const Alloc& alloc = Alloc())
{
  assert(op);
  const auto& registry = getPreconditionerRegistry<Op, Alloc>();
  auto type_name = config.get("type", DUNE_COPASI_DEFAULT_PRECONDITIONER);
  if (config.get("verbosity", 1) > 0)
    spdlog::info("Creating preconditioner with type '{}'", type_name);
  if (registry.contains(type_name))
    return registry.create(type_name, op, config, alloc);
  else
    throw format_exception(IOError{},
                           "The key '{}' is not a known preconditioner type. Allowed types are {}",
                           type_name,
                           registry.keys());
}

} // Dune::Copasi::ISTL

#endif // DUNE_COPASI_SOLVER_ISTL_FACTORY_PRECONDITIONER_HH
