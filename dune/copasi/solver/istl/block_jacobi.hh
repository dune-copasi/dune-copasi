#ifndef DUNE_COPASI_SOLVER_ISTL_BLOCK_JACOBI_HH
#define DUNE_COPASI_SOLVER_ISTL_BLOCK_JACOBI_HH

#include <dune/copasi/solver/istl/util.hh>
#include <dune/copasi/common/exceptions.hh>

#include <dune/pdelab/common/algebra.hh>

#include <dune/istl/preconditioner.hh>
#include <dune/istl/solver.hh>

#include <memory>

namespace Dune::Copasi::ISTL {

template<Concept::AssembledLinearOperator O, class Alloc = std::allocator<void>>
class BlockJacobi final : public Preconditioner<typename O::domain_type, typename O::range_type>
{
  using M = typename O::matrix_type;
  using X = typename O::domain_type;
  using Y = typename O::range_type;

public:
  BlockJacobi(const std::shared_ptr<const O>& op,
              const ParameterTree& config,
              const Alloc& alloc = Alloc())
    : _alloc(alloc)
    , _op(op)
    , _config{ config }
    , _parallel{ _config.get("parallel", false) }
    , _iterations{ _config.get("iterations", std::size_t{ 1 }) }
    , _weight{ _config.get("relaxation", 1.) }
  {
  }

  void pre(X& x, Y& b) override
  {
    if (_parallel)
      this->pre(PDELab::Execution::par, x, b);
    else
      this->pre(PDELab::Execution::seq, x, b);
  }

  void pre(auto policy, X& /* x */, Y& /* b */)
  {
    const ParameterTree& sub_solver = _config.sub("diagonal_solver");
    if constexpr (IsNumber<typename M::block_type>{}) {
      if (sub_solver.getSubKeys().size() != 0)
        spdlog::warn("'diagona_solver' contains sub-sections that will be ignored");
      if (sub_solver.get("type", "DenseInverse") != "DenseInverse")
        spdlog::warn("Block jacobi at this level is scalar, meaning that the diagonal solver can only be 'DenseInverse'");
    }

    const M& mat = _op->getmat();
    // store inverse of diagonal entries
    _diag_inv.resize(mat.M());
    PDELab::forEach(policy, mat, [&](const auto& row, auto i) {
      auto row_it = row.begin();
      for (; row_it != row.end() and row_it.index() != i; ++row_it) {
      };
      if (row_it == row.end())
        throw format_exception(MathError{}, "BlockJacobi can only be done on matrices with diagonal entries");
      if constexpr (IsNumber<typename M::block_type>{}) {
        _diag_inv[i] = 1. / (*row_it);
      } else {
        auto i_str = fmt::format("{}", i);
        const ParameterTree& sub_solver_i =
          (sub_solver.hasSub(i_str)) ? sub_solver.sub(i_str) : sub_solver;

        using MatrixAdapter = Dune::
          MatrixAdapter<typename M::block_type, typename X::value_type, typename Y::value_type>;
        using MatrixAdapterAlloc =
          typename std::allocator_traits<Alloc>::template rebind_alloc<MatrixAdapter>;
        MatrixAdapterAlloc maalloc(_alloc);
        auto diag_op = std::allocate_shared<MatrixAdapter>(maalloc, *row_it);
        _diag_inv[i] = makeInverseOperator(diag_op, sub_solver_i, _alloc, Impl::ADLTag{});
      }
    });
  }

  void apply(X& v, const Y& rhs) override
  {
    if (_parallel)
      this->apply(PDELab::Execution::par, v, rhs);
    else
      this->apply(PDELab::Execution::seq, v, rhs);
  }

  void apply(auto policy, X& v, const Y& rhs)
  {
    Y b = rhs;
    // we need a copy to avoid overriding old iterates, but is not needed if there is only one block
    std::optional<X> x_tmp;
    const M& mat = _op->getmat();

    using Xi = typename X::value_type;
    [[maybe_unused]] std::optional<Xi> vi;
    X& x = (v.size() > 1) ? x_tmp.emplace(v) : v;
    x = v;
    for (std::size_t it = 0; it != _iterations; ++it) {
      if constexpr (IsNumber<typename M::block_type>{}) {
        PDELab::forEach(policy, mat, [&b, &x, v, weight = _weight, diag_inv = _diag_inv.data()](const auto& row, auto i) mutable {
          // compute (b-Ax)
          for (auto col = row.begin(); col != row.end(); ++col)
            b[i] -= (*col) * v[col.index()];
          // solve v for Dv=(b-Ax); then x=x+v
          x[i] += weight * diag_inv[i] * b[i];
        });
      } else {
        PDELab::forEach(policy, mat, [&b, &x, v, vi, weight = _weight, diag_inv = _diag_inv.data()](const auto& row, auto i) mutable {
          // compute (b-Ax)
          InverseOperatorResult res;
          for (auto col = row.begin(); col != row.end(); ++col)
            col->mmv(v[col.index()], b[i]);
          // solve v for Dv=(b-Ax); then x=x+v
          auto& vir = vi ? *vi : vi.emplace();
          // copy 'x' in 'vi' is just to get the sizes right
          vir = x[i];
          diag_inv[i]->apply(vir, b[i], res);
          PDELab::axpy(x[i], weight, vir);
        });
      }
      if (v.size() > 1)
        v = x;
    }
  }

  void post(X&) override {}

  SolverCategory::Category category() const override { return SolverCategory::sequential; }

  static auto factory()
  {
    return [](const std::shared_ptr<LinearOperator<X, Y>>& op,
               const ParameterTree& config,
               const Alloc& alloc) -> std::shared_ptr<Preconditioner<X, Y>> {
      // make an allocator for this preconditioner
      using PrecAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<BlockJacobi>;
      PrecAlloc palloc(alloc);

      // cast operator to something the preconditioner expects
      auto aop = std::dynamic_pointer_cast<O>(op);
      if (not aop)
        throw format_exception(InvalidStateException{}, "Linear operator does not hold a matrix!");
      // construct preconditioner instance
      return std::allocate_shared<BlockJacobi>(palloc, aop, config);
    };
  }

private:
  using DiagInverse = std::conditional_t<
    IsNumber<typename M::block_type>::value,
    typename M::block_type,
    std::shared_ptr<InverseOperator<typename X::value_type, typename Y::value_type>>>;

  using DiagInverseAlloc =
    typename std::allocator_traits<Alloc>::template rebind_alloc<DiagInverse>;

  Alloc _alloc;
  std::shared_ptr<const O> _op;
  std::vector<DiagInverse, DiagInverseAlloc> _diag_inv;
  ParameterTree _config;
  bool _parallel;
  std::size_t _iterations;
  double _weight;
};

} // namespace Dune::Copasi::ISTL

#endif // DUNE_COPASI_SOLVER_ISTL_BLOCK_JACOBI_HH
