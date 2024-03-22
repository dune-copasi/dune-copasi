#ifndef DUNE_COPASI_SOLVER_ISTL_DENSE_INVERSE_HH
#define DUNE_COPASI_SOLVER_ISTL_DENSE_INVERSE_HH

#include <dune/istl/solver.hh>

namespace Dune::Copasi::ISTL {

template<class Inverse, class X, class Y>
struct DenseInverse final : public InverseOperator<X, Y>
{

  DenseInverse(Inverse&& inverse)
    : _inverse{ std::move(inverse) }
  {
  }

  void apply(X& x, Y& b, InverseOperatorResult& res) override
  {
    _inverse.mv(b,x);
    res.iterations = 1;
    res.converged = 1;
  }

  void apply(X& x, Y& b, [[maybe_unused]] double reduction, InverseOperatorResult& res) override
  {
    this->apply(x, b, res);
  }

  SolverCategory::Category category() const override
  {
    return SolverCategory::Category::sequential;
  }

  Inverse _inverse;
};

} // namespace Dune::Copasi::ISTL

#endif // DUNE_COPASI_SOLVER_ISTL_DENSE_INVERSE_HH
