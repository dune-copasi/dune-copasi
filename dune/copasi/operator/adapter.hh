#ifndef DUNE_COPASI_OPERATOR_ADAPTER_HH
#define DUNE_COPASI_OPERATOR_ADAPTER_HH

#include <dune/copasi/operator/operator.hh>

#include <memory>
#include <type_traits>
#include <functional>

namespace Dune::Copasi::inline Experimental {

/**
 * @brief Adaptes an operator @f$f@f$ to work with generic functors
 *
 * @tparam Domain
 * @tparam Range
 */
template<std::copy_constructible Domain, class Range>
class OperatorAdapter : public Operator<Domain,Range>
{
public:

  template<class F>
  requires (std::is_invocable_r_v<PDELab::ErrorCondition,F,Operator<Domain,Range>&,Range&,Domain&>)
  OperatorAdapter(F&& apply_op) {
    setApply(std::forward<F>(apply_op), PriorityTag<1>{});
  }

  virtual PDELab::ErrorCondition apply(const Domain& domain, Range& range) override {
    return _capply_op(*this, domain, range);
  }

  virtual PDELab::ErrorCondition apply(Domain& domain, Range& range) override {
    return _apply_op(*this, domain, range);
  }

private:

  template<class F>
  requires (std::copy_constructible<Range> && std::is_invocable_r_v<PDELab::ErrorCondition,F,Operator<Domain,Range>&,Domain&,Range&>)
  void setApply(F&& apply_op, PriorityTag<1>) {
    _apply_op = std::forward<F>(apply_op);
    _capply_op = [&](Operator<Domain,Range>& self, const Domain& domain, Range& range) mutable {
      Range domain_copy = domain;
      return _apply_op(self, domain_copy, range);
    };
  }

  template<class F>
  requires std::is_invocable_r_v<PDELab::ErrorCondition,F,Operator<Domain,Range>&,const Domain&,Range&>
  void setApply(F&& apply_op, PriorityTag<0>) {
    _apply_op = apply_op;
    _capply_op = std::forward<F>(apply_op);
  }

  std::function<PDELab::ErrorCondition(Operator<Domain,Range>&,Domain&,Range&)> _apply_op;
  std::function<PDELab::ErrorCondition(Operator<Domain,Range>&,const Domain&,Range&)> _capply_op;
};

} // namespace Dune::Copasi::inline Experimental

#endif // DUNE_COPASI_OPERATOR_ADAPTER_HH
