#ifndef DUNE_COPASI_GRID_FUNCTION_COMPARE_HH
#define DUNE_COPASI_GRID_FUNCTION_COMPARE_HH

#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/function/minus.hh>

#include <dune/grid/common/rangegenerators.hh>

#include <dune/common/parametertree.hh>
#include <dune/common/float_cmp.hh>

#include <cmath>

template<class GF, class T, class TernaryOp>
T
grid_function_reduce(const GF& gf, unsigned int int_order, T init, TernaryOp ternary_op)
{
  using Range = typename GF::Traits::RangeType;

  const auto& gv = gf.getGridView();

  // Loop over the grid
  for (const auto& cell : Dune::elements(gv,Dune::Partitions::interior))
  {
    // Get geometry
    auto geo = cell.geometry();

    // Integrate on the cell
    for (const auto& point : quadratureRule(geo, int_order))
    {
      // Evaluate function
      Range y;
      const auto& x = point.position();
      gf.evaluate(cell,x,y);

      auto factor = point.weight() * geo.integrationElement(x);

      // reduce operator
      init = ternary_op(init,y,factor);
    }
  }

  return init;
}


template<class GF_A, class GF_B>
void grid_function_compare(const Dune::ParameterTree& param, GF_A& gf_a, GF_B& gf_b)
{
  if (not param.hasKey("l1_error") and
      not param.hasKey("l2_error") and
      not param.hasKey("linf_error"))
    DUNE_THROW(Dune::IOError,"At least one norm has to be compared.");

  Dune::PDELab::MinusGridFunctionAdapter<GF_A,GF_B> gf_diff(gf_a,gf_b);

  using RangeField = typename GF_A::Traits::RangeFieldType;

  // reduce operation for l-1 norm
  auto accumulate_op = [](auto acc, auto y, auto factor){return acc+std::abs(y)*factor;};

  // reduce operation for l-2 norm
  auto two_norm2_op = [](auto two_norm2, auto y, auto factor){return two_norm2+y.two_norm2()*factor;};

  // reduce operation for l-inf norm
  auto max_diff_op = [](auto max_y, auto y, auto factor){return std::max(max_y,std::abs(y));};

  // reduce operation for l-1, l-2, and l-inf norms
  auto reduce_op = [&](auto val, auto y, auto factor)
  {
    val[0] = accumulate_op(val[0],y,factor);
    val[1] = two_norm2_op(val[1],y,factor);
    val[2] = max_diff_op(val[2],y,factor);
    return val;
  };

  Dune::FieldVector<RangeField,3>  zero{0.,0.,0.};
  auto norms = grid_function_reduce(gf_diff,5,zero,reduce_op);
  norms[1] = std::sqrt(norms[1]);

  const auto& l1_norm = norms[0];
  const auto& l2_norm = norms[1];
  const auto& linf_norm = norms[2];

  if (param.hasKey("l1_error")){
    const auto max_l1_norm = param.template get<RangeField>("l1_error");
    if (Dune::FloatCmp::gt(l1_norm,max_l1_norm))
      DUNE_THROW(Dune::MathError, "l-1 error is " << l1_norm << " whle the maximum allowed is " << max_l1_norm);
  }

  if (param.hasKey("l2_error")){
    const auto max_l2_norm = param.template get<RangeField>("l2_error");
    if (Dune::FloatCmp::gt(l2_norm,max_l2_norm))
      DUNE_THROW(Dune::MathError, "l-2 error is " << l2_norm << " whle the maximum allowed is " << max_l2_norm);
  }

  if (param.hasKey("linf_error")){
    const auto max_linf_norm = param.template get<RangeField>("linf_error");
    if (Dune::FloatCmp::gt(linf_norm,max_linf_norm))
      DUNE_THROW(Dune::MathError, "l-inf error is " << linf_norm << " whle the maximum allowed is " << max_linf_norm);
  }
}

#endif // DUNE_COPASI_GRID_FUNCTION_COMPARE_HH