#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <dune/copasi/finite_element/dynamic_local_basis.hh>
#include <dune/copasi/finite_element/dynamic_local_coefficients.hh>
#include <dune/copasi/finite_element/dynamic_local_finite_element.hh>

#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lagrange.hh>

#include <dune/logging/logging.hh>

#include <dune/common/dynvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <cassert>
#include <set>

template<class Basis>
bool
test_power_local_basis(const Basis& basis, std::size_t power_size)
{
  using Domain = typename Basis::Traits::DomainType;
  using Range = typename Basis::Traits::RangeType;
  using Jacobian = typename Basis::Traits::JacobianType;

  Dune::Copasi::DynamicPowerLocalBasis<Basis> power_basis(basis, power_size);

  assert(basis.order() == power_basis.order());

  // evaluate always in the middle of the domain
  Domain in(0.5);
  std::vector<Range> out, power_out;
  std::vector<Jacobian> jout, power_jout;

  basis.evaluateFunction(in, out);
  power_basis.evaluateFunction(in, power_out);

  basis.evaluateJacobian(in, jout);
  power_basis.evaluateJacobian(in, power_jout);

  for (std::size_t i = 0; i < basis.size(); ++i)
    for (std::size_t j = 0; j < power_size; ++j) {
      // output of the blocked basis must be blocked by the size of the original
      // basis
      assert(out[i] == power_out[basis.size() * j + i]);
      assert(jout[i] == power_jout[basis.size() * j + i]);
    }

  return false;
}

template<class Coefficients>
bool
test_power_local_coefficients(const Coefficients& coefficients,
                              std::size_t power_size)
{
  Dune::Copasi::DynamicPowerLocalCoefficients<Coefficients> power_coefficients(
    power_size);

  assert(coefficients.size() * power_size == power_coefficients.size());

  std::set<Dune::LocalKey> unique_key;
  for (std::size_t i = 0; i < power_coefficients.size(); ++i) {
    auto key = power_coefficients.localKey(i);
    auto t = unique_key.insert(key);

    // ensure that the key is unique
    assert(t.second);
  }

  // check that local keys are ordered
  if (unique_key.size() == 0)
    return false;
  auto it = std::next(unique_key.begin());
  while (it != unique_key.end()) {
    const auto& a = *std::prev(it);
    const auto& b = *it;

    // ensure index keys are consecutive if codim and sub entity are the same
    if (a.subEntity() == b.subEntity() and a.codim() == b.codim())
      assert(a.index() + 1 == b.index());
    ++it;
  }

  return false;
}

// This class defines a local finite element function.
// It is determined by a local finite element and
// representing the local basis and a coefficient vector.
// This provides the evaluate method needed by the interpolate()
// method.
template<class FE>
class FEFunction
  : public Dune::Function<
      typename FE::Traits::LocalBasisType::Traits::DomainType,
      typename FE::Traits::LocalBasisType::Traits::RangeType>
{
  using Domain = typename FE::Traits::LocalBasisType::Traits::DomainType;
  using Range = typename FE::Traits::LocalBasisType::Traits::RangeType;

  const FE& fe;

public:
  using CT = typename FE::Traits::LocalBasisType::Traits::RangeFieldType;

  std::vector<CT> coeff;

  FEFunction(const FE& fe_)
    : fe(fe_)
  {
    resetCoefficients();
  }

  void resetCoefficients()
  {
    coeff.resize(fe.localBasis().size());
    for (std::size_t i = 0; i < coeff.size(); ++i)
      coeff[i] = 0;
  }

  void setRandom(double max)
  {
    coeff.resize(fe.localBasis().size());
    for (std::size_t i = 0; i < coeff.size(); ++i)
      coeff[i] = ((1.0 * std::rand()) / RAND_MAX - 0.5) * 2.0 * max;
  }

  void evaluate(const Domain& x, Range& y) const
  {
    std::vector<Range> yy;
    fe.localBasis().evaluateFunction(x, yy);

    y = 0.0;
    for (std::size_t i = 0; i < yy.size(); ++i)
      y.axpy(coeff[i], yy[i]);
  }
};

template<class F>
class DynamicFEFunction
  : public Dune::Function<typename F::DomainType,
                          Dune::DynamicVector<typename F::RangeType>>
{
  using Domain = typename F::DomainType;
  using Range = Dune::DynamicVector<typename F::RangeType>;

public:
  const F& _f;
  std::vector<double> _scales;
  DynamicFEFunction(const F& f, const std::vector<double>& scales)
    : _f(f)
    , _scales(scales)
  {}

  void evaluate(const Domain& x, Range& y) const
  {
    y.resize(_scales.size());
    typename F::RangeType y_base;
    _f.evaluate(x, y_base);
    for (std::size_t i = 0; i < _scales.size(); i++)
      y[i] = y_base * _scales[i];
  }
};

template<class F, class Interpolation>
bool
test_power_local_interpolation(const F& f,
                               const Interpolation& interpolation,
                               std::size_t power_size)
{
  Dune::Copasi::DynamicPowerLocalInterpolation<Interpolation>
    power_interpolation(power_size);

  bool failed = false;

  std::vector<double> coeff;
  interpolation.interpolate(f, coeff);

  std::vector<double> power_coeff;
  if (power_size == 1) {
    power_interpolation.interpolate(f, power_coeff);
    for (std::size_t j = 0; j < coeff.size(); j++)
      failed |= Dune::FloatCmp::ne(coeff[j], power_coeff[j]);
  }

  std::vector<double> scales(power_size);
  std::iota(scales.begin(), scales.end(), 0);
  DynamicFEFunction<F> dyn_f(f, scales);

  power_interpolation.interpolate(dyn_f, power_coeff);

  for (std::size_t i = 0; i < power_size; i++)
    for (std::size_t j = 0; j < coeff.size(); j++)
      failed |= Dune::FloatCmp::ne(coeff[j] * scales[i],
                                   power_coeff[coeff.size() * i + j]);

  return failed;
}

template<class FiniteElement>
bool
test_power_local_finite_element(const FiniteElement& finite_element,
                                std::size_t power_size)
{
  bool failed = false;

  // test constructor
  Dune::Copasi::DynamicPowerLocalFiniteElement<FiniteElement> power_fe(
    power_size);

  failed |= test_power_local_basis(finite_element.localBasis(), power_size);
  failed |= test_power_local_coefficients(finite_element.localCoefficients(),
                                          power_size);
  FEFunction<FiniteElement> f(finite_element);
  f.setRandom(1.);
  failed |= test_power_local_interpolation(
    f, finite_element.localInterpolation(), power_size);

  return failed;
}

int
main(int argc, char** argv)
{
  bool failed = false;

  try {
    // initialize mpi
    auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);
    auto comm = mpi_helper.getCollectiveCommunication();

    // initialize loggers
    Dune::Logging::Logging::init(comm);

    using RF = double;
    using DF = double;

    Dune::Pk2DLocalFiniteElement<DF, RF, 1> finite_element_1;
    for (int i = 0; i < 10; ++i)
      failed |= test_power_local_finite_element(finite_element_1, i);

    Dune::Pk2DLocalFiniteElement<DF, RF, 2> finite_element_2;
    for (int i = 0; i < 10; ++i)
      failed |= test_power_local_finite_element(finite_element_2, i);

    Dune::Pk2DLocalFiniteElement<DF, RF, 3> finite_element_3;
    for (int i = 0; i < 10; ++i)
      failed |= test_power_local_finite_element(finite_element_3, i);

    return failed;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
