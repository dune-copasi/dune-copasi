#include <dune/copasi/finite_element/dynamic_power.hh>
#include <dune/copasi/finite_element/dynamic_power/local_basis.hh>
#include <dune/copasi/finite_element/dynamic_power/local_coefficients.hh>

#include <dune/localfunctions/common/localkey.hh>
#include <dune/localfunctions/lagrange.hh>

#include <dune/common/dynvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <gtest/gtest.h>

template<class Basis>
void
test_power_local_basis(const Basis& basis, std::size_t power_size)
{
  using Domain = typename Basis::Traits::DomainType;
  using Range = typename Basis::Traits::RangeType;
  using Jacobian = typename Basis::Traits::JacobianType;

  Dune::Copasi::DynamicPowerLocalBasis<Basis> power_basis(basis, power_size);

  EXPECT_EQ(basis.order(), power_basis.order());

  // evaluate always in the middle of the domain
  Domain in(0.5);
  std::vector<Range> out, power_out;
  std::vector<Jacobian> jout, power_jout;

  basis.evaluateFunction(in, out);
  power_basis.evaluateFunction(in, power_out);

  basis.evaluateJacobian(in, jout);
  power_basis.evaluateJacobian(in, power_jout);

  EXPECT_EQ(power_size * out.size(), power_out.size());
  EXPECT_EQ(power_size * jout.size(), power_jout.size());

  for (std::size_t i = 0; i < basis.size(); ++i) {
    for (std::size_t j = 0; j < power_size; ++j) {
      // output of the blocked basis must be blocked by the size of the original
      // basis
      EXPECT_EQ(out[i], power_out[basis.size() * j + i]);
      EXPECT_EQ(jout[i], power_jout[basis.size() * j + i]);
    }
  }
}

template<class Coefficients>
void
test_power_local_coefficients(const Coefficients& coefficients, std::size_t power_size)
{
  Dune::Copasi::DynamicPowerLocalCoefficients<Coefficients> power_coefficients(power_size);

  EXPECT_EQ(coefficients.size() * power_size, power_coefficients.size());

  std::set<Dune::LocalKey> unique_key;
  for (std::size_t i = 0; i < power_coefficients.size(); ++i) {
    // ensure that the key is unique
    auto key = power_coefficients.localKey(i);
    auto t = unique_key.insert(key);
    EXPECT_TRUE(t.second) << "Key is not unique";
  }

  // check that local keys are ordered
  if (unique_key.size() == 0)
    return;
  auto it = std::next(unique_key.begin());
  while (it != unique_key.end()) {
    const auto& a = *std::prev(it);
    const auto& b = *it;

    // ensure index keys are consecutive if codim and sub entity are the same
    if (a.subEntity() == b.subEntity() and a.codim() == b.codim())
      EXPECT_EQ(a.index() + 1, b.index());
    ++it;
  }
}

template<class C = double>
auto
get_random_coeff(const std::size_t& size, const double& max)
{
  std::vector<C> coeff;
  coeff.resize(size);
  for (std::size_t i = 0; i < coeff.size(); ++i)
    coeff[i] = ((1.0 * std::rand()) / RAND_MAX - 0.5) * 2.0 * max;
  return coeff;
}

template<class DomainType, class F, class Interpolation>
void
test_power_local_interpolation(const F& f,
                               const Interpolation& interpolation,
                               std::size_t power_size)
{
  Dune::Copasi::DynamicPowerLocalInterpolation<Interpolation, DomainType> power_interpolation(
    power_size);

  std::vector<double> coeff;
  interpolation.interpolate(f, coeff);

  std::vector<double> power_coeff;
  if (power_size == 1) {
    power_interpolation.interpolate(f, power_coeff);
    for (std::size_t j = 0; j < coeff.size(); j++)
      EXPECT_FLOAT_EQ(coeff[j], power_coeff[j]);
  }

  std::vector<double> scales(power_size);
  std::iota(scales.begin(), scales.end(), 0);

  auto dyn_f = [&](const auto& x) {
    auto y_base = f(x);
    using Range = Dune::DynamicVector<decltype(y_base)>;
    Range y(scales.size());
    for (std::size_t i = 0; i < scales.size(); i++)
      y[i] = y_base * scales[i];
    return y;
  };

  power_interpolation.interpolate(dyn_f, power_coeff);

  for (std::size_t i = 0; i < power_size; i++)
    for (std::size_t j = 0; j < coeff.size(); j++)
      EXPECT_FLOAT_EQ(coeff[j] * scales[i], power_coeff[coeff.size() * i + j]);
}

template<class FiniteElement>
void
test_power_local_finite_element(const FiniteElement& finite_element, std::size_t power_size)
{
  // test constructor
  using FE = Dune::Copasi::DynamicPowerLocalFiniteElement<FiniteElement>;
  auto power_fe = FE{ power_size };

  test_power_local_basis(finite_element.localBasis(), power_size);
  test_power_local_coefficients(finite_element.localCoefficients(), power_size);

  auto coeff = get_random_coeff(finite_element.localBasis().size(), 1.);
  auto f = [&](const auto& x) {
    using Range = typename FE::Traits::LocalBasisType::Traits::RangeType;
    std::vector<Range> yy;
    finite_element.localBasis().evaluateFunction(x, yy);

    Range y = 0.0;
    for (std::size_t i = 0; i < yy.size(); ++i)
      y.axpy(coeff[i], yy[i]);
    return y;
  };

  using Domain = typename FE::Traits::LocalBasisType::Traits::DomainType;
  test_power_local_interpolation<Domain>(f, finite_element.localInterpolation(), power_size);
}

template<std::size_t dim, std::size_t k>
class Pk2DLocalFiniteElementFixture : public ::testing::Test
{
protected:
  using FE = Dune::LagrangeSimplexLocalFiniteElement<double, double, dim, k>;

  void SetUp() override { _finite_element = std::make_unique<FE>(); }

  std::unique_ptr<FE> _finite_element;
};

// declare a dummy class that forwards it's template argument as the fixture to test
template<typename Fixture>
class ForwardFixture : public Fixture
{};

// declare the templated fixture
TYPED_TEST_SUITE_P(ForwardFixture);

// declare test to make on each space
TYPED_TEST_P(ForwardFixture, TestPowerLocalFiniteElement)
{
  for (int i = 0; i < 10; ++i)
    test_power_local_finite_element(*this->_finite_element, i);
}

// register the test
REGISTER_TYPED_TEST_SUITE_P(ForwardFixture, TestPowerLocalFiniteElement);

using LocalFiniteElements = ::testing::Types<Pk2DLocalFiniteElementFixture<1, 1>,
                                             Pk2DLocalFiniteElementFixture<1, 2>,
                                             Pk2DLocalFiniteElementFixture<1, 3>,
                                             Pk2DLocalFiniteElementFixture<2, 1>,
                                             Pk2DLocalFiniteElementFixture<2, 2>,
                                             Pk2DLocalFiniteElementFixture<2, 3>,
                                             Pk2DLocalFiniteElementFixture<3, 1>,
                                             Pk2DLocalFiniteElementFixture<3, 2>,
                                             Pk2DLocalFiniteElementFixture<3, 3>>;

// instantiate the test for each type
INSTANTIATE_TYPED_TEST_SUITE_P(PowerLocalFiniteElement, ForwardFixture, LocalFiniteElements);
