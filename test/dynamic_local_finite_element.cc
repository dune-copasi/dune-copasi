#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <dune/copasi/dynamic_local_coefficients.hh>
#include <dune/copasi/dynamic_local_basis.hh>

#include <dune/localfunctions/lagrange/pk.hh>
#include <dune/localfunctions/common/localkey.hh>

#include <dune/logging/logging.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <iostream>
#include <set>
#include <cassert>

template<class Basis>
bool test_power_local_basis(const Basis& basis, std::size_t power_size)
{
  using Domain = typename Basis::Traits::DomainType;
  using Range = typename Basis::Traits::RangeType;
  using Jacobian = typename Basis::Traits::JacobianType;

  Dune::Copasi::DynamicPowerLocalBasis<Basis> power_basis(basis,power_size);

  assert(basis.order() == power_basis.order());

  // evaluate always in the middle of the domain
  Domain in(0.5);
  std::vector<Range> out, power_out;
  std::vector<Jacobian> jout, power_jout;

  basis.evaluateFunction(in,out);
  power_basis.evaluateFunction(in,power_out);

  basis.evaluateJacobian(in,jout);
  power_basis.evaluateJacobian(in,power_jout);

  for (std::size_t i = 0; i < basis.size(); ++i)
    for (std::size_t j = 0; j < power_size; ++j) 
    {
      // output of the blocked basis must be blocked by the size of the original basis
      assert(out[i] == power_out[basis.size()*j+i]);
      assert(jout[i] == power_jout[basis.size()*j+i]);
    }

  return false;
}

template<class Coefficients>
bool test_power_local_coefficients(const Coefficients& coefficients, std::size_t power_size)
{
  Dune::Copasi::DynamicPowerLocalCoefficients<Coefficients> power_coefficients(power_size);
  
  assert( coefficients.size()*power_size == power_coefficients.size() );

  std::set<Dune::LocalKey> unique_key;
  for (int i = 0; i < power_coefficients.size(); ++i)
  {
    auto key = power_coefficients.localKey(i);
    auto t = unique_key.insert(key);

    // ensure that the key is unique
    assert(t.second);
  }

  // check that local keys are order
  auto it = std::next(unique_key.begin());
  while (it != unique_key.end())
  {
    const auto& a = *std::prev(it);
    const auto& b = *it;

    // ensure index keys are consecutive if codim and sub entity are the same
    if (a.subEntity()==b.subEntity() and a.codim()==b.codim())
      assert( a.index()+1 == b.index());
    ++it;
  }

  return false;
}

int main(int argc, char** argv)
{
  bool failed = false;

  try{
    // initialize mpi
    auto& mpi_helper = Dune::MPIHelper::instance(argc, argv);
    auto comm = mpi_helper.getCollectiveCommunication();

    // initialize loggers
    Dune::Logging::Logging::init(comm);

    using RF = double;
    using DF = double;

    constexpr uint k = 3;
    Dune::Pk2DLocalBasis<DF,RF,k> basis;
    Dune::Pk2DLocalCoefficients<k> coefficients;

    for (int i = 1; i < 10; ++i) 
    {
      failed |= test_power_local_basis<Dune::Pk2DLocalBasis<DF,RF,k>>(basis,3);
      failed |= test_power_local_coefficients<Dune::Pk2DLocalCoefficients<k>>(coefficients,i);
    }
    
    return failed;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
