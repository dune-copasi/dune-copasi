#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/copasi/common/parameter_base.hh>
#include <dune/copasi/common/parameterization.hh>

template<int I>
class Parameter : public Dune::Copasi::ParameterBase<Parameter<I>> 
{
public:

  template<class EtityType, class DomainType, class ParametrizationType>
  inline
  auto evaluate ( const EtityType& e,
                  const DomainType& x,
                  const ParametrizationType& p) const
  {
    return I;
  }
};

int main(int argc, char **argv)
{
  bool failed = false;
  try {
    int i;
    Dune::Copasi::Parametererization<Parameter<0>,Parameter<1>,Parameter<3>> p;
    assert(0==p.evaluate<Parameter<0>>(i,0,p));
    assert(1==p.evaluate<Parameter<1>>(i,1,p));
    assert(3==p.evaluate<Parameter<3>>(i,3,p));

    return failed;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return failed;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return failed;
  }
}