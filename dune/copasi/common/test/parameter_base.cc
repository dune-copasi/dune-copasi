#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/exceptions.hh>

#include <dune/copasi/common/parameter_base.hh>

class Parameter : public Dune::Copasi::ParameterBase<Parameter> 
{
public:

  template<class EtityType, class DomainType, class ParametrizationType>
  inline
  auto evaluate ( const EtityType& e,
                  const DomainType& x,
                  const ParametrizationType& p) const
  {
    return x;
  }
};

int main(int argc, char **argv)
{
  bool failed = false;

  try {
    int i;
    Parameter p;

    assert(1 == p.evaluate(i,1,i));
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