#ifndef DUNE_COPASI_LOCAL_FINITE_ELEMENT_PK_HH
#define DUNE_COPASI_LOCAL_FINITE_ELEMENT_PK_HH

#include <dune/copasi/common/factory.hh>

#include <dune/localfunctions/lagrange/pk.hh>

#include <dune/geometry/type.hh>

namespace Dune::Copasi {

template<class D, class R, int d, int k>
struct Factory<Dune::PkLocalFiniteElement<D,R,d,k>>
{
  template<class Ctx>
  static auto create(Ctx&& ctx)
  {
    return std::make_unique<Dune::PkLocalFiniteElement<D,R,d,k>>();
  }
};

} // namespace Dune::Copasi

#endif // DUNE_COPASI_LOCAL_FINITE_ELEMENT_PK_HH